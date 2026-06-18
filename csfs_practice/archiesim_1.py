#!/u/home/h/haroldzw/miniforge3/envs/python311/bin/python
import demes
import numpy as np
import msprime
import argparse
import json
def simulate_under_prufer(nsite,deme,sample,seed=None,nrep=1,gen_time=29,rec_rate=1.e-8,mu=1.25e-8):
    print(seed)
    graph=demes.load(deme)
    demography=msprime.Demography.from_demes(graph)
    with open(sample, "r") as json_file:
        # Use json.load() to parse the JSON data from the file
        samples = json.load(json_file)
    ancestry_reps = msprime.sim_ancestry(
        recombination_rate=rec_rate ,
        sequence_length=nsite,
        num_replicates = nrep,
        samples=samples,
        demography = demography,
        record_migrations=True,  # Needed for tracking segments.
        random_seed=seed
    )
    # graph=demography.to_demes()
    # demes.dump(graph,"test.yaml")
    for ts in ancestry_reps:
        mutated_ts = msprime.sim_mutations(ts, random_seed=seed,rate=mu)
        yield mutated_ts
"""
def run_simulation(sequence_length, random_seed=None):
    time_units = 1000 / 25  # Conversion factor for kya to generations
    demography = msprime.Demography()
    # The same size for all populations; highly unrealistic!
    Ne = 10**4
    demography.add_population(name="Africa", initial_size=Ne)
    demography.add_population(name="Eurasia", initial_size=Ne)
    demography.add_population(name="Neanderthal", initial_size=Ne)

    # 2% introgression 50 kya
    demography.add_mass_migration(
        time=50 * time_units, source='Eurasia', dest='Neanderthal', proportion=0.02)
    # Eurasian 'merges' backwards in time into Africa population, 70 kya
    demography.add_mass_migration(
        time=70 * time_units, source='Eurasia', dest='Africa', proportion=1)
    # Neanderthal 'merges' backwards in time into African population, 300 kya
    demography.add_mass_migration(
        time=300 * time_units, source='Neanderthal', dest='Africa', proportion=1)

    ts = msprime.sim_ancestry(
        recombination_rate=1e-8,
        sequence_length=sequence_length,  
        samples=[
            msprime.SampleSet(100, ploidy=1, population='Africa'),
            msprime.SampleSet(100, ploidy=1, population='Eurasia'),
            # Neanderthal sample taken 30 kya
            msprime.SampleSet(1, ploidy=1, time=30 * time_units, population='Neanderthal'),
        ],
        demography = demography,
        record_migrations=True,  # Needed for tracking segments.
        random_seed=random_seed,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=25)
    print(f"Simulation of {sequence_length/10**6}Mb run, using record_migrations=True")
    print(
        "NB: time diff from Neanderthal split to admixture event is",
        f"{300 * time_units - 50 * time_units:.0f} gens",
        f"({(300 * time_units - 50 * time_units) / 2 / Ne} coalescence units)"
    )
    return ts
"""
def findpop(tag,ts):
    for pop in ts.populations():
        if pop.metadata['name']==tag:
            return pop.id
    raise ValueError("the tag does not exist in the population")
def get_true_introgressed_segments(tag,ts,source,dest,pop=None):
    # if just NEA into EUR, then source is EUR, dest is NEA
    # if preOOA into AFR, but we want to see ghost ancestry in EUR
    # then source is AFR, dest is Ghost, pop is EUR
    tag=str(tag)
    def combine_segs(segs, get_segs = False):
        """
        Taken from:
        https://github.com/LauritsSkov/Introgression-detection/tree/master/Simulating%20data
        """
        merged = np.empty([0, 2])
        if len(segs) == 0:
            if get_segs:
                return([])
            else:
                return(0)
        sorted_segs = segs[np.argsort(segs[:, 0]), :]
        for higher in sorted_segs:
            if len(merged) == 0:
                merged = np.vstack([merged, higher])            
            else:
                lower = merged[-1, :]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1, :] = (lower[0], upper_bound) 
                else:
                    merged = np.vstack([merged, higher])
        if get_segs:
            return(merged)
        else:
            return(np.sum(merged[:, 1] - merged[:, 0])/seq_len)
    # Keep track of which segments are actually introgressed (in this case from pop 3 into pop 1)
    sourcep=findpop(source,ts)
    destp=findpop(dest,ts)
    if pop==None:
        pop=source
        popp=sourcep
    else:
        popp=findpop(pop,ts)
    Testpopulation = ts.get_samples(popp)    
    de_seg = {i: [] for i in Testpopulation}
    for mr in ts.migrations():
        if mr.source == sourcep and mr.dest == destp:
            for tree in ts.trees(leaf_lists=True):
                if mr.left > tree.get_interval()[0]:
                    continue
                if mr.right <= tree.get_interval()[0]:
                    break
                for l in tree.leaves(mr.node):
                    if l in Testpopulation:
                        #print l, mr
                        de_seg[l].append(tree.get_interval())
    true_de_segs = [combine_segs(np.array(de_seg[i]), True) for i in sorted(de_seg.keys())]
    return true_de_segs

def write_true_introgressed_segs(tag,ts, source,dest,true_de_segs):
    tag=str(tag)
    with open(tag+source+dest+'.anc', 'w') as fo:
        for variant in ts.variants():
            row = ''
            pos = int(variant.site.position)
            for haplotype, archaic_segments in enumerate(true_de_segs):
                keep = '0'
                for archaic_segment in archaic_segments:
                    start = int(archaic_segment[0])
                    end = int(archaic_segment[1])
                    if (start <= pos) and (pos <= end):
                        keep = '1'
                row = row + keep
            fo.write(row + '\n')

def write_eigenstrat(name,tag,ts,nsites):
    chrm="1"
    name=str(name)
    # write snp
    
    Population=ts.get_samples(findpop(tag,ts))
    with open(name+tag+'.snp', 'w') as fo:
        for variant in ts.variants():
            pos = int(variant.site.position)
            dist = str(pos / nsites)
            pos = str(pos)
            row = '\t'.join([chrm + ':' + pos, chrm, dist, pos, 'A', 'G'])
            fo.write(row + '\n')
    # write geno
    with open(name+tag+'.geno', 'w') as fo:
        for variant in ts.variants():
            row = ''.join([str(x) for x in variant.genotypes[Population]]) + '\n'
            fo.write(row)
    # write ind
    with open(name+tag+'.ind', 'w') as fo:
        for i, x in enumerate(Population):
            row = "0:" + str(i) + "\tU\t"+tag+"\n"
            fo.write(row)
            
            
def main(args):
    TS = simulate_under_prufer(args.nsites,args.demes,args.samples,args.seed,args.rep,rec_rate=args.rec_rate,mu=args.mut_rate)
    ts=next(TS)
    # ts = run_simulation(args.nsites, 1)
    # i=0
    # for ts in TS:
      # if i % 100 == 0:
        # print("Loop number:", i)
    write_eigenstrat(args.tag,args.population,ts, args.nsites)
    write_eigenstrat(args.tag,args.reference,ts, args.nsites)
    write_eigenstrat(args.tag,args.introgression,ts, args.nsites)
    if args.target==None:
        true_de_segs = get_true_introgressed_segments(args.tag,ts,args.population,args.introgression)
        write_true_introgressed_segs(args.tag,ts,args.population,args.introgression,true_de_segs)
    else:
        write_eigenstrat(args.tag,args.target,ts, args.nsites)
        true_de_segs = get_true_introgressed_segments(args.tag,ts,args.target,args.introgression,args.population)
        write_true_introgressed_segs(args.tag,ts,args.population,args.introgression,true_de_segs)
        # true_de_segs2 = get_true_introgressed_segments(args.tag,ts,args.target,args.introgression,args.target)
        # write_true_introgressed_segs(args.tag,ts,args.target,args.introgression,true_de_segs2)
      # i+=1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    parser.add_argument("-t", "--tag",
                        type=str, default="out",
                        help="tag")
    parser.add_argument("-d", "--demes",
                        type=str,
                        help="demes file for demography")
    parser.add_argument("-s", "--samples",
                        type=str, 
                        help="sample json file")
    parser.add_argument("-p", "--population",
                        type=str, 
                        help="Population to investigate introgression ind")
    parser.add_argument("-f", "--reference",
                    type=str, 
                    help="Population used as reference")
    parser.add_argument("-ns", "--nsites",
                        type=int, default=1000000,
                        help="No. sites per locus (window size in ArchIE)")
    parser.add_argument("--rec_rate",
                        type=float, default=1e-8,
                        help="No. sites per locus (window size in ArchIE)")
    parser.add_argument("--mut_rate",
                        type=float, default=1.25e-8,
                        help="No. sites per locus (window size in ArchIE)")
    parser.add_argument("-i", "--introgression",
                        type=str,
                        help="introgresing population")
    parser.add_argument("-r", "--rep",
                        type=int, default=1,
                        help="repetitions")
    parser.add_argument("-x", "--seed",
                        type=int,
                        help="Random seed")
    parser.add_argument( "--target",
                        type=str,
                        help="introgression target")
    main(parser.parse_args())