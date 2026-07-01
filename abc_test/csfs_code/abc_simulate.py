"""
ABC simulation script for ghost archaic introgression.
Draws parameters from prior, runs msprime simulation, computes CSFS, saves results.

Usage (called by array job):
    python abc_simulate.py --task-id 1 --n-sims 200 --out-dir abc_results
"""

import argparse
import csv
import sys
import os
import numpy as np
import msprime
import demes

sys.path.insert(0, os.path.dirname(__file__))
from compute_csfs import compute_csfs
from filter_derived_alleles import filter_derived_shared

# ---------------------------------------------------------------------------
# Prior ranges (uniform)
# ---------------------------------------------------------------------------
# All times in generations (gen_time = 20 years)
PRIORS = {
    "split_time":       (6_201,   100_000),   # in generations
    "admix_time":       (100,       6_200),    # in generations
    "admix_frac":       (0.02,      0.50),
    "ghost_ne":         (10_000,   50_000),
}

# Fixed parameters (from Prufer et al. 2015 / your YAML)
FIXED = {
    "chimp_ne":         10_000,
    "ancestral_ne":     10_000,
    "africa_ne":        10_000,
    "neanderthal_ne":   10_000,
    "chimp_split":      300_000,   # generations
    "human_nea_split":   20_000,   # generations
    "nea_extinction":    2_000,    # generations
    "nsite":           1_000_000,
    "mu":               1.25e-8,
    "rec_rate":         1e-8,
    "n_africa_individuals":     50,  # diploid Africans sampled
    "n_neanderthal_haplotypes": 1,    # single haploid Neanderthal genome
}
FIXED["n_africa_haplotypes"] = 2 * FIXED["n_africa_individuals"]  # = 100


# ---------------------------------------------------------------------------
# Build demes graph from parameters
# ---------------------------------------------------------------------------
def build_graph(split_time, admix_time, admix_frac, ghost_ne):
    """Build demes graph with given ghost population parameters."""
    b = demes.Builder(time_units="generations")

    b.add_deme("Chimp", epochs=[
        dict(end_time=0, start_size=FIXED["chimp_ne"])
    ])
    b.add_deme("Ancestral", ancestors=["Chimp"],
               start_time=FIXED["chimp_split"], epochs=[
        dict(end_time=FIXED["human_nea_split"], start_size=FIXED["ancestral_ne"])
    ])
    b.add_deme("Ghost", ancestors=["Ancestral"],
               start_time=split_time, epochs=[
        dict(end_time=0, start_size=ghost_ne)
    ])
    b.add_deme("Neanderthal", ancestors=["Ancestral"],
               start_time=FIXED["human_nea_split"], epochs=[
        dict(end_time=FIXED["nea_extinction"], start_size=FIXED["neanderthal_ne"])
    ])
    b.add_deme("Africa", ancestors=["Ancestral"],
               start_time=FIXED["human_nea_split"], epochs=[
        dict(end_time=0, start_size=FIXED["africa_ne"])
    ])

    b.add_pulse(sources=["Ghost"], dest="Africa",
                time=admix_time, proportions=[admix_frac])

    return b.resolve()


# ---------------------------------------------------------------------------
# Run one simulation and compute CSFS
# ---------------------------------------------------------------------------

def simulate_and_compute_csfs(params, rng, sim_rng):
    """Draw params, simulate, compute CSFS. Returns (csfs, params) or None on failure."""
    split_time = params["split_time"]
    admix_time = params["admix_time"]
    admix_frac = params["admix_frac"]
    ghost_ne   = params["ghost_ne"]

    if split_time <= FIXED["human_nea_split"]:
        return None

    try:
        graph = build_graph(split_time, admix_time, admix_frac, ghost_ne)
        demography = msprime.Demography.from_demes(graph)

        # Africa: n_africa_individuals diploid people -> n_africa_haplotypes
        # chromosomes. Neanderthal: a single haploid genome (ploidy=1 override),
        # not a diploid individual, so we get exactly 1 sampled chromosome.
        samples = [
            msprime.SampleSet(FIXED["n_africa_individuals"],
                               population="Africa", ploidy=2),
            msprime.SampleSet(FIXED["n_neanderthal_haplotypes"],
                               population="Neanderthal", ploidy=1),
        ]

        ts = msprime.sim_ancestry(
            samples=samples,
            demography=demography,
            sequence_length=FIXED["nsite"],
            recombination_rate=FIXED["rec_rate"],
            record_migrations=True,
            random_seed=int(sim_rng.integers(1, 2**31)),
        )
        ts = msprime.sim_mutations(
            ts,
            rate=FIXED["mu"],
            random_seed=int(sim_rng.integers(1, 2**31)),
        )

        # One dosage column per sampled individual (Africans: sum of their
        # 2 haplotypes, 0/1/2; Neanderthal: its single haplotype, 0/1).
        afr_geno, afr_snps = _ts_to_dosage_arrays(ts, _popid(ts, "Africa"))
        nea_geno, nea_snps = _ts_to_dosage_arrays(ts, _popid(ts, "Neanderthal"))

        if afr_geno.shape[0] == 0:
            return None

        # Filter and compute CSFS directly in memory
        geno_array, snp_info, kept_inds, _ = _filter_in_memory(
            afr_geno, afr_snps, nea_geno, nea_snps
        )

        if geno_array.shape[0] == 0:
            return None

        csfs, n_chrom = compute_csfs(
            geno_array, snp_info, kept_inds,
            n_african=FIXED["n_africa_individuals"],
        )

        return csfs

    except Exception as e:
        print(f"Simulation failed: {e}")
        return None


def _popid(ts, name):
    for pop in ts.populations():
        if pop.metadata["name"] == name:
            return pop.id
    raise ValueError(f"Population {name} not found")


def _ts_to_dosage_arrays(ts, pop_id):
    """
    Convert tree sequence to (dosage_matrix, snp_list) for a population.

    One column per *individual* (not per haplotype): dosage = number of
    alt-allele copies summed across that individual's nodes (0-2 for a
    diploid, 0-1 for a haploid, per SampleSet(ploidy=...) at sampling time).
    This matches the per-individual dosage format compute_csfs expects.
    """
    node_groups = [ind.nodes for ind in ts.individuals() if ind.population == pop_id]

    rows = []
    snps = []
    for variant in ts.variants():
        geno = variant.genotypes
        dosage = np.array([geno[nodes].sum() for nodes in node_groups], dtype=np.int8)
        rows.append(dosage)
        pos = int(variant.site.position)
        snps.append({
            "id": f"1:{pos}", "chrom": "1",
            "gen_pos": str(pos / FIXED["nsite"]),
            "phys_pos": str(pos),
            "ref": "A", "alt": "G",
        })
    if rows:
        return np.vstack(rows), snps
    return np.empty((0, len(node_groups)), dtype=np.int8), []


def _filter_in_memory(afr_geno, afr_snps, nea_geno, nea_snps):
    """Filter SNPs where derived allele present in both Africa and Neanderthal."""
    n_afr = afr_geno.shape[1]
    n_nea = nea_geno.shape[1]
    passing_afr, passing_nea, passing_snps = [], [], []

    for i, snp in enumerate(afr_snps):
        ref, alt = snp["ref"], snp["alt"]
        ancestral = ref  # ref is always ancestral in simulation

        afr_row = afr_geno[i]
        nea_row = nea_geno[i]

        # derived = alt (G); present if any value > 0
        if np.any(afr_row > 0) and np.any(nea_row > 0):
            passing_afr.append(afr_row)
            passing_nea.append(nea_row)
            passing_snps.append({**snp, "ancestral": ancestral})

    if not passing_snps:
        empty = np.empty((0, n_afr + n_nea), dtype=np.int8)
        return empty, [], [], None

    afr_block = np.vstack(passing_afr)
    nea_block = np.vstack(passing_nea)
    geno_array = np.hstack([afr_block, nea_block])

    # kept_individuals metadata
    kept_inds = (
        [(f"0:{i}", "U", "Africa") for i in range(n_afr)] +
        [(f"0:{i}", "U", "Neanderthal") for i in range(n_nea)]
    )

    return geno_array, passing_snps, kept_inds, None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id",  type=int, required=True)
    parser.add_argument("--n-sims",   type=int, default=200)
    parser.add_argument("--out-dir", type=str, default="abc_results")
    parser.add_argument("--nsite", type=int, default=1000000)
    args = parser.parse_args()

    FIXED["nsite"] = args.nsite

    os.makedirs(args.out_dir, exist_ok=True)

    # Per-task files (one shard per array-job task, matching the old .npy
    # layout) so concurrent tasks never write to the same file. Concatenate
    # the shards afterward (e.g. `csvstack` or pandas.concat over the glob)
    # to get the two combined CSVs.
    params_path = os.path.join(args.out_dir, f"params_task{args.task_id}.csv")
    csfs_path   = os.path.join(args.out_dir, f"csfs_task{args.task_id}.csv")

    param_names = list(PRIORS.keys())

    # Separate RNGs for parameter draws and simulations
    rng     = np.random.default_rng(seed=args.task_id * 1000)
    sim_rng = np.random.default_rng(seed=args.task_id * 1000 + 1)

    n_success = 0
    n_attempt = 0
    csfs_header_written = False

    with open(params_path, "w", newline="") as params_f, \
         open(csfs_path,   "w", newline="") as csfs_f:

        params_writer = csv.writer(params_f)
        csfs_writer   = csv.writer(csfs_f)

        params_writer.writerow(["task_id", "sim_id"] + param_names)

        while n_success < args.n_sims:
            n_attempt += 1

            # Draw parameters from uniform priors
            params = {
                key: float(rng.uniform(lo, hi))
                for key, (lo, hi) in PRIORS.items()
            }
            # Split time and ghost Ne should be integers
            params["split_time"] = int(params["split_time"])
            params["ghost_ne"]   = int(params["ghost_ne"])

            csfs = simulate_and_compute_csfs(params, rng, sim_rng)
            if csfs is None:
                continue

            if not csfs_header_written:
                csfs_writer.writerow(
                    ["task_id", "sim_id"] + [f"csfs_k{k+1}" for k in range(len(csfs))]
                )
                csfs_header_written = True

            params_writer.writerow(
                [args.task_id, n_success] + [params[key] for key in param_names]
            )
            csfs_writer.writerow([args.task_id, n_success] + list(csfs))

            n_success += 1

            if n_success % 20 == 0:
                params_f.flush()
                csfs_f.flush()
                print(f"Task {args.task_id}: {n_success}/{args.n_sims} complete "
                      f"({n_attempt} attempts)")

    print(f"Task {args.task_id}: saved {n_success} simulations "
          f"({n_attempt} attempts total) to {params_path} and {csfs_path}")


if __name__ == "__main__":
    main()
