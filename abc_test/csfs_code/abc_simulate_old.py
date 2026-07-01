"""
ABC simulation script for ghost archaic introgression.
Draws parameters from prior, runs msprime simulation, computes CSFS, saves results.

Usage (called by array job):
    python abc_simulate.py --task-id 1 --n-sims 200 --out-dir abc_results
"""

import argparse
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
    "split_time":       (5_000,   100_000),   # 100ka to 2Ma in generations
    "admix_time":       (100,       6_200),    # ~2ka to 124ka in generations
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
    "n_africa":         100,       # haplotypes
    "n_neanderthal":    1,         # haplotype
}


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

    # Validate: split_time must be > human_nea_split, admix_time < nea_extinction
    if split_time <= FIXED["human_nea_split"]:
        return None
    if admix_time >= FIXED["nea_extinction"]:
        return None
    if admix_time >= split_time:
        return None

    try:
        graph = build_graph(split_time, admix_time, admix_frac, ghost_ne)
        demography = msprime.Demography.from_demes(graph)

        samples = {
            "Africa":      FIXED["n_africa"],
            "Neanderthal": FIXED["n_neanderthal"],
        }

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

        # Write temporary eigenstrat files
        afr_pop  = ts.samples(population=_popid(ts, "Africa"))
        nea_pop  = ts.samples(population=_popid(ts, "Neanderthal"))

        afr_geno, afr_snps = _ts_to_arrays(ts, afr_pop)
        nea_geno, nea_snps = _ts_to_arrays(ts, nea_pop)

        if afr_geno.shape[0] == 0:
            return None

        # Filter and compute CSFS directly in memory
        geno_array, snp_info, kept_inds, _ = _filter_in_memory(
            afr_geno, afr_snps, nea_geno, nea_snps,
            FIXED["n_africa"], FIXED["n_neanderthal"]
        )

        if geno_array.shape[0] == 0:
            return None

        csfs, n_chrom = compute_csfs(
            geno_array, snp_info, kept_inds,
            n_african=FIXED["n_africa"],
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


def _ts_to_arrays(ts, pop_samples):
    """Convert tree sequence to (geno_matrix, snp_list) for a population."""
    rows = []
    snps = []
    for variant in ts.variants():
        geno = variant.genotypes[pop_samples]
        rows.append(geno.astype(np.int8))
        pos = int(variant.site.position)
        snps.append({
            "id": f"1:{pos}", "chrom": "1",
            "gen_pos": str(pos / FIXED["nsite"]),
            "phys_pos": str(pos),
            "ref": "A", "alt": "G",
        })
    if rows:
        return np.vstack(rows), snps
    return np.empty((0, len(pop_samples)), dtype=np.int8), []


def _filter_in_memory(afr_geno, afr_snps, nea_geno, nea_snps, n_afr, n_nea):
    """Filter SNPs where derived allele present in both Africa and Neanderthal."""
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

    # Separate RNGs for parameter draws and simulations
    rng     = np.random.default_rng(seed=args.task_id * 1000)
    sim_rng = np.random.default_rng(seed=args.task_id * 1000 + 1)

    all_params = []
    all_csfs   = []
    n_success  = 0
    n_attempt  = 0

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

        all_params.append([
            params["split_time"],
            params["admix_time"],
            params["admix_frac"],
            params["ghost_ne"],
        ])
        all_csfs.append(csfs)
        n_success += 1

        if n_success % 20 == 0:
            print(f"Task {args.task_id}: {n_success}/{args.n_sims} complete "
                  f"({n_attempt} attempts)")

    all_params = np.array(all_params, dtype=np.float32)
    all_csfs   = np.array(all_csfs,   dtype=np.float32)

    np.save(f"{args.out_dir}/params_task{args.task_id}.npy", all_params)
    np.save(f"{args.out_dir}/csfs_task{args.task_id}.npy",   all_csfs)
    print(f"Task {args.task_id}: saved {n_success} simulations "
          f"({n_attempt} attempts total)")


if __name__ == "__main__":
    main()