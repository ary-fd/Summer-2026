"""
Generate benchmark test cases for ABC pipeline evaluation.

Draws N_TEST parameter sets from the same priors used for the reference
table, simulates a single "pseudo-observed" dataset for each, and computes
its CSFS. Output is used as the target for ABC benchmarking (prediction
error + SBC), NOT for the reference table itself.

Usage:
    python generate_test_cases.py --n-test 300 --out-dir abc_results_10M --seed 999
"""

import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

# Reuse everything from the existing simulation pipeline so test cases are
# generated with IDENTICAL logic/assumptions to the reference table.
from abc_simulate import PRIORS, FIXED, simulate_and_compute_csfs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-test", type=int, default=300)
    parser.add_argument("--out-dir", type=str, default="abc_results_10M")
    parser.add_argument("--nsite", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=999)
    args = parser.parse_args()

    FIXED["nsite"] = args.nsite
    os.makedirs(args.out_dir, exist_ok=True)

    # Use a seed range disjoint from the reference table's task seeds
    # (those used task_id * 1000 for task_id 1..375, i.e. up to ~375000).
    # Seeding far outside that range avoids any accidental overlap.
    rng     = np.random.default_rng(seed=args.seed * 1000)
    sim_rng = np.random.default_rng(seed=args.seed * 1000 + 1)

    all_params = []
    all_csfs   = []
    n_success  = 0
    n_attempt  = 0

    while n_success < args.n_test:
        n_attempt += 1

        params = {
            key: float(rng.uniform(lo, hi))
            for key, (lo, hi) in PRIORS.items()
        }
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
            print(f"Test cases: {n_success}/{args.n_test} complete "
                  f"({n_attempt} attempts)")

    all_params = np.array(all_params, dtype=np.float64)
    all_csfs   = np.array(all_csfs,   dtype=np.float64)

    # Write as CSV with headers matching params.csv / csfs_sim.csv conventions
    param_header = "split_time,admix_time,admix_frac,ghost_ne"
    np.savetxt(
        os.path.join(args.out_dir, "params_test.csv"),
        all_params, delimiter=",", header=param_header, comments="", fmt="%.10g"
    )

    n_bins = all_csfs.shape[1]
    csfs_header = ",".join(str(i) for i in range(n_bins))
    np.savetxt(
        os.path.join(args.out_dir, "csfs_test.csv"),
        all_csfs, delimiter=",", header=csfs_header, comments="", fmt="%.10g"
    )

    print(f"Done: saved {n_success} test cases "
          f"({n_attempt} attempts total) to {args.out_dir}/"
          f"params_test.csv and csfs_test.csv")


if __name__ == "__main__":
    main()