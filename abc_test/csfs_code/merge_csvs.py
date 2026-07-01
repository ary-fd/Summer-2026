"""
Merge per-task CSV shards into combined params.csv and csfs_sim.csv.

Usage:
    python csfs_code/merge_csvs.py --out-dir abc_results_10M --n-tasks 375
"""

import argparse
import os
import glob
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir",  type=str, required=True)
    parser.add_argument("--n-tasks",  type=int, required=True)
    args = parser.parse_args()

    params_shards = []
    csfs_shards   = []
    missing = []

    for task_id in range(1, args.n_tasks + 1):
        p = os.path.join(args.out_dir, f"params_task{task_id}.csv")
        c = os.path.join(args.out_dir, f"csfs_task{task_id}.csv")

        if not os.path.exists(p) or not os.path.exists(c):
            missing.append(task_id)
            continue

        params_shards.append(pd.read_csv(p))
        csfs_shards.append(pd.read_csv(c))

    if missing:
        print(f"WARNING: missing shards for tasks: {missing}")

    params_all = pd.concat(params_shards, ignore_index=True)
    csfs_all   = pd.concat(csfs_shards,   ignore_index=True)

    # Drop task_id and sim_id before saving — run_abc.R expects only
    # parameter columns in params.csv and only CSFS columns in csfs_sim.csv
    params_out = params_all.drop(columns=["task_id", "sim_id"])
    csfs_out   = csfs_all.drop(columns=["task_id", "sim_id"])

    params_out.to_csv(os.path.join(args.out_dir, "params.csv"),   index=False)
    csfs_out.to_csv(  os.path.join(args.out_dir, "csfs_sim.csv"), index=False)

    print(f"Merged {len(params_out)} simulations from "
          f"{args.n_tasks - len(missing)} tasks into "
          f"{args.out_dir}/params.csv and {args.out_dir}/csfs_sim.csv")

if __name__ == "__main__":
    main()