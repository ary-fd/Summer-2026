import numpy as np
import glob
import pandas as pd

# Load all simulation results
param_files = sorted(glob.glob("abc_results_10M/params_task*.npy"))
csfs_files  = sorted(glob.glob("abc_results_10M/csfs_task*.npy"))

all_params = np.vstack([np.load(f) for f in param_files])
all_csfs   = np.vstack([np.load(f) for f in csfs_files])

# Normalize CSFS rows to sum to 1
totals = all_csfs.sum(axis=1, keepdims=True)
all_csfs_norm = all_csfs / totals

# Save as CSV
pd.DataFrame(all_params, columns=["split_time","admix_time","admix_frac","ghost_ne"])\
  .to_csv("abc_results_10M/params.csv", index=False)

pd.DataFrame(all_csfs_norm)\
  .to_csv("abc_results_10M/csfs_sim.csv", index=False)

# Also save observed CSFS
obs = np.load("csfs/averaged_csfs.npy")
obs_norm = obs / obs.sum()
pd.DataFrame([obs_norm]).to_csv("abc_results_10M/csfs_observed.csv", index=False)

print(f"Saved {all_params.shape[0]} simulations")