import numpy as np
import glob

csfs_files = sorted(glob.glob("csfs/rep*_csfs.npy"))
csfs_list = [np.load(f) for f in csfs_files]

# average then normalize
csfs_mean = np.mean(csfs_list, axis=0)
csfs_norm = csfs_mean / csfs_mean.sum()

np.save("csfs/averaged_csfs.npy", csfs_norm)