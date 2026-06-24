"""
Train neural network on ABC simulation results and predict parameters
from observed CSFS.

Usage:
    python train_nn.py --results-dir abc_results --observed csfs/averaged_csfs.npy
"""

import argparse
import os
import numpy as np
import glob
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
import joblib


PARAM_NAMES = ["split_time", "admix_time", "admix_frac", "ghost_ne"]
TRUE_PARAMS  = {
    "split_time": 31250,
    "admix_time": 2150,
    "admix_frac": 0.11,
    "ghost_ne":   25000,
}


def load_results(results_dir):
    param_files = sorted(glob.glob(f"{results_dir}/params_task*.npy"))
    csfs_files  = sorted(glob.glob(f"{results_dir}/csfs_task*.npy"))

    print(f"Found {len(param_files)} task result files")

    all_params = np.vstack([np.load(f) for f in param_files])
    all_csfs   = np.vstack([np.load(f) for f in csfs_files])

    print(f"Total simulations: {all_params.shape[0]}")
    print(f"CSFS shape: {all_csfs.shape}")
    return all_params, all_csfs


def normalize_csfs(csfs_matrix):
    """Normalize each CSFS row to sum to 1."""
    totals = csfs_matrix.sum(axis=1, keepdims=True)
    totals = np.where(totals == 0, 1, totals)
    return csfs_matrix / totals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=str, default="abc_results")
    parser.add_argument("--observed",    type=str, required=True,
                        help="Path to observed averaged CSFS .npy file")
    parser.add_argument("--out-dir",     type=str, default="abc_results")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Load all simulation results
    print("Loading simulation results ...")
    all_params, all_csfs = load_results(args.results_dir)

    # Normalize CSFS
    all_csfs_norm = normalize_csfs(all_csfs)

    # Scale parameters
    param_scaler = StandardScaler()
    all_params_scaled = param_scaler.fit_transform(all_params)

    # Train neural network
    print("Training neural network ...")
    nn = MLPRegressor(
        hidden_layer_sizes=(256, 256, 128),
        activation="relu",
        max_iter=500,
        random_state=42,
        verbose=True,
        early_stopping=True,
        validation_fraction=0.1,
    )
    nn.fit(all_csfs_norm, all_params_scaled)
    print("Training complete.")

    # Save model
    joblib.dump(nn,           f"{args.out_dir}/nn_model.pkl")
    joblib.dump(param_scaler, f"{args.out_dir}/param_scaler.pkl")
    print(f"Saved model to {args.out_dir}/nn_model.pkl")

    # Predict from observed CSFS
    print("\nPredicting parameters from observed CSFS ...")
    observed_csfs = np.load(args.observed)
    observed_norm = observed_csfs / observed_csfs.sum()
    observed_norm = observed_norm.reshape(1, -1)

    predicted_scaled = nn.predict(observed_norm)
    predicted = param_scaler.inverse_transform(predicted_scaled)[0]

    print("\nResults:")
    print(f"{'Parameter':<15} {'True':>12} {'Predicted':>12}")
    print("-" * 42)
    for name, true_val, pred_val in zip(PARAM_NAMES, TRUE_PARAMS.values(), predicted):
        print(f"{name:<15} {true_val:>12.2f} {pred_val:>12.2f}")

    np.save(f"{args.out_dir}/predicted_params.npy", predicted)


if __name__ == "__main__":
    main()