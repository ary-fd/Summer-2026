"""
Filter EIGENSTRAT SNPs for derived alleles shared between Africans and Neanderthals.

Logic:
- Chimp allele = ancestral; any differing allele = derived
- Keep SNPs where the derived allele is present in >= 1 African AND >= 1 Neanderthal
- Chimp must be non-missing at the site (required to call ancestral state)

Returns
-------
geno_array : np.ndarray, shape (n_snps, n_individuals), dtype int8
    Alt allele dosage (0/1/2); -1 = missing. Columns = African + Neanderthal individuals.
snp_info   : list of dicts with keys id, chrom, gen_pos, phys_pos, ref, alt, ancestral
individuals : list of (sample_id, sex, population) for the kept columns
"""

import argparse
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_ind(ind_file):
    individuals = []
    with open(ind_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                individuals.append((parts[0], parts[1], parts[2]))
    return individuals


def parse_snp(snp_file):
    snps = []
    with open(snp_file) as f:
        for line in f:
            parts = line.split()
            snps.append({
                "id":       parts[0],
                "chrom":    parts[1],
                "gen_pos":  parts[2],
                "phys_pos": parts[3],
                "ref":      parts[4] if len(parts) > 4 else "?",
                "alt":      parts[5] if len(parts) > 5 else "?",
            })
    return snps


def group_indices(individuals, african_pops, neanderthal_pops, chimp_pops):
    african_idx, neanderthal_idx, chimp_idx = [], [], []
    for i, (_, _, pop) in enumerate(individuals):
        if pop in african_pops:
            african_idx.append(i)
        elif pop in neanderthal_pops:
            neanderthal_idx.append(i)
        elif pop in chimp_pops:
            chimp_idx.append(i)
    return african_idx, neanderthal_idx, chimp_idx


# ---------------------------------------------------------------------------
# Per-site helpers
# ---------------------------------------------------------------------------

GENO_TO_INT = {"0": 0, "1": 1, "2": 2, "9": -1}


def _geno_ints(row, indices):
    """Return int8 array of alt-dosage values for the given column indices."""
    return np.array([GENO_TO_INT.get(row[i], -1) for i in indices], dtype=np.int8)


def ancestral_allele(row, chimp_idx, ref, alt):
    """
    Infer ancestral allele from chimp genotypes.

    All non-missing chimp calls must agree on a single allele.
    Returns the ancestral allele string, or None if undetermined.
    """
    alleles = set()
    for idx in chimp_idx:
        code = row[idx]
        if code == "9":
            continue
        if code == "0":
            alleles.add(ref)
        elif code == "2":
            alleles.add(alt)
        elif code == "1":
            alleles.update([ref, alt])
    return alleles.pop() if len(alleles) == 1 else None


def has_derived_allele(row, indices, ref, alt, ancestral):
    """Return True if any individual in indices carries the derived allele."""
    derived_code = "2" if ancestral == ref else "0"
    het_code = "1"
    for idx in indices:
        c = row[idx]
        if c == derived_code or c == het_code:
            return True
    return False


# ---------------------------------------------------------------------------
# Main filter function
# ---------------------------------------------------------------------------

def filter_derived_shared(
    geno_file,
    snp_file,
    ind_file,
    african_pops=("Yoruba", "Mbuti", "Dinka", "Mandenka"),
    neanderthal_pops=("Neanderthal", "Vindija", "Altai"),
    chimp_pops=("Chimp", "PanTro"),
    verbose=True,
):
    """
    Filter EIGENSTRAT SNPs for derived alleles shared by Africans and Neanderthals.

    Parameters
    ----------
    geno_file, snp_file, ind_file : str
        Paths to EIGENSTRAT input files.
    african_pops, neanderthal_pops, chimp_pops : sequence of str
        Population labels in the .ind file for each group.
    verbose : bool
        Print progress to stdout.

    Returns
    -------
    geno_array : np.ndarray, shape (n_snps_passing, n_kept_individuals), dtype int8
        Alt allele dosage (0/1/2); -1 = missing.
        Column order: all African individuals, then all Neanderthal individuals.
    snp_info : list of dict
        One entry per passing SNP with keys:
        id, chrom, gen_pos, phys_pos, ref, alt, ancestral
    kept_individuals : list of (sample_id, sex, population)
        Metadata for each column in geno_array.
    """
    african_pops    = set(african_pops)
    neanderthal_pops = set(neanderthal_pops)
    chimp_pops      = set(chimp_pops)

    if verbose:
        print("Reading .ind ...")
    individuals = parse_ind(ind_file)
    african_idx, neanderthal_idx, chimp_idx = group_indices(
        individuals, african_pops, neanderthal_pops, chimp_pops
    )

    if verbose:
        print(f"  Africans:     {len(african_idx)}")
        print(f"  Neanderthals: {len(neanderthal_idx)}")
        print(f"  Chimps:       {len(chimp_idx)}")

    if not african_idx:
        raise ValueError("No African individuals found — check african_pops")
    if not neanderthal_idx:
        raise ValueError("No Neanderthal individuals found — check neanderthal_pops")
    if not chimp_idx:
        raise ValueError("No Chimp individuals found — check chimp_pops")

    if verbose:
        print("Reading .snp ...")
    snps = parse_snp(snp_file)
    if verbose:
        print(f"  Total SNPs: {len(snps)}")

    # Column order in output: African first, then Neanderthal
    keep_idx = african_idx + neanderthal_idx

    if verbose:
        print("Filtering .geno ...")

    passing_rows = []
    passing_snps = []

    with open(geno_file) as f:
        for snp, line in zip(snps, f):
            row = line.rstrip("\n")
            ref, alt = snp["ref"], snp["alt"]

            if ref in ("?", "X") or alt in ("?", "X"):
                continue

            anc = ancestral_allele(row, chimp_idx, ref, alt)
            if anc is None:
                continue

            if (has_derived_allele(row, african_idx, ref, alt, anc) and
                    has_derived_allele(row, neanderthal_idx, ref, alt, anc)):
                passing_rows.append(_geno_ints(row, keep_idx))
                passing_snps.append({**snp, "ancestral": anc})

    if verbose:
        print(f"  SNPs passing filter: {len(passing_snps)}")

    if not passing_rows:
        geno_array = np.empty((0, len(keep_idx)), dtype=np.int8)
    else:
        geno_array = np.vstack(passing_rows)  # shape: (n_snps, n_individuals)

    kept_individuals = [individuals[i] for i in keep_idx]
    return geno_array, passing_snps, kept_individuals


# ---------------------------------------------------------------------------
# CLI wrapper
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter EIGENSTRAT SNPs for derived alleles shared by "
                    "Africans and Neanderthals. Prints array shape and SNP count."
    )
    parser.add_argument("--geno", required=True)
    parser.add_argument("--snp",  required=True)
    parser.add_argument("--ind",  required=True)
    parser.add_argument("--african-pops",     nargs="+",
                        default=["Yoruba", "Mbuti", "Dinka", "Mandenka"])
    parser.add_argument("--neanderthal-pops", nargs="+",
                        default=["Neanderthal", "Vindija", "Altai"])
    parser.add_argument("--chimp-pops",       nargs="+",
                        default=["Chimp", "PanTro"])
    parser.add_argument("--save-npy", metavar="PATH",
                        help="Optional: save geno_array to a .npy file")
    args = parser.parse_args()

    geno_array, snp_info, kept_inds = filter_derived_shared(
        geno_file        = args.geno,
        snp_file         = args.snp,
        ind_file         = args.ind,
        african_pops     = args.african_pops,
        neanderthal_pops = args.neanderthal_pops,
        chimp_pops       = args.chimp_pops,
    )

    print(f"\ngeno_array shape : {geno_array.shape}  (SNPs × individuals)")
    print(f"dtype            : {geno_array.dtype}")
    print(f"Kept individuals : {len(kept_inds)}")
    print(f"  {[ind[0] for ind in kept_inds]}")

    if args.save_npy:
        np.save(args.save_npy, geno_array)
        print(f"Saved array to {args.save_npy}")


if __name__ == "__main__":
    main()
