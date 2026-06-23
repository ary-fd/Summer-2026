"""
Filter EIGENSTRAT SNPs for derived alleles shared between Africans and Neanderthals.

Takes SEPARATE EIGENSTRAT files per population (as produced by write_eigenstrat()
in the simulation script) plus the .anc introgression indicator file produced by
write_true_introgressed_segs().

File formats
------------
{tag}{pop}.geno  : one row per SNP, one character per individual (0/1/2/9)
{tag}{pop}.snp   : snp_id, chrom, gen_pos, phys_pos, ref, alt   (tab-separated)
{tag}{pop}.ind   : sample_id, sex, population                    (tab-separated)
{tag}{pop1}{pop2}.anc : one row per SNP; each character = 0 or 1 per African
                        HAPLOTYPE indicating true introgression at that site.

Ancestral allele
----------------
In the simulation all SNPs have ref='A' (ancestral) and alt='G' (derived),
so the ancestral allele is taken as ref from the .snp file.
The derived allele is alt.

Returns
-------
geno_array : np.ndarray, shape (n_snps_passing, n_afr + n_nea), dtype int8
    Alt-allele dosage (0/1/2); -1 = missing.
    Columns: African individuals first, then Neanderthal individuals.
snp_info   : list of dicts
    id, chrom, gen_pos, phys_pos, ref, alt, ancestral
kept_individuals : list of (sample_id, sex, population)
    Column metadata for geno_array.
anc_matrix : np.ndarray or None, shape (n_snps_passing, n_afr_haplotypes), dtype int8
    True introgression indicator for passing SNPs (from the .anc file).
    None if no .anc file was provided.
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
                "ref":      parts[4] if len(parts) > 4 else "A",
                "alt":      parts[5] if len(parts) > 5 else "G",
            })
    return snps


def parse_geno(geno_file):
    """Return list of row strings (one per SNP)."""
    with open(geno_file) as f:
        return [line.rstrip("\n") for line in f]


def parse_anc(anc_file):
    """
    Parse the .anc introgression indicator file.

    Each row corresponds to a SNP; each character is '0' or '1' per
    African haplotype (not per individual — the simulation writes one
    character per haploid sample / per haplotype in the population).

    Returns a list of strings, one per SNP.
    """
    with open(anc_file) as f:
        return [line.rstrip("\n") for line in f]


# ---------------------------------------------------------------------------
# Per-site helpers
# ---------------------------------------------------------------------------

GENO_TO_INT = {"0": 0, "1": 1, "2": 2, "9": -1}


def _row_to_array(row, n_individuals):
    """Convert a .geno row string to an int8 array of length n_individuals."""
    return np.array(
        [GENO_TO_INT.get(row[i], -1) for i in range(n_individuals)],
        dtype=np.int8,
    )


def has_derived_allele(row, n_individuals, ancestral, ref):
    """
    Return True if any individual in this .geno row carries the derived allele.

    ancestral == ref  →  derived = alt  →  any dosage > 0 means carrying derived
    ancestral == alt  →  derived = ref  →  any dosage < 2 means carrying derived
    """
    if ancestral == ref:
        return any(row[i] in ("1", "2") for i in range(n_individuals))
    else:
        return any(row[i] in ("0", "1") for i in range(n_individuals))


# ---------------------------------------------------------------------------
# Main filter function
# ---------------------------------------------------------------------------

def filter_derived_shared(
    african_geno,
    african_snp,
    african_ind,
    neanderthal_geno,
    neanderthal_snp,
    neanderthal_ind,
    anc_file=None,
    verbose=True,
):
    """
    Filter SNPs for derived alleles shared by Africans and Neanderthals.

    Parameters
    ----------
    african_geno, african_snp, african_ind : str
        Paths to the African EIGENSTRAT files.
    neanderthal_geno, neanderthal_snp, neanderthal_ind : str
        Paths to the Neanderthal EIGENSTRAT files.
    anc_file : str or None
        Path to the .anc introgression indicator file (produced by
        write_true_introgressed_segs). Optional — if None, anc_matrix
        in the return value will also be None.
    verbose : bool

    Returns
    -------
    geno_array : np.ndarray, shape (n_passing, n_afr + n_nea), dtype int8
        Alt-allele dosage (0/1/2; -1=missing).
        Columns: African individuals first, then Neanderthal.
    snp_info : list of dict
        One entry per passing SNP: id, chrom, gen_pos, phys_pos, ref, alt, ancestral
    kept_individuals : list of (sample_id, sex, population)
        Column metadata matching geno_array.
    anc_matrix : np.ndarray or None, shape (n_passing, n_afr_haplotypes), dtype int8
        True introgression indicator rows for passing SNPs.
        None if anc_file was not provided.
    """
    if verbose:
        print("Reading African files ...")
    afr_inds  = parse_ind(african_ind)
    afr_snps  = parse_snp(african_snp)
    afr_geno  = parse_geno(african_geno)
    n_afr     = len(afr_inds)

    if verbose:
        print(f"  {n_afr} African individuals, {len(afr_snps)} SNPs")

    if verbose:
        print("Reading Neanderthal files ...")
    nea_inds  = parse_ind(neanderthal_ind)
    nea_snps  = parse_snp(neanderthal_snp)
    nea_geno  = parse_geno(neanderthal_geno)
    n_nea     = len(nea_inds)

    if verbose:
        print(f"  {n_nea} Neanderthal individuals, {len(nea_snps)} SNPs")

    # Validate SNP counts match
    if len(afr_snps) != len(nea_snps):
        raise ValueError(
            f"African .snp has {len(afr_snps)} SNPs but Neanderthal .snp has "
            f"{len(nea_snps)} SNPs — files must cover the same sites."
        )
    n_snps = len(afr_snps)

    anc_rows = None
    if anc_file is not None:
        if verbose:
            print("Reading .anc file ...")
        anc_rows = parse_anc(anc_file)
        if len(anc_rows) != n_snps:
            raise ValueError(
                f".anc file has {len(anc_rows)} rows but expected {n_snps}."
            )

    if verbose:
        print("Filtering SNPs ...")

    passing_afr   = []
    passing_nea   = []
    passing_snps  = []
    passing_anc   = []

    for i in range(n_snps):
        snp = afr_snps[i]
        ref = snp["ref"]
        alt = snp["alt"]

        if ref in ("?", "X") or alt in ("?", "X"):
            continue

        # Ancestral allele = ref (msprime convention: ref is ancestral)
        ancestral = ref

        afr_row = afr_geno[i]
        nea_row = nea_geno[i]

        if (has_derived_allele(afr_row, n_afr, ancestral, ref) and
                has_derived_allele(nea_row, n_nea, ancestral, ref)):

            passing_afr.append(_row_to_array(afr_row, n_afr))
            passing_nea.append(_row_to_array(nea_row, n_nea))
            passing_snps.append({**snp, "ancestral": ancestral})
            if anc_rows is not None:
                passing_anc.append(anc_rows[i])

    n_passing = len(passing_snps)
    if verbose:
        print(f"  SNPs passing filter: {n_passing} / {n_snps}")

    if n_passing == 0:
        geno_array = np.empty((0, n_afr + n_nea), dtype=np.int8)
        anc_matrix = np.empty((0, 0), dtype=np.int8) if anc_rows is not None else None
    else:
        afr_block = np.vstack(passing_afr)   # (n_passing, n_afr)
        nea_block = np.vstack(passing_nea)   # (n_passing, n_nea)
        geno_array = np.hstack([afr_block, nea_block])  # (n_passing, n_afr + n_nea)

        if anc_rows is not None:
            # Each row in passing_anc is a string of '0'/'1' per African haplotype
            anc_matrix = np.array(
                [[int(c) for c in row] for row in passing_anc],
                dtype=np.int8,
            )
        else:
            anc_matrix = None

    kept_individuals = afr_inds + nea_inds
    return geno_array, passing_snps, kept_individuals, anc_matrix


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter SNPs for derived alleles shared by Africans and "
                    "Neanderthals using separate EIGENSTRAT files per population."
    )
    parser.add_argument("--african-geno", required=True)
    parser.add_argument("--african-snp",  required=True)
    parser.add_argument("--african-ind",  required=True)
    parser.add_argument("--neanderthal-geno", required=True)
    parser.add_argument("--neanderthal-snp",  required=True)
    parser.add_argument("--neanderthal-ind",  required=True)
    parser.add_argument("--anc", default=None,
                        help="True introgression indicator .anc file (optional)")
    parser.add_argument("--save-npy", metavar="PREFIX",
                        help="Save geno_array to PREFIX_geno.npy "
                             "and anc_matrix to PREFIX_anc.npy")
    args = parser.parse_args()

    geno_array, snp_info, kept_inds, anc_matrix = filter_derived_shared(
        african_geno      = args.african_geno,
        african_snp       = args.african_snp,
        african_ind       = args.african_ind,
        neanderthal_geno  = args.neanderthal_geno,
        neanderthal_snp   = args.neanderthal_snp,
        neanderthal_ind   = args.neanderthal_ind,
        anc_file          = args.anc,
    )

    print(f"\ngeno_array shape : {geno_array.shape}  (SNPs × individuals)")
    print(f"dtype            : {geno_array.dtype}")
    if anc_matrix is not None:
        print(f"anc_matrix shape : {anc_matrix.shape}  (SNPs × African haplotypes)")

    if args.save_npy:
        np.save(args.save_npy + "_geno.npy", geno_array)
        print(f"Saved geno_array to {args.save_npy}_geno.npy")
        if anc_matrix is not None:
            np.save(args.save_npy + "_anc.npy", anc_matrix)
            print(f"Saved anc_matrix to {args.save_npy}_anc.npy")


if __name__ == "__main__":
    main()
