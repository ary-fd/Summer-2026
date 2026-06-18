"""
Compute the Conditional Site Frequency Spectrum (CSFS) as defined in
Durvasula & Sankararaman 2020 (Sci. Adv. 6, eaax5097).

Definition (from paper, Materials & Methods):
    CSFS_{pop1, pop2}[k] = number of SNPs at which the derived allele
    is present on k chromosomes in a sample of n total chromosomes in
    pop1 (Africans), while a single randomly sampled chromosome from
    the archaic outgroup pop2 (Neanderthal) carries the derived allele.

    k ranges over {1, 2, ..., n-1}  (excludes fixed sites).

Under a neutral model with no admixture, the CSFS is expected to be
uniform. Deviation (U-shape) is evidence of ghost archaic introgression.

Input
-----
The numpy array, snp_info list, and kept_individuals list returned by
filter_derived_shared() in filter_derived_alleles.py.

Usage
-----
    from filter_derived_alleles import filter_derived_shared
    from compute_csfs import compute_csfs, normalize_csfs

    geno, snps, inds = filter_derived_shared(...)
    csfs, n_chrom = compute_csfs(geno, snps, inds)
    csfs_norm = normalize_csfs(csfs)
"""

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _infer_n_african(kept_individuals):
    """
    Infer how many of the kept individuals are African vs Neanderthal.

    filter_derived_shared() guarantees the order is [afr_inds..., nea_inds...].
    We detect the boundary by finding where the population label changes to a
    known archaic label. Falls back to raising a clear error with instructions.
    """
    archaic_pops = {"Neanderthal", "Vindija", "Altai", "Denisova"}
    for i, (_, _, pop) in enumerate(kept_individuals):
        if pop in archaic_pops:
            return i  # first archaic = boundary; everything before is African
    # All individuals are non-archaic
    return len(kept_individuals)


def _derived_dosage(dosage_row, ancestral, ref, alt):
    """
    Convert haploid alt-dosage row to derived-allele dosage row.
    
    Haploid encoding: 0 = ancestral allele, 1 = alt allele, -1 = missing.
    If alt is derived: derived_dosage = dosage as-is.
    If ref is derived: derived_dosage = 1 - dosage (flip 0 and 1).
    Missing (-1) stays -1.
    """
    if ancestral == ref:
        # alt is derived — dosage already counts derived copies
        return dosage_row.copy()
    else:
        # ref is derived — flip: 0→1, 1→0, missing stays -1
        out = np.where(dosage_row == -1, -1, 1 - dosage_row)
        return out.astype(np.int8)


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def compute_csfs(
    geno_array,
    snp_info,
    kept_individuals,
    n_african=None,
    min_african_chrom=None,
    condition_on_archaic="any",
):
    """
    Compute the CSFS from the filtered genotype array.

    Parameters
    ----------
    geno_array : np.ndarray, shape (n_snps, n_individuals), dtype int8
        Output of filter_derived_shared(). Alt-allele dosage; -1 = missing.
        Columns: African individuals first, then Neanderthal.
    snp_info : list of dict
        Output of filter_derived_shared(). Must include 'ref', 'alt', 'ancestral'.
    kept_individuals : list of (sample_id, sex, population)
        Column metadata for geno_array.
    n_african : int, optional
        Number of African individuals (= first n_african columns).
        Inferred from kept_individuals if not provided.
    min_african_chrom : int, optional
        Minimum non-missing African chromosomes required to use a SNP.
        Defaults to all African chromosomes (no missing data allowed).
    condition_on_archaic : {"any", "hom", "random"}
        How to handle heterozygous archaic sites when determining whether
        the archaic carries the derived allele:
          "any"    – include SNP if any archaic chromosome is derived (dosage >= 1).
                     This is what the filter already enforces; all SNPs pass.
          "hom"    – include SNP only if archaic is homozygous derived (dosage == 2).
          "random" – for each heterozygous archaic individual, randomly sample one
                     chromosome; include SNP if that chromosome is derived.

    Returns
    -------
    csfs : np.ndarray, shape (n_chrom - 1,), dtype int64
        csfs[k-1] = count of SNPs with k derived African chromosomes,
        for k = 1, 2, ..., n_chrom - 1.
        (k=0 and k=n_chrom are excluded as they are not polymorphic.)
    n_chrom : int
        Total diploid African chromosome count used (= 2 * n_african).
        If per-SNP missing data is allowed, this is the modal or max value;
        see Returns notes below.

    Notes
    -----
    When min_african_chrom is None (default), only SNPs with zero missing
    African data are used, and n_chrom is constant across all SNPs.
    When min_african_chrom < 2*n_african, the per-SNP chromosome count
    varies; the CSFS is then indexed by derived allele *count* (not frequency),
    and n_chrom returns the maximum observed. For frequency-based analyses
    (e.g., plotting as proportion), use normalize_csfs() which divides each
    bin by total SNP count.
    """
    if n_african is None:
        n_african = _infer_n_african(kept_individuals)

    n_neanderthal = len(kept_individuals) - n_african
    n_chrom_max   = n_african

    if min_african_chrom is None:
        min_african_chrom = n_chrom_max  # require no missing data by default

    afr_cols = geno_array[:, :n_african]          # (n_snps, n_afr)
    nea_cols = geno_array[:, n_african:]          # (n_snps, n_nea)

    csfs = np.zeros(n_chrom_max - 1, dtype=np.int64)  # bins k=1..n_chrom_max-1

    rng = np.random.default_rng(seed=42)

    for i, snp in enumerate(snp_info):
        ref = snp["ref"]
        alt = snp["alt"]
        anc = snp["ancestral"]

        # ---- Archaic condition ----------------------------------------
        nea_row_raw = nea_cols[i]  # alt-dosage for Neanderthal individuals

        # Convert to derived dosage for archaic
        nea_derived = _derived_dosage(nea_row_raw, anc, ref, alt)
        nea_nonmiss = nea_derived[nea_derived != -1]

        if len(nea_nonmiss) == 0:
            continue  # no archaic data

        if condition_on_archaic == "any":
            # Already guaranteed by the filter; just a safety check
            if not np.any(nea_nonmiss > 0):
                continue
        elif condition_on_archaic == "hom":
            if not np.any(nea_nonmiss == 2):
                continue
        elif condition_on_archaic == "random":
            # For each archaic individual, randomly sample one chromosome
            passes = False
            for d in nea_nonmiss:
                if d == 2:
                    passes = True; break
                if d == 1 and rng.random() < 0.5:
                    passes = True; break
            if not passes:
                continue

        # ---- African derived allele count -----------------------------
        afr_row_raw = afr_cols[i]
        afr_derived = _derived_dosage(afr_row_raw, anc, ref, alt)

        missing_mask = afr_derived == -1
        n_chrom_this = 2 * np.sum(~missing_mask.reshape(n_african))
        # Each individual contributes 2 chromosomes; missing individuals
        # (all -1 for that individual) reduce the count.
        # Simpler: count non-missing allele slots directly.
        afr_flat = afr_derived  # shape (n_african,) — diploid dosage values
        # count of non-missing chromosomes = 2 * non-missing individuals
        # (EIGENSTRAT dosage is per-individual, not per-chromosome)
        n_nonmiss_inds = np.sum(afr_flat != -1)
        n_chrom_this = n_nonmiss_inds  # not 2 * n_nonmiss_inds

        if n_chrom_this < min_african_chrom:
            continue

        # Derived allele count k (sum of dosage values, ignoring missing)
        k = int(np.sum(afr_flat[afr_flat != -1]))

        # Exclude fixed sites (k=0 or k=n_chrom_this)
        if k == 0 or k == n_chrom_this:
            continue

        # Place in CSFS — if n_chrom varies, only bin up to n_chrom_max-1
        if 1 <= k <= n_chrom_max - 1:
            csfs[k - 1] += 1

    return csfs, n_chrom_max


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def normalize_csfs(csfs):
    """Return CSFS as proportions (sum to 1). Avoids divide-by-zero."""
    total = csfs.sum()
    if total == 0:
        return csfs.astype(float)
    return csfs / total


def csfs_summary(csfs, n_chrom):
    """
    Print a summary table of the CSFS.

    Parameters
    ----------
    csfs : np.ndarray, shape (n_chrom - 1,)
    n_chrom : int
    """
    total = csfs.sum()
    props = normalize_csfs(csfs)
    freq_bins = np.arange(1, n_chrom) / n_chrom  # derived allele frequency

    print(f"{'k':>5} {'freq':>8} {'count':>10} {'proportion':>12}")
    print("-" * 40)
    for k, (cnt, prop, freq) in enumerate(zip(csfs, props, freq_bins), start=1):
        print(f"{k:>5} {freq:>8.3f} {cnt:>10} {prop:>12.6f}")
    print("-" * 40)
    print(f"{'Total':>5} {'':>8} {total:>10} {'1.000000':>12}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import sys
    sys.path.insert(0, ".")
    from filter_derived_alleles import filter_derived_shared

    parser = argparse.ArgumentParser(
        description="Filter SNPs and compute CSFS (Durvasula & Sankararaman 2020)."
    )
    parser.add_argument("--african-geno", required=True)
    parser.add_argument("--african-snp",  required=True)
    parser.add_argument("--african-ind",  required=True)
    parser.add_argument("--neanderthal-geno", required=True)
    parser.add_argument("--neanderthal-snp",  required=True)
    parser.add_argument("--neanderthal-ind",  required=True)
    parser.add_argument("--anc", default=None,
                        help="True introgression .anc file (optional, passed through)")
    parser.add_argument("--condition-on-archaic",
                        choices=["any", "hom", "random"], default="any",
                        help="How to treat het Neanderthal sites (default: any)")
    parser.add_argument("--save-npy", metavar="PREFIX",
                        help="Save csfs array to PREFIX_csfs.npy")
    args = parser.parse_args()

    geno_array, snp_info, kept_inds, anc_matrix = filter_derived_shared(
        african_geno     = args.african_geno,
        african_snp      = args.african_snp,
        african_ind      = args.african_ind,
        neanderthal_geno = args.neanderthal_geno,
        neanderthal_snp  = args.neanderthal_snp,
        neanderthal_ind  = args.neanderthal_ind,
        anc_file         = args.anc,
    )

    # African individuals are always the first block in kept_inds;
    # n_african is inferred automatically from the population label boundary.
    n_african = None  # let compute_csfs call _infer_n_african

    print("\nComputing CSFS ...")
    csfs, n_chrom = compute_csfs(
        geno_array,
        snp_info,
        kept_inds,
        n_african=n_african,
        condition_on_archaic=args.condition_on_archaic,
    )

    print(f"\nCSFS (n_chrom = {n_chrom}, total SNPs = {csfs.sum()})\n")
    csfs_summary(csfs, n_chrom)

    if args.save_npy:
        np.save(args.save_npy + "_csfs.npy", csfs)
        print(f"\nSaved CSFS to {args.save_npy}_csfs.npy")
