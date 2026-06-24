"""
Plot the Conditional Site Frequency Spectrum (CSFS) as a bar chart

X-axis: derived allele frequency bin (k / n_chrom)
Y-axis: proportion of SNPs in each bin

Usage
-----
    python plot_csfs.py --csfs csfs_csfs.npy --out csfs_plot.pdf

Or as a library:
    from plot_csfs import plot_csfs
    plot_csfs(csfs_array, n_chrom, out="csfs_plot.pdf")
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def plot_csfs(
    csfs,
    n_chrom,
    out="csfs_plot.pdf",
    title="Conditional Site Frequency Spectrum",
    label="Observed CSFS",
    color="#2166ac",
    show_uniform=True,
    figsize=(7, 4),
):
    """
    Plot the CSFS as a bar histogram.

    Parameters
    ----------
    csfs : np.ndarray, shape (n_chrom - 1,)
        Raw counts from compute_csfs(). Converted to proportions internally.
    n_chrom : int
        Total number of African chromosomes (= 2 * n_individuals).
    out : str
        Output file path. Extension determines format (pdf, png, svg).
    title : str
    label : str
        Legend label for the observed CSFS bars.
    color : str
        Bar color.
    show_uniform : bool
        Overlay a dashed line at the uniform expectation (1 / (n_chrom - 1)).
    figsize : tuple
    """
    total = csfs.sum()
    if total == 0:
        raise ValueError("CSFS is all zeros — nothing to plot.")

    proportions = csfs / total
    # x positions: derived allele frequency for each bin k = 1..n_chrom-1
    freqs = np.arange(1, n_chrom) / n_chrom
    bin_width = freqs[1] - freqs[0] if len(freqs) > 1 else 1.0

    fig, ax = plt.subplots(figsize=figsize)

    ax.bar(
        freqs,
        proportions,
        width=bin_width * 0.85,
        color=color,
        alpha=0.8,
        label=label,
        zorder=2,
    )

    if show_uniform:
        uniform = 1.0 / (n_chrom - 1)
        ax.axhline(
            uniform,
            color="black",
            linestyle="--",
            linewidth=1.2,
            label="Uniform expectation",
            zorder=3,
        )

    ax.set_xlabel("Derived allele frequency in Africans", fontsize=12)
    ax.set_ylabel("Proportion of SNPs", fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.4f"))
    ax.legend(fontsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.5, zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print(f"Saved plot to {out}")
    plt.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot CSFS bar chart from a saved .npy array."
    )
    parser.add_argument("--csfs", required=True,
                        help="Path to csfs_csfs.npy produced by compute_csfs.py")
    parser.add_argument("--n-chrom", type=int, default=None,
                        help="Total African chromosome count (2 × n_individuals). "
                             "Inferred from array length if not provided.")
    parser.add_argument("--out", default="csfs_plot.pdf",
                        help="Output file (pdf/png/svg). Default: csfs_plot.pdf")
    parser.add_argument("--title", default="Conditional Site Frequency Spectrum")
    parser.add_argument("--label", default="Observed CSFS")
    parser.add_argument("--color", default="#2166ac")
    parser.add_argument("--no-uniform", action="store_true",
                        help="Omit the uniform expectation line")
    args = parser.parse_args()

    csfs = np.load(args.csfs)

    # n_chrom = len(csfs) + 1  (csfs has n_chrom-1 bins)
    n_chrom = args.n_chrom if args.n_chrom is not None else len(csfs) + 1

    print(f"Loaded CSFS: {len(csfs)} bins, {csfs.sum()} total SNPs, n_chrom={n_chrom}")

    plot_csfs(
        csfs,
        n_chrom,
        out=args.out,
        title=args.title,
        label=args.label,
        color=args.color,
        show_uniform=not args.no_uniform,
    )
