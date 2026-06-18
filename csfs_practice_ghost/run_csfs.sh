#!/bin/bash
#$ -cwd
#$ -j y
#$ -N run_csfs
#$ -l h_data=8G,h_rt=01:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/csfs_practice_ghost/logs/$JOB_NAME.o$JOB_ID
set -x

cd /u/scratch/a/aryadini/Summer-2026/csfs_practice_ghost
source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

# --- paths -----------------------------------------------------------------
EIGENSTRAT=/u/scratch/a/aryadini/Summer-2026/csfs_practice_ghost/eigenstrat
SCRIPTS=/u/scratch/a/aryadini/Summer-2026/csfs_practice_ghost
OUT=/u/scratch/a/aryadini/Summer-2026/csfs_practice_ghost

# simulation parameters (must match what was used to generate the files)
TAG=ghost_sim
POP=Africa
INTROGRESS=Ghost
NEANDERTHAL=Neanderthal

# derived filenames from write_eigenstrat() and write_true_introgressed_segs()
AFR_GENO=${EIGENSTRAT}/${TAG}${POP}.geno
AFR_SNP=${EIGENSTRAT}/${TAG}${POP}.snp
AFR_IND=${EIGENSTRAT}/${TAG}${POP}.ind

NEA_GENO=${EIGENSTRAT}/${TAG}${NEANDERTHAL}.geno
NEA_SNP=${EIGENSTRAT}/${TAG}${NEANDERTHAL}.snp
NEA_IND=${EIGENSTRAT}/${TAG}${NEANDERTHAL}.ind

ANC=${EIGENSTRAT}/${TAG}${POP}${INTROGRESS}.anc

# --- step 1: filter for shared derived alleles, return numpy array ---------
echo "Step 1: filtering derived alleles shared between ${POP} and ${NEANDERTHAL}"
python ${SCRIPTS}/filter_derived_alleles.py \
    --african-geno     ${AFR_GENO} \
    --african-snp      ${AFR_SNP} \
    --african-ind      ${AFR_IND} \
    --neanderthal-geno ${NEA_GENO} \
    --neanderthal-snp  ${NEA_SNP} \
    --neanderthal-ind  ${NEA_IND} \
    --anc              ${ANC} \
    --save-npy         ${OUT}/filtered 2>&1

echo "Step 2: computing CSFS"
python ${SCRIPTS}/compute_csfs.py \
    --african-geno     ${AFR_GENO} \
    --african-snp      ${AFR_SNP} \
    --african-ind      ${AFR_IND} \
    --neanderthal-geno ${NEA_GENO} \
    --neanderthal-snp  ${NEA_SNP} \
    --neanderthal-ind  ${NEA_IND} \
    --anc              ${ANC} \
    --save-npy         ${OUT}/csfs 2>&1

echo "Step 3: plotting CSFS"
python ${SCRIPTS}/plot_csfs.py \
    --csfs  ${OUT}/csfs_csfs.npy \
    --out   ${OUT}/csfs_plot.pdf \
    --title "CSFS: ${POP} conditioned on ${NEANDERTHAL}" \
    --label "${POP} | ${NEANDERTHAL} derived" 2>&1

echo "Done. Output files:"
echo "  ${OUT}/filtered_geno.npy   — filtered genotype array (SNPs x individuals)"
echo "  ${OUT}/filtered_anc.npy    — true introgression matrix (SNPs x haplotypes)"
echo "  ${OUT}/csfs_csfs.npy       — CSFS vector (length n_chrom-1)"
echo "  ${OUT}/csfs_plot.pdf       — CSFS histogram"
