#!/bin/bash

# ============================================================
# 01_prepare_treemix_input.sh
#
# Prepare population allele-frequency data for TreeMix.
#
# Input:
#   - input.vcf.gz
#   - popmap.txt
#
# Output:
#   - biallelic_snps.withID.vcf.gz
#   - input_plink_withID.*
#   - pop_freq.frq.strat
#   - treemix_input.txt.gz
#
# TreeMix analysis:
#   47 geographic populations
#   521,270 SNPs
# ============================================================

#SBATCH -A bharpur
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=30:00:00
#SBATCH --job-name=prepare_treemix
#SBATCH -o prepare_treemix.out
#SBATCH -e prepare_treemix.err


# ------------------------------------------------------------
# Load software
# ------------------------------------------------------------

module load biocontainers
module load bcftools/1.17
module load vcftools
module load gsl


# ------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------

WORKDIR="/depot/bharpur/data/projects/fikere/rfmix/result/flow"

cd "${WORKDIR}" || exit 1


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

VCF="input.vcf.gz"
POPMAP="popmap.txt"

# Path to PLINK
PLINK="/depot/bharpur/data/projects/fikere/plink"

# TreeMix input conversion script
FREQ_SCRIPT="./freq_to_treemix_faster.py"


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if [[ ! -f "${VCF}" ]]; then
    echo "ERROR: VCF file not found: ${VCF}"
    exit 1
fi

if [[ ! -f "${POPMAP}" ]]; then
    echo "ERROR: Population map not found: ${POPMAP}"
    exit 1
fi

if [[ ! -f "${FREQ_SCRIPT}" ]]; then
    echo "ERROR: TreeMix conversion script not found: ${FREQ_SCRIPT}"
    exit 1
fi


# ------------------------------------------------------------
# 1. Retain biallelic SNPs
# ------------------------------------------------------------

echo "Step 1: Extracting biallelic SNPs..."

vcftools \
    --gzvcf "${VCF}" \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --out biallelic_snps


# ------------------------------------------------------------
# 2. Prepare VCF header
# ------------------------------------------------------------

echo "Step 2: Preparing VCF header..."

bcftools view -h \
    biallelic_snps.recode.vcf |
    grep '^##contig' > contigs.txt


bcftools annotate \
    --header-lines contigs.txt \
    biallelic_snps.recode.vcf \
    -Oz \
    -o biallelic_snps.header.vcf.gz


bcftools index \
    biallelic_snps.header.vcf.gz


# ------------------------------------------------------------
# 3. Add unique SNP IDs
# ------------------------------------------------------------

echo "Step 3: Adding SNP IDs..."

bcftools annotate \
    --set-id '%CHROM_%POS' \
    biallelic_snps.header.vcf.gz \
    -Oz \
    -o biallelic_snps.withID.vcf.gz


bcftools index \
    biallelic_snps.withID.vcf.gz


# ------------------------------------------------------------
# 4. Convert VCF to PLINK
# ------------------------------------------------------------

echo "Step 4: Converting VCF to PLINK..."

"${PLINK}" \
    --vcf biallelic_snps.withID.vcf.gz \
    --make-bed \
    --out input_plink_withID


# ------------------------------------------------------------
# 5. Calculate population allele frequencies
# ------------------------------------------------------------

echo "Step 5: Calculating population allele frequencies..."

"${PLINK}" \
    --bfile input_plink_withID \
    --freq \
    --within "${POPMAP}" \
    --out pop_freq


# ------------------------------------------------------------
# 6. Convert allele frequencies to TreeMix format
# ------------------------------------------------------------

echo "Step 6: Creating TreeMix input..."

python "${FREQ_SCRIPT}" \
    pop_freq.frq.strat \
    treemix_input.txt


# ------------------------------------------------------------
# 7. Compress TreeMix input
# ------------------------------------------------------------

gzip -f treemix_input.txt


# ------------------------------------------------------------
# 8. Report SNP/population counts
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "TreeMix input preparation completed"
echo "============================================================"

echo "PLINK SNP count:"
wc -l input_plink_withID.bim

echo "Population count:"
cut -f 3 "${POPMAP}" | sort -u | wc -l

echo "TreeMix input:"
ls -lh treemix_input.txt.gz

echo "============================================================"
