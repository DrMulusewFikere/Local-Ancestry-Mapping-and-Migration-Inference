#!/bin/bash

# ============================================================
# 02_run_treemix_models.sh
#
# Run TreeMix with multiple migration-edge settings.
#
# Migration models:
#   m = 3, 4, 10, 15, 20, 30
#
# TreeMix:
#   -k 500
#
# The actual number of migration events inferred by TreeMix
# is determined from the .llik output and may be lower than
# the specified -m value.
# ============================================================

#SBATCH -A bharpur
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=60:00:00
#SBATCH --job-name=treemix_models
#SBATCH -o treemix_models.out
#SBATCH -e treemix_models.err


# ------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------

WORKDIR="/depot/bharpur/data/projects/fikere/rfmix/result/flow"

cd "${WORKDIR}" || exit 1


# ------------------------------------------------------------
# TreeMix executable
# ------------------------------------------------------------

TREEMIX="${WORKDIR}/treemix-1.13/src/treemix"


# ------------------------------------------------------------
# Input
# ------------------------------------------------------------

INPUT="treemix_input.txt.gz"


# ------------------------------------------------------------
# TreeMix settings
# ------------------------------------------------------------

K=500

# Migration-edge models to evaluate
MIGRATION_MODELS=(3 4 10 15 20 30)


# ------------------------------------------------------------
# Check input
# ------------------------------------------------------------

if [[ ! -f "${INPUT}" ]]; then
    echo "ERROR: TreeMix input not found:"
    echo "${INPUT}"
    exit 1
fi

if [[ ! -x "${TREEMIX}" ]]; then
    echo "ERROR: TreeMix executable not found:"
    echo "${TREEMIX}"
    exit 1
fi


# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

mkdir -p treemix_models


# ------------------------------------------------------------
# Run TreeMix models
# ------------------------------------------------------------

for M in "${MIGRATION_MODELS[@]}"
do

    echo ""
    echo "============================================================"
    echo "Running TreeMix: m = ${M}"
    echo "============================================================"

    "${TREEMIX}" \
        -i "${INPUT}" \
        -o "treemix_models/output_m${M}" \
        -k "${K}" \
        -m "${M}"

    echo ""
    echo "Completed m = ${M}"

done


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "TreeMix model runs completed"
echo "============================================================"

ls -lh treemix_models/*.llik

echo ""
echo "Migration models evaluated:"
printf '%s\n' "${MIGRATION_MODELS[@]}"

echo ""
echo "============================================================"
