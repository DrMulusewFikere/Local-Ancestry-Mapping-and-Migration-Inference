#!/bin/bash

# ============================================================
# 03_summarize_model_fit.sh
#
# Summarize TreeMix model fit across migration-edge settings.
#
# Input:
#   treemix_models/output_m*.llik
#
# Output:
#   results/model_selection/treemix_model_comparison.tsv
# ============================================================


# ------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------

WORKDIR="/depot/bharpur/data/projects/fikere/rfmix/result/flow"

cd "${WORKDIR}" || exit 1


# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

mkdir -p results/model_selection


OUTPUT="results/model_selection/treemix_model_comparison.tsv"


# ------------------------------------------------------------
# Header
# ------------------------------------------------------------

printf "Requested_m\tStarting_migration_events\tStarting_logLik\tActual_migration_events\tFinal_logLik\n" \
    > "${OUTPUT}"


# ------------------------------------------------------------
# Process likelihood files
# ------------------------------------------------------------

for FILE in treemix_models/output_m*.llik
do

    if [[ ! -f "${FILE}" ]]; then
        continue
    fi


    # Extract requested m from filename
    BASENAME=$(basename "${FILE}")
    REQUESTED_M=$(echo "${BASENAME}" | sed 's/output_m//' | sed 's/.llik//')


    # Starting migration events
    START_M=$(grep "Starting ln(likelihood)" "${FILE}" |
        sed -E 's/.*with ([0-9]+) migration events:.*/\1/')


    # Starting likelihood
    START_LIK=$(grep "Starting ln(likelihood)" "${FILE}" |
        sed -E 's/.*events: ([-+0-9.eE]+).*/\1/')


    # Actual/final migration events
    FINAL_M=$(grep "Exiting ln(likelihood)" "${FILE}" |
        sed -E 's/.*with ([0-9]+) migration events:.*/\1/')


    # Final likelihood
    FINAL_LIK=$(grep "Exiting ln(likelihood)" "${FILE}" |
        sed -E 's/.*events: ([-+0-9.eE]+).*/\1/')


    # Write result
    printf "%s\t%s\t%s\t%s\t%s\n" \
        "${REQUESTED_M}" \
        "${START_M}" \
        "${START_LIK}" \
        "${FINAL_M}" \
        "${FINAL_LIK}" \
        >> "${OUTPUT}"

done


# ------------------------------------------------------------
# Sort by requested migration setting
# ------------------------------------------------------------

{
    head -n 1 "${OUTPUT}"
    tail -n +2 "${OUTPUT}" |
        sort -n -k1,1
} > "${OUTPUT}.tmp"

mv "${OUTPUT}.tmp" "${OUTPUT}"


# ------------------------------------------------------------
# Display results
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "TreeMix model comparison"
echo "============================================================"

column -t -s $'\t' "${OUTPUT}"

echo ""
echo "Results saved to:"
echo "${OUTPUT}"

echo "============================================================"
