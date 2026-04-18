#!/bin/bash
# ============================================================
# 3 ⚠ outlier seed-variance study (thin_6, thin_9, particulate_12)
# Each case runs 5 independent seeds → 15 total simulations
# Each run creates its own directories:
#   post_<caseID>_s<N>/       atom/mesh/contact dumps
#   restart_<caseID>_s<N>/    restart snapshots
#   plate_<caseID>_s<N>.stl   dynamic plate
#   log_<caseID>_s<N>.out     stdout
# ============================================================

set -u

# --- Configuration ---
# LIGGGHTS launch command — processors auto-distributed (no -processors in script)
# --oversubscribe: allow more MPI ranks than physical cores (needed on WSL/laptops).
# Default: 10 ranks on 10 threads. Adjust if fewer/more threads available.
LIGGGHTS_CMD="mpirun --oversubscribe -np 10 liggghts"

# Skip already-completed runs (check for contact dumps)
SKIP_IF_DONE=true

# Case definitions: (case_id, script, seeds)
# Seeds: first is original, rest are independent primes
declare -A CASE_SCRIPTS=(
    [thin6]="thin6_seed.liggghts"
    [thin9]="thin9_seed.liggghts"
    [part12]="particulate12_seed.liggghts"
)

declare -A CASE_SEEDS=(
    [thin6]="119089 230491 347111 463019 581249"
    [thin9]="42677 156821 273949 391883 509731"
    [part12]="80363 201571 329047 458399 589453"
)

# Execution order (run all seeds of one case before moving to next)
CASE_ORDER=(thin6 thin9 part12)

# --- Runner ---
OVERALL_START=$(date +%s)
TOTAL=0
DONE=0
FAILED=0

for CASE_ID in "${CASE_ORDER[@]}"; do
    SCRIPT="${CASE_SCRIPTS[$CASE_ID]}"
    SEEDS="${CASE_SEEDS[$CASE_ID]}"
    SEED_ARR=($SEEDS)

    if [ ! -f "${SCRIPT}" ]; then
        echo "!! ERROR: ${SCRIPT} not found. Skipping case ${CASE_ID}."
        continue
    fi

    echo ""
    echo "###################################################"
    echo "###  CASE: ${CASE_ID}  (script: ${SCRIPT})"
    echo "###  Seeds: ${SEEDS}"
    echo "###################################################"

    for i in "${!SEED_ARR[@]}"; do
        SEED=${SEED_ARR[$i]}
        IDX=$((i+1))
        RUN_ID="${CASE_ID}_s${IDX}"
        POST_DIR="post_${RUN_ID}"
        LOG="log_${RUN_ID}.out"
        TOTAL=$((TOTAL+1))

        echo ""
        echo "--- [${TOTAL}/15] ${CASE_ID} seed ${IDX}/5 = ${SEED} → ${RUN_ID}"
        echo "    Start: $(date '+%Y-%m-%d %H:%M:%S')"

        if [ "${SKIP_IF_DONE}" = "true" ] && ls "${POST_DIR}"/contact_*.liggghts 1>/dev/null 2>&1; then
            echo "    SKIP: ${POST_DIR} already has contact files"
            DONE=$((DONE+1))
            continue
        fi

        START=$(date +%s)
        ${LIGGGHTS_CMD} -in "${SCRIPT}" -var seed "${SEED}" -var run_id "${RUN_ID}" 2>&1 | tee "${LOG}"
        STATUS=${PIPESTATUS[0]}
        END=$(date +%s)
        DUR=$((END - START))
        H=$((DUR / 3600))
        M=$(((DUR % 3600) / 60))

        if [ ${STATUS} -eq 0 ]; then
            echo "    DONE (${H}h ${M}m)"
            DONE=$((DONE+1))
        else
            echo "    FAILED exit=${STATUS}  (${H}h ${M}m)  — check ${LOG}"
            FAILED=$((FAILED+1))
        fi
    done
done

# --- Summary ---
OVERALL_END=$(date +%s)
OVERALL_DUR=$((OVERALL_END - OVERALL_START))
OH=$((OVERALL_DUR / 3600))
OM=$(((OVERALL_DUR % 3600) / 60))

echo ""
echo "###################################################"
echo "###  ALL 3 CASES × 5 SEEDS PROCESSED"
echo "###  Total time: ${OH}h ${OM}m"
echo "###  Done: ${DONE}/${TOTAL}   Failed: ${FAILED}"
echo "###################################################"
echo ""
echo "Output directory status:"
for CASE_ID in "${CASE_ORDER[@]}"; do
    SEED_ARR=(${CASE_SEEDS[$CASE_ID]})
    for i in "${!SEED_ARR[@]}"; do
        IDX=$((i+1))
        RUN_ID="${CASE_ID}_s${IDX}"
        if ls "post_${RUN_ID}"/contact_*.liggghts 1>/dev/null 2>&1; then
            LAST=$(ls -t post_${RUN_ID}/contact_*.liggghts | head -1)
            echo "  ✓ ${RUN_ID}  (seed=${SEED_ARR[$i]})  → ${LAST}"
        else
            echo "  ✗ ${RUN_ID}  (seed=${SEED_ARR[$i]})  INCOMPLETE"
        fi
    done
done
