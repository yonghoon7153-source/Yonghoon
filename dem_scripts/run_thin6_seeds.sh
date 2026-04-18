#!/bin/bash
# ============================================================
# thin_6 seed-variance study: run 5 independent seeds sequentially
# Each seed creates separate directories:
#   post_thin6_s1/, restart_thin6_s1/, plate_thin6_s1.stl, log_thin6_s1.out
#   post_thin6_s2/, ... (and so on)
# ============================================================

set -e  # exit on error in the loop control (not in LIGGGHTS itself)

# --- Configuration ---
SCRIPT="thin6_seed.liggghts"
# LIGGGHTS launch command — processors auto-distributed (no -processors in script)
LIGGGHTS_CMD="mpirun -np 10 liggghts"

# 5 independent insertion seeds (primes, max statistical independence)
# First is original thin_6 seed; rest are new.
SEEDS=(119089 230491 347111 463019 581249)

# Skip already-completed runs (checks for contact file)
SKIP_IF_DONE=true

# --- Main loop ---
for i in "${!SEEDS[@]}"; do
    SEED=${SEEDS[$i]}
    IDX=$((i+1))
    RUN_ID="thin6_s${IDX}"
    POST_DIR="post_${RUN_ID}"
    LOG="log_${RUN_ID}.out"

    echo ""
    echo "========================================================"
    echo "  [${IDX}/5] Running thin_6 with seed=${SEED} → ${RUN_ID}"
    echo "  Start: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================================"

    # Skip if already done
    if [ "$SKIP_IF_DONE" = "true" ] && ls "${POST_DIR}"/contact_*.liggghts 1>/dev/null 2>&1; then
        echo "  SKIP: ${POST_DIR} already has contact files"
        continue
    fi

    # Clean any partial outputs from prior aborted run (optional — commented out by default)
    # rm -rf "${POST_DIR}" "restart_${RUN_ID}" "plate_${RUN_ID}.stl"

    START=$(date +%s)
    ${LIGGGHTS_CMD} -in ${SCRIPT} -var seed ${SEED} -var run_id ${RUN_ID} 2>&1 | tee "${LOG}"
    STATUS=${PIPESTATUS[0]}
    END=$(date +%s)
    DURATION=$((END - START))
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))

    if [ $STATUS -eq 0 ]; then
        echo "  [${IDX}/5] DONE seed=${SEED}  (${HOURS}h ${MINUTES}m)"
    else
        echo "  [${IDX}/5] FAILED seed=${SEED} (exit ${STATUS}, ${HOURS}h ${MINUTES}m)"
        echo "  Inspect ${LOG} and either resume from restart_${RUN_ID}/ or re-run."
    fi
done

echo ""
echo "========================================================"
echo "  All seeds processed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================================"
echo ""
echo "Output layout:"
for i in "${!SEEDS[@]}"; do
    IDX=$((i+1))
    RUN_ID="thin6_s${IDX}"
    if ls "post_${RUN_ID}"/contact_*.liggghts 1>/dev/null 2>&1; then
        LAST_CONTACT=$(ls -t post_${RUN_ID}/contact_*.liggghts | head -1)
        echo "  ${RUN_ID}  [seed=${SEEDS[$i]}]  → ${LAST_CONTACT}"
    else
        echo "  ${RUN_ID}  [seed=${SEEDS[$i]}]  → INCOMPLETE"
    fi
done
