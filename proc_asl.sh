#!/bin/bash
# Automated ASL Pipeline - Optimized for Single or Multi-Patient Processing


# Run: ./proc_asl_adulte_all.sh SUBJECT_NAME      or      ./proc_asl_adulte_all.sh
# --- Fixed Physical Parameters (3T Standard) ---
PLD=1.8
BOLUS=1.8
TR_M0=10.0
SLOPE=0.1

# --- Global output file ---
GLOBAL_SUMMARY="./resumen_cbf_grupal.txt"

# Initialize global summary file with header only if it doesn't already exist
if [ ! -f "$GLOBAL_SUMMARY" ]; then
    echo "ID,CBF_GM_Mean,CBF_WM_Mean" > "$GLOBAL_SUMMARY"
fi

# --- Fixed input filenames per patient ---
ASL_IN="pcasl3d_singleTI.nii"
M0_SEP="pcasl3d_M0.nii"
T1_RAW="denoised_mp2rage.nii"

# FILTER CONTROL: If an argument is passed (e.g., ./script.sh MEGAN), process only that subject
if [ -n "$1" ]; then
    # Strip trailing slashes if the user included them by mistake
    SINGLE_SUBJECT=$(basename "$1")
    if [ -d "$SINGLE_SUBJECT" ]; then
        FOLDER_LIST="${SINGLE_SUBJECT}/"
        echo "--> Single-patient mode active: Processing [$SINGLE_SUBJECT] only"
    else
        echo "ERROR: Folder '$SINGLE_SUBJECT' does not exist in the current directory."
        exit 1
    fi
else
    FOLDER_LIST="*/"
    echo "--> Multi-patient mode active: Processing entire directory"
fi

# Main processing loop
for folder in $FOLDER_LIST; do
    # Extract clean ID from folder
    ID=$(basename "$folder")
    
    # Ignore previous results folders or system structures
    if [[ "$ID" == results_* ]] || [[ "$ID" == mean_CBF ]] || [[ "$ID" == Reporte_Visual_CBF ]]; then
        continue
    fi

    INPUT_DIR="./${ID}"
    OUT="./results_${ID}"
    
    echo "===================================================="
    echo "CHECKING STRUCTURE FOR PATIENT: $ID"
    echo "===================================================="

    # SMART FILE ASSIGNMENT (Detects if files are .nii or .nii.gz)
    ASL_IN_FILE=$(ls "${INPUT_DIR}"/pcasl3d_singleTI.nii* 2>/dev/null)
    M0_SEP_FILE=$(ls "${INPUT_DIR}"/pcasl3d_M0.nii* 2>/dev/null)
    T1_RAW_FILE=$(ls "${INPUT_DIR}"/denoised_mp2rage.nii* 2>/dev/null)

    # Verify that all 3 key files were found
    if [ -z "$ASL_IN_FILE" ] || [ -z "$M0_SEP_FILE" ] || [ -z "$T1_RAW_FILE" ]; then
        echo "Skipping $ID: Missing pcasl3d_singleTI, pcasl3d_M0, or denoised_mp2rage."
        continue
    fi

    echo "-> All set. Starting processing for $ID..."
    mkdir -p "$OUT"

    # 1. T1 BRAIN EXTRACTION AND SEGMENTATION
    echo "[1/5] Anatomy and segmentation..."
    #bet "$T1_RAW_FILE" "${OUT}/T1_brain.nii.gz" -f 0.2 -B -m
    #fast -t 1 -n 3 -o "${OUT}/T1_seg" "${OUT}/T1_brain.nii.gz"

    # 2. ASL PREPARATION AND REGISTRATIONS
    echo "[2/5] Extracting volumes and computing registrations..."
    fslroi "$ASL_IN_FILE" "${OUT}/asl_pairs.nii.gz" 1 20
    fslroi "$ASL_IN_FILE" "${OUT}/asl_ref.nii.gz" 0 1

    # Register M0 to ASL space
    flirt -in "$M0_SEP_FILE" -ref "${OUT}/asl_ref.nii.gz" -out "${OUT}/m0_reg.nii.gz" -dof 6

    # Register structural PVE maps to ASL space
    flirt -in "${OUT}/T1_brain_mask.nii.gz" -ref "${OUT}/asl_ref.nii.gz" \
          -out "${OUT}/asl_mask_asl.nii.gz" -interp nearestneighbour -dof 6

    flirt -in "${OUT}/T1_seg_pve_1.nii.gz" -ref "${OUT}/asl_ref.nii.gz" \
          -out "${OUT}/GM_pve_asl.nii.gz" -interp trilinear -dof 6

    flirt -in "${OUT}/T1_seg_pve_2.nii.gz" -ref "${OUT}/asl_ref.nii.gz" \
          -out "${OUT}/WM_pve_asl.nii.gz" -interp trilinear -dof 6

    # Binarize with the new strict WM threshold (0.95)
    fslmaths "${OUT}/GM_pve_asl.nii.gz" -thr 0.5 -bin "${OUT}/GM_mask_asl.nii.gz"
    fslmaths "${OUT}/WM_pve_asl.nii.gz" -thr 0.95 -bin "${OUT}/WM_mask_asl.nii.gz"

    # 3. OXFORD_ASL QUANTIFICATION
    echo "[3/5] Running oxford_asl..."
    oxford_asl -i "${OUT}/asl_pairs.nii.gz" \
               -o "${OUT}/oxford_asl" \
               -c "${OUT}/m0_reg.nii.gz" \
               -m "${OUT}/asl_mask_asl.nii.gz" \
               --tr $TR_M0 --plds $PLD --bolus $BOLUS \
               --iaf tc --casl --wp --mc --slicedt 0 \
               --t1b 1.65 --alpha 0.85 --cmethod voxel --bat 1.2

    # 4. GENERATE FINAL QUANTITATIVE MAP (SLOPE=1)
    echo "[4/5] Saving absolute CBF map..."
    RAW_CBF="${OUT}/oxford_asl/native_space/perfusion_calib.nii.gz"
    FINAL_CBF="${OUT}/${ID}_CBF_QUANT.nii.gz"

    if [ -f "$RAW_CBF" ]; then
        fslmaths "$RAW_CBF" -mul $SLOPE "$FINAL_CBF"

        # 5. COMPUTE STATISTICS AND REPORT
        GM_VAL=$(fslstats "$FINAL_CBF" -k "${OUT}/GM_mask_asl.nii.gz" -M)
        WM_VAL=$(fslstats "$FINAL_CBF" -k "${OUT}/WM_mask_asl.nii.gz" -M)
        
        echo "----------------------------------------------------"
        echo "RESULTS FOR $ID:"
        echo "CBF mean GM: $GM_VAL mL/100g/min"
        echo "CBF mean WM: $WM_VAL mL/100g/min"
        echo "----------------------------------------------------"
        
        # Save individual text report in the results folder
        LOG_FILE="${OUT}/${ID}_reporte_final.txt"
        echo "Date: $(date)" > "$LOG_FILE"
        echo "Patient ID: $ID" >> "$LOG_FILE"
        echo "Parameters: PLD=$PLD, Bolus=$BOLUS" >> "$LOG_FILE"
        echo "Applied WM Threshold: 0.95 (PVE)" >> "$LOG_FILE"
        echo "Mean GM CBF: $GM_VAL mL/100g/min" >> "$LOG_FILE"
        echo "Mean WM CBF: $WM_VAL mL/100g/min" >> "$LOG_FILE"

        # Append data to the centralized summary file in CSV format
        echo "${ID},${GM_VAL},${WM_VAL}" >> "$GLOBAL_SUMMARY"
    else
        echo "ERROR: oxford_asl failed to generate maps for $ID."
    fi
    echo ""
done

echo "===================================================="
echo "PROCESS COMPLETED"
echo "Results logged to: $GLOBAL_SUMMARY"
echo "===================================================="
