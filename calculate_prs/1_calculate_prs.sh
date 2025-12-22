#!/bin/bash

################################################################################
# POLYGENIC RISK SCORE (PRS) CALCULATION PIPELINE
#
# This script implements best practices for PRS calculation using PLINK2
# Following guidelines from:
# - Choi et al. (2020) "A guide to performing PRS analyses"
# - Lambert et al. (2022) "UK Biobank PRS practical guide"
# - PRS Tutorial (choishingwan.github.io/PRS-Tutorial/)
#
# Pipeline steps:
# 1. Input validation
# 2. Base data (GWAS) QC
# 3. Target data QC
# 4. Allele harmonization
# 5. LD clumping
# 6. PRS calculation across multiple p-value thresholds
# 7. Results summary and QC report
#
# Author: Mantou
# Date: 2025-12-17
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Catch errors in pipes

# ============================================================================
# ERROR HANDLING FUNCTION (trap set up later after LOG_FILE is defined)
# ============================================================================

# Function to handle errors and provide detailed context
error_handler() {
    local line_number=$1
    local command=$2
    local error_code=$3
    
    echo ""
    echo "========================================"
    echo "ERROR: Script failed!"
    echo "========================================"
    echo "Line number: $line_number"
    echo "Command: $command"
    echo "Exit code: $error_code"
    echo "Time: $(date)"
    echo "========================================"
    echo ""
    
    # Show context (5 lines before and after the error)
    if [[ -f "${BASH_SOURCE[1]}" ]]; then
        echo "Context around line $line_number:"
        echo "----------------------------------------"
        awk -v line="$line_number" '
            NR >= line-5 && NR <= line+5 {
                prefix = (NR == line) ? ">>> " : "    "
                printf "%s%4d: %s\n", prefix, NR, $0
            }
        ' "${BASH_SOURCE[1]}"
        echo "----------------------------------------"
        echo ""
    fi
    
    if [[ -n "${LOG_FILE:-}" ]]; then
        echo "Check the log file for more details: $LOG_FILE"
    fi
    echo ""
    
    exit "$error_code"
}

# Function to handle successful completion
success_handler() {
    echo ""
    echo "========================================"
    echo "✓ Script completed successfully!"
    echo "========================================"
    echo ""
}

# ============================================================================
# LOAD CONFIGURATION
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${SCRIPT_DIR}/prs_config.sh"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE"
    echo "Please create prs_config.sh first."
    exit 1
fi

source "$CONFIG_FILE"

# ============================================================================
# SETUP
# ============================================================================

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Setup logging
LOG_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "$LOG_FILE")
exec 2>&1

# Now that LOG_FILE is set, activate error trap
trap 'error_handler ${LINENO} "$BASH_COMMAND" $?' ERR

echo "=================================="
echo "PRS CALCULATION PIPELINE"
echo "=================================="
echo "Project: $PROJECT_NAME"
echo "Started: $(date)"
echo "Output directory: $OUTPUT_DIR"
echo "Using $THREADS threads"
echo "=================================="
echo ""

# ============================================================================
# CHECK DEPENDENCIES
# ============================================================================

echo "[$(date +%H:%M:%S)] Checking dependencies..."

command -v plink2 >/dev/null 2>&1 || {
    echo "ERROR: plink2 not found. Please install PLINK2."
    exit 1
}

PLINK2_VERSION=$(plink2 --version | head -1)
echo "Found: $PLINK2_VERSION"

if command -v Rscript >/dev/null 2>&1; then
    echo "Found: R (for optional analyses)"
else
    echo "WARNING: R not found. Some optional analyses will be skipped."
fi

echo ""

# ============================================================================
# FUNCTION: INPUT VALIDATION
# ============================================================================

validate_inputs() {
    echo "[$(date +%H:%M:%S)] Validating input files..."
    
    local errors=0
    
    # Check target genotype files
    if [[ "$TARGET_FORMAT" == "bed" ]]; then
        for ext in bed bim fam; do
            if [[ ! -f "${TARGET_GENO}.${ext}" ]]; then
                echo "ERROR: Missing ${TARGET_GENO}.${ext}"
                ((errors++))
            fi
        done
    elif [[ "$TARGET_FORMAT" == "pgen" ]]; then
        for ext in pgen pvar psam; do
            if [[ ! -f "${TARGET_GENO}.${ext}" ]]; then
                echo "ERROR: Missing ${TARGET_GENO}.${ext}"
                ((errors++))
            fi
        done
    else
        echo "ERROR: TARGET_FORMAT must be 'bed' or 'pgen'"
        ((errors++))
    fi
    
    # Check GWAS summary statistics
    if [[ ! -f "$GWAS_SUMSTATS" ]]; then
        echo "ERROR: GWAS summary statistics not found: $GWAS_SUMSTATS"
        ((errors++))
    else
        echo "GWAS file: $GWAS_SUMSTATS"
        echo "Preview (first 5 lines):"
        head -5 "$GWAS_SUMSTATS"
    fi
    
    if [[ $errors -gt 0 ]]; then
        echo "ERROR: Found $errors input validation errors. Please fix and retry."
        exit 1
    fi
    
    echo "✓ All input files validated"
    echo ""
}

# ============================================================================
# FUNCTION: BASE DATA QC (GWAS SUMMARY STATISTICS)
# ============================================================================

qc_base_data() {
    echo "[$(date +%H:%M:%S)] Quality control of base data (GWAS)..."
    
    local input="$GWAS_SUMSTATS"
    local output="${PROJECT_NAME}.QC"
    
    # Create R script for base data QC
    cat > qc_base.R << 'RSCRIPT'
#!/usr/bin/env Rscript

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
has_header <- as.logical(args[3])
use_or <- as.logical(args[4])
snp_col <- as.integer(args[5])
a1_col <- as.integer(args[6])
beta_col <- as.integer(args[7])
pval_col <- as.integer(args[8])
base_maf <- as.numeric(args[9])
remove_ambiguous <- as.logical(args[10])
remove_duplicates <- as.logical(args[11])

cat("Reading GWAS summary statistics...\n")

# Count and skip comment lines (lines starting with #)
comment_lines <- 0
con <- file(input_file, "r")
while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (!grepl("^#", line)) {
        break
    }
    comment_lines <- comment_lines + 1
}
close(con)

if (comment_lines > 0) {
    cat(sprintf("Found %d comment lines starting with '#' - skipping these\n", comment_lines))
}

# Read data, skipping comment lines
if (has_header) {
    dat <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, 
                     comment.char = "#", check.names = FALSE)
} else {
    dat <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE,
                     comment.char = "#", skip = comment_lines)
    # Add generic column names
    colnames(dat) <- paste0("V", 1:ncol(dat))
}

cat(sprintf("Loaded %d variants\n", nrow(dat)))

# Get column names
cols <- colnames(dat)
snp_name <- cols[snp_col]
a1_name <- cols[a1_col]
beta_name <- cols[beta_col]
pval_name <- cols[pval_col]

cat(sprintf("SNP column: %s\n", snp_name))
cat(sprintf("Effect allele column: %s\n", a1_name))
cat(sprintf("Effect size column: %s\n", beta_name))
cat(sprintf("P-value column: %s\n", pval_name))

# Rename for easier handling
dat$SNP <- dat[[snp_name]]
dat$A1 <- dat[[a1_name]]
dat$BETA <- dat[[beta_name]]
dat$P <- dat[[pval_name]]

# If using OR, convert to log(OR)
if (use_or) {
    cat("Converting OR to log(OR)...\n")
    dat$BETA <- log(dat$BETA)
}

n_original <- nrow(dat)

# Remove missing values
cat("Removing variants with missing data...\n")
dat <- dat[!is.na(dat$SNP) & !is.na(dat$A1) & !is.na(dat$BETA) & !is.na(dat$P), ]
cat(sprintf("Removed %d variants with missing data\n", n_original - nrow(dat)))

# Remove duplicates
if (remove_duplicates) {
    cat("Removing duplicate SNPs...\n")
    n_before <- nrow(dat)
    dat <- dat[!duplicated(dat$SNP), ]
    cat(sprintf("Removed %d duplicate SNPs\n", n_before - nrow(dat)))
}

# Remove ambiguous SNPs (A/T, G/C)
if (remove_ambiguous) {
    cat("Removing ambiguous SNPs (A/T, G/C)...\n")
    # Need A2 column for this - if not available, infer from A1
    # For simplicity, remove based on allele combinations
    if ("A2" %in% colnames(dat)) {
        ambiguous <- (dat$A1 == "A" & dat$A2 == "T") |
                     (dat$A1 == "T" & dat$A2 == "A") |
                     (dat$A1 == "G" & dat$A2 == "C") |
                     (dat$A1 == "C" & dat$A2 == "G")
        n_ambig <- sum(ambiguous)
        dat <- dat[!ambiguous, ]
        cat(sprintf("Removed %d ambiguous SNPs\n", n_ambig))
    } else {
        cat("WARNING: A2 column not found, cannot filter ambiguous SNPs\n")
    }
}

# Filter on MAF if provided in data
if ("MAF" %in% colnames(dat) && base_maf > 0) {
    cat(sprintf("Filtering on MAF > %.4f...\n", base_maf))
    n_before <- nrow(dat)
    dat <- dat[dat$MAF > base_maf, ]
    cat(sprintf("Removed %d variants with MAF < %.4f\n", n_before - nrow(dat), base_maf))
}

# Filter on INFO if provided
if ("INFO" %in% colnames(dat)) {
    base_info <- as.numeric(args[12])
    if (!is.na(base_info) && base_info > 0) {
        cat(sprintf("Filtering on INFO > %.2f...\n", base_info))
        n_before <- nrow(dat)
        dat <- dat[dat$INFO > base_info, ]
        cat(sprintf("Removed %d variants with INFO < %.2f\n", n_before - nrow(dat), base_info))
    }
}

# Remove variants with P-value issues
cat("Removing variants with invalid P-values...\n")
n_before <- nrow(dat)
dat <- dat[dat$P > 0 & dat$P <= 1, ]
cat(sprintf("Removed %d variants with invalid P-values\n", n_before - nrow(dat)))

cat(sprintf("\nFinal QC'd dataset: %d variants\n", nrow(dat)))

# Write output
write.table(dat, output_file, quote = FALSE, row.names = FALSE, sep = "\t")
cat(sprintf("Saved to: %s\n", output_file))

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat(sprintf("Min P-value: %.2e\n", min(dat$P)))
cat(sprintf("Median P-value: %.2e\n", median(dat$P)))
cat(sprintf("P < 5e-8: %d variants (%.1f%%)\n", 
    sum(dat$P < 5e-8), 100 * sum(dat$P < 5e-8) / nrow(dat)))
cat(sprintf("P < 0.05: %d variants (%.1f%%)\n", 
    sum(dat$P < 0.05), 100 * sum(dat$P < 0.05) / nrow(dat)))

RSCRIPT

    # Run R script
    Rscript qc_base.R \
        "$input" \
        "$output" \
        "$GWAS_HAS_HEADER" \
        "$USE_OR" \
        "$SNP_COL" \
        "$A1_COL" \
        "$BETA_COL" \
        "$PVAL_COL" \
        "$BASE_MAF" \
        "$REMOVE_AMBIGUOUS" \
        "$REMOVE_DUPLICATES" \
        "$BASE_INFO"
    
    # Clean up temp R script (controlled by KEEP_TEMP_FILES)
    if [[ "$KEEP_TEMP_FILES" == false ]]; then
        rm -f qc_base.R
    fi
    
    echo ""
}

# ============================================================================
# FUNCTION: TARGET DATA QC
# ============================================================================

qc_target_data() {
    echo "[$(date +%H:%M:%S)] Quality control of target data..."
    
    local prefix="${PROJECT_NAME}.target"
    
    # Validate input files exist BEFORE trying to use them
    echo "Validating target genotype files..."
    if [[ "$TARGET_FORMAT" == "bed" ]]; then
        for ext in bed bim fam; do
            if [[ ! -f "${TARGET_GENO}.${ext}" ]]; then
                echo "ERROR: Missing file: ${TARGET_GENO}.${ext}"
                echo "Check TARGET_GENO path in config: $TARGET_GENO"
                exit 1
            fi
        done
        echo "✓ Found .bed/.bim/.fam files"
    else
        for ext in pgen pvar psam; do
            if [[ ! -f "${TARGET_GENO}.${ext}" ]]; then
                echo "ERROR: Missing file: ${TARGET_GENO}.${ext}"
                echo "Check TARGET_GENO path in config: $TARGET_GENO"
                exit 1
            fi
        done
        echo "✓ Found .pgen/.pvar/.psam files"
    fi
    echo ""
    
    # Determine input flag
    local input_flag
    if [[ "$TARGET_FORMAT" == "bed" ]]; then
        input_flag="--bfile $TARGET_GENO"
    else
        input_flag="--pfile $TARGET_GENO"
    fi
    
    # Handle sample subsetting
    local subset_flag=""
    if [[ -n "$KEEP_SAMPLES" ]] && [[ -f "$KEEP_SAMPLES" ]]; then
        echo "Keeping samples from: $KEEP_SAMPLES"
        n_keep=$(wc -l < "$KEEP_SAMPLES")
        echo "Number of samples to keep: $n_keep"
        subset_flag="--keep $KEEP_SAMPLES"
    fi
    
    if [[ -n "$REMOVE_SAMPLES" ]] && [[ -f "$REMOVE_SAMPLES" ]]; then
        echo "Removing samples from: $REMOVE_SAMPLES"
        n_remove=$(wc -l < "$REMOVE_SAMPLES")
        echo "Number of samples to remove: $n_remove"
        subset_flag="$subset_flag --remove $REMOVE_SAMPLES"
    fi
    
    echo "Step 1: Sample subset + Variant missingness filter..."
    plink2 \
        $input_flag \
        $subset_flag \
        --geno $TARGET_GENO_MISSING \
        --make-bed \
        --out ${prefix}.geno \
        --threads $THREADS
    
    if [[ ! -f ${prefix}.geno.bed ]]; then
        echo "ERROR: Step 1 failed - no output files created"
        echo "Check if:"
        echo "  - Genotype files exist at: $TARGET_GENO"
        echo "  - KEEP_SAMPLES file is valid (if used): $KEEP_SAMPLES"
        echo "  - PLINK2 log: ${prefix}.geno.log"
        exit 1
    fi
    
    echo "Step 2: Individual missingness filter..."
    plink2 \
        --bfile ${prefix}.geno \
        --mind $TARGET_IND_MISSING \
        --make-bed \
        --out ${prefix}.mind \
        --threads $THREADS
    
    if [[ ! -f ${prefix}.mind.bed ]]; then
        echo "ERROR: Step 2 failed - no samples passed missingness filter"
        echo "Check: TARGET_IND_MISSING threshold (currently: $TARGET_IND_MISSING)"
        exit 1
    fi
    
    echo "Step 3: MAF filter..."
    plink2 \
        --bfile ${prefix}.mind \
        --maf $TARGET_MAF \
        --make-bed \
        --out ${prefix}.maf \
        --threads $THREADS
    
    if [[ ! -f ${prefix}.maf.bed ]]; then
        echo "ERROR: Step 3 failed - no variants passed MAF filter"
        echo "Check: TARGET_MAF threshold (currently: $TARGET_MAF)"
        exit 1
    fi
    
    echo "Step 4: Hardy-Weinberg equilibrium filter..."
    local hwe_flag="--hwe $TARGET_HWE"
    if [[ "$HWE_KEEP_FEWHET" == true ]]; then
        hwe_flag="$hwe_flag keep-fewhet"
    fi
    
    plink2 \
        --bfile ${prefix}.maf \
        $hwe_flag \
        --make-bed \
        --out ${prefix}.qc \
        --threads $THREADS
    
    # Optional: Sex check
    if [[ "$SEX_CHECK" == true ]]; then
        echo "Step 5: Sex check..."
        plink2 \
            --bfile ${prefix}.qc \
            --check-sex \
            --out ${prefix}.sexcheck \
            --threads $THREADS
        
        # Remove samples with sex discordance
        awk '$5=="PROBLEM" {print $1, $2}' ${prefix}.sexcheck.sexcheck > ${prefix}.sex_discordance.txt
        
        if [[ -s ${prefix}.sex_discordance.txt ]]; then
            echo "Found $(wc -l < ${prefix}.sex_discordance.txt) samples with sex discordance"
            plink2 \
                --bfile ${prefix}.qc \
                --remove ${prefix}.sex_discordance.txt \
                --make-bed \
                --out ${prefix}.qc.sex \
                --threads $THREADS
            mv ${prefix}.qc.sex.bed ${prefix}.qc.bed
            mv ${prefix}.qc.sex.bim ${prefix}.qc.bim
            mv ${prefix}.qc.sex.fam ${prefix}.qc.fam
        fi
    fi
    
    # Optional: Relatedness check
    if [[ "$RELATEDNESS_CHECK" == true ]]; then
        echo "Step 6: Relatedness check..."
        plink2 \
            --bfile ${prefix}.qc \
            --king-cutoff $RELATEDNESS_THRESHOLD \
            --out ${prefix}.king \
            --threads $THREADS
        
        if [[ -f ${prefix}.king.king.cutoff.out.id ]]; then
            plink2 \
                --bfile ${prefix}.qc \
                --remove ${prefix}.king.king.cutoff.out.id \
                --make-bed \
                --out ${prefix}.qc.unrelated \
                --threads $THREADS
            mv ${prefix}.qc.unrelated.bed ${prefix}.qc.bed
            mv ${prefix}.qc.unrelated.bim ${prefix}.qc.bim
            mv ${prefix}.qc.unrelated.fam ${prefix}.qc.fam
        fi
    fi
    
    # Generate summary statistics
    echo ""
    echo "=== TARGET DATA QC SUMMARY ==="
    
    # Robust summary based on .bim/.fam, with logs only as optional context
    for step in geno mind maf qc; do
        bim_file="${prefix}.${step}.bim"
        fam_file="${prefix}.${step}.fam"
        log_file="${prefix}.${step}.log"
        
        if [[ -f "$bim_file" && -f "$fam_file" ]]; then
            n_variants=$(wc -l < "$bim_file" 2>/dev/null || echo "?")
            n_samples=$(wc -l < "$fam_file" 2>/dev/null || echo "?")
            echo "After ${step}: ${n_variants} variants, ${n_samples} samples (from files)"
        elif [[ -f "$log_file" ]]; then
            echo "After ${step}: .bed/.bim/.fam missing, but log exists – this step may have failed or produced no output"
        else
            echo "After ${step}: no log and no output files – this step seems to have FAILED QC or not run"
        fi
    done
    echo ""
    
    # Clean up intermediate genotype files (controlled by KEEP_INTERMEDIATE_GENOTYPES)
    if [[ "$KEEP_INTERMEDIATE_GENOTYPES" == false ]]; then
        echo "Cleaning up intermediate genotype files..."
        for step in geno mind maf; do
            if [[ -f ${prefix}.${step}.bed ]]; then
                rm -f ${prefix}.${step}.bed ${prefix}.${step}.bim ${prefix}.${step}.fam
                echo "  Removed: ${prefix}.${step}.{bed,bim,fam}"
            fi
        done
        echo "  Kept: ${prefix}.qc.{bed,bim,fam} (final QC'd genotypes)"
        echo ""
    fi
}

# ============================================================================
# FUNCTION: LD CLUMPING
# ============================================================================

perform_clumping() {
    echo "[$(date +%H:%M:%S)] Performing LD clumping..."
    
    local base_data="${PROJECT_NAME}.QC"
    local target_data="${PROJECT_NAME}.target.qc"
    local output="${PROJECT_NAME}.clumped"
    
    # Determine reference for LD calculation
    local ld_ref
    if [[ -n "$LD_REFERENCE" ]] && [[ -f "${LD_REFERENCE}.bed" ]]; then
        ld_ref="--bfile $LD_REFERENCE"
        echo "Using custom LD reference: $LD_REFERENCE"
    else
        ld_ref="--bfile $target_data"
        echo "Using target data for LD reference"
    fi
    
    # Perform clumping
    plink2 \
        $ld_ref \
        --clump $base_data \
        --clump-p1 $CLUMP_P1 \
        --clump-p2 $CLUMP_P2 \
        --clump-r2 $CLUMP_R2 \
        --clump-kb $CLUMP_KB \
        --clump-snp-field SNP \
        --clump-field P \
        --out $output \
        --threads $THREADS
    
    # Extract valid SNP IDs
    if [[ -f "${output}.clumps" ]]; then
        # PLINK2 format
        awk 'NR!=1 {print $3}' ${output}.clumps > ${output}.valid.snp
    elif [[ -f "${output}.clumped" ]]; then
        # PLINK1.9 format
        awk 'NR!=1 {print $3}' ${output}.clumped > ${output}.valid.snp
    else
        echo "ERROR: Clumping output not found"
        exit 1
    fi
    
    n_clumped=$(wc -l < ${output}.valid.snp)
    echo "Retained $n_clumped independent SNPs after clumping"
    
    # Report on missing SNPs (SNPs in GWAS but not in target genotypes)
    if [[ -f "${output}.clumps.missing_id" ]]; then
        n_missing=$(wc -l < ${output}.clumps.missing_id)
        n_gwas=$(wc -l < $base_data)
        ((n_gwas--))  # Remove header
        pct_missing=$(awk -v m=$n_missing -v t=$n_gwas 'BEGIN {printf "%.1f", (m/t)*100}')
        
        echo ""
        echo "=== SNP OVERLAP SUMMARY ==="
        echo "Total SNPs in GWAS: $n_gwas"
        echo "SNPs missing from target data: $n_missing ($pct_missing%)"
        echo "SNPs available for clumping: $((n_gwas - n_missing))"
        echo "Independent SNPs after clumping: $n_clumped"
        
        if (( n_missing > 0 )); then
            echo ""
            echo "Note: Missing SNPs are listed in: ${output}.clumps.missing_id"
            echo "This is normal - SNPs may be filtered during QC or unavailable in target data."
        fi
    fi
    
    echo ""
}

# ============================================================================
# FUNCTION: CALCULATE PRS
# ============================================================================

calculate_prs() {
    echo "[$(date +%H:%M:%S)] Calculating polygenic risk scores..."
    
    local base_data="${PROJECT_NAME}.QC"
    local target_data="${PROJECT_NAME}.target.qc"
    local valid_snps="${PROJECT_NAME}.clumped.valid.snp"
    
    # Create range file for p-value thresholds
    # Format: range_name lower_bound upper_bound
    local range_file="${PROJECT_NAME}.range_list"
    > $range_file  # Clear file
    
    echo "P-value thresholds:"
    for threshold in "${PRS_THRESHOLDS[@]}"; do
        # Parse: "max_p min_p max_p" format from config
        read -r max_p min_p max_p2 <<< "$threshold"
        
        # Create a clean range name (replace dots and dashes for filename safety)
        range_name=$(echo "$max_p" | sed 's/\./_/g; s/e-/e_/g; s/-/_/g')
        
        # Write to range file: range_name min_p max_p
        echo "$range_name $min_p $max_p" >> $range_file
        
        echo "  P < $max_p → ${PROJECT_NAME}.prs.${range_name}.sscore"
    done
    
    echo ""
    echo "Range file contents:"
    cat $range_file
    echo ""
    
    # Create SNP-pvalue file (using the standardized SNP and P columns from QC)
    local snp_pval="${PROJECT_NAME}.SNP.pvalue"
    echo "Creating SNP-pvalue file..."
    
    # The QC'd file has columns: SNP, A1, BETA, P (plus others)
    # Extract SNP (col 1 after QC) and P (col 4 after QC)
    # But safer to use awk with column names
    awk 'NR==1 {
        for(i=1; i<=NF; i++) {
            if($i=="SNP") snp_col=i;
            if($i=="P") p_col=i;
        }
        next;
    }
    {
        print $snp_col, $p_col
    }' $base_data > $snp_pval
    
    n_snp_pval=$(wc -l < $snp_pval)
    echo "SNP-pvalue file contains $n_snp_pval SNPs"
    echo "Preview:"
    head -3 $snp_pval
    echo ""
    
    # Extract needed columns: SNP, A1, BETA with header
    echo "Creating score file with SNP, effect allele, and effect size..."
    awk 'NR==1 {
        for(i=1; i<=NF; i++) {
            if($i=="SNP") snp_col=i;
            if($i=="A1") a1_col=i;
            if($i=="BETA") beta_col=i;
        }
        print "SNP", "A1", "BETA";
        next;
    }
    {
        print $snp_col, $a1_col, $beta_col
    }' $base_data > ${PROJECT_NAME}.score_file.txt
    
    # Verify score file
    n_snps=$(wc -l < ${PROJECT_NAME}.score_file.txt)
    ((n_snps--))  # Remove header
    echo "Score file contains $n_snps SNPs"
    echo "Preview:"
    head -3 ${PROJECT_NAME}.score_file.txt
    echo ""
    
    # Build score command
    if [[ "$SCORE_METHOD" == "sum" ]]; then
        score_flags="--score ${PROJECT_NAME}.score_file.txt 1 2 3 header sum"
    else
        score_flags="--score ${PROJECT_NAME}.score_file.txt 1 2 3 header"
    fi
    
    if [[ "$MISSING_GENOTYPE" == "no_mean_impute" ]]; then
        score_flags="$score_flags no-mean-imputation"
    fi
    
    # Calculate PRS
    echo "Running PLINK2 --score..."
    plink2 \
        --bfile $target_data \
        $score_flags \
        --q-score-range $range_file $snp_pval \
        --extract $valid_snps \
        --out ${PROJECT_NAME}.prs \
        --threads $THREADS
    
    echo ""
    echo "PRS calculation complete!"
    echo ""
    
    # List output files (with more informative error handling)
    if ls ${PROJECT_NAME}.prs.*.sscore 1> /dev/null 2>&1; then
        echo "=== Output PRS files ==="
        ls -lh ${PROJECT_NAME}.prs.*.sscore
        echo ""
        echo "Total files created: $(ls ${PROJECT_NAME}.prs.*.sscore | wc -l)"
        echo ""
        
        # Show sample from first file
        first_file=$(ls ${PROJECT_NAME}.prs.*.sscore | head -1)
        echo "Preview of $first_file:"
        head -3 "$first_file"
        echo "... (showing first 3 rows)"
    else
        echo "WARNING: No .sscore files were created!"
        echo ""
        echo "Checking for alternative output patterns..."
        ls -lh ${PROJECT_NAME}.prs.* || echo "No output files found at all!"
        echo ""
        echo "Please check the PLINK2 log above for errors."
    fi
    echo ""
}

# ============================================================================
# FUNCTION: GENERATE QC REPORT
# ============================================================================

generate_report() {
    echo "[$(date +%H:%M:%S)] Generating QC report..."
    
    local report="${PROJECT_NAME}_QC_REPORT.txt"
    
    cat > $report << REPORT
================================================================================
POLYGENIC RISK SCORE CALCULATION - QC REPORT
================================================================================

Project: $PROJECT_NAME
Date: $(date)
Output Directory: $OUTPUT_DIR

================================================================================
CONFIGURATION
================================================================================

Target Genotype Data: $TARGET_GENO
GWAS Summary Stats: $GWAS_SUMSTATS

QC Thresholds:
  - Target genotype missingness: $TARGET_GENO_MISSING
  - Target individual missingness: $TARGET_IND_MISSING
  - Target MAF: $TARGET_MAF
  - Target HWE p-value: $TARGET_HWE
  - Base MAF: $BASE_MAF
  - Base INFO: $BASE_INFO

LD Clumping:
  - P1 threshold: $CLUMP_P1
  - P2 threshold: $CLUMP_P2
  - R² threshold: $CLUMP_R2
  - Window size: $CLUMP_KB kb

================================================================================
RESULTS SUMMARY
================================================================================

REPORT

    # Add base data QC summary
    if [[ -f "${PROJECT_NAME}.QC" ]]; then
        n_base=$(wc -l < ${PROJECT_NAME}.QC)
        ((n_base--))  # Remove header
        echo "Base Data (GWAS):" >> $report
        echo "  - Variants after QC: $n_base" >> $report
        echo "" >> $report
    fi
    
    # Add target data QC summary (from .bim/.fam, not log text)
    if [[ -f "${PROJECT_NAME}.target.qc.bim" && -f "${PROJECT_NAME}.target.qc.fam" ]]; then
        n_variants=$(wc -l < ${PROJECT_NAME}.target.qc.bim)
        n_samples=$(wc -l < ${PROJECT_NAME}.target.qc.fam)
        echo "Target Data:" >> $report
        echo "  - Final: $n_variants variants, $n_samples samples" >> $report
        echo "" >> $report
    elif [[ -f "${PROJECT_NAME}.target.qc.log" ]]; then
        echo "Target Data:" >> $report
        echo "  - WARNING: .bim/.fam files missing; target QC may have failed (see log ${PROJECT_NAME}.target.qc.log)" >> $report
        echo "" >> $report
    fi
    
    # Add clumping summary
    if [[ -f "${PROJECT_NAME}.clumped.valid.snp" ]]; then
        n_clumped=$(wc -l < ${PROJECT_NAME}.clumped.valid.snp)
        echo "LD Clumping:" >> $report
        echo "  - Independent SNPs retained: $n_clumped" >> $report
        
        # Add missing SNPs information if available
        if [[ -f "${PROJECT_NAME}.clumped.clumps.missing_id" ]]; then
            n_missing=$(wc -l < ${PROJECT_NAME}.clumped.clumps.missing_id)
            n_base=$(wc -l < ${PROJECT_NAME}.QC)
            ((n_base--))  # Remove header
            pct_missing=$(awk -v m=$n_missing -v t=$n_base 'BEGIN {printf "%.1f", (m/t)*100}')
            
            echo "  - SNPs in GWAS but missing from target: $n_missing ($pct_missing%)" >> $report
            echo "  - SNPs available for clumping: $((n_base - n_missing))" >> $report
        fi
        
        echo "" >> $report
    fi
    
    # Add PRS calculation summary
    echo "PRS Calculation:" >> $report
    echo "  - Method: $SCORE_METHOD" >> $report
    echo "  - P-value thresholds: ${#PRS_THRESHOLDS[@]}" >> $report
    for prs_file in ${PROJECT_NAME}.prs.*.sscore; do
        if [[ -f "$prs_file" ]]; then
            threshold=$(basename $prs_file .sscore | sed 's/.*prs\.//')
            n_samples=$(wc -l < $prs_file)
            ((n_samples--))  # Remove header
            echo "    * $threshold: $n_samples samples" >> $report
        fi
    done
    
    cat >> $report << REPORT

================================================================================
OUTPUT FILES
================================================================================

PRS Scores (*.sscore files):
  - Each file contains PRS for one p-value threshold
  - Columns: #IID, ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG

QC Files:
  - ${PROJECT_NAME}.QC: QC'd GWAS summary statistics
  - ${PROJECT_NAME}.target.qc.*: QC'd target genotype data
  - ${PROJECT_NAME}.clumped.valid.snp: List of independent SNPs

================================================================================
NEXT STEPS
================================================================================

1. Choose the best-fit PRS threshold:
   - Test association with phenotype
   - Select threshold with highest R²/AUC

2. Adjust for population stratification:
   - Calculate principal components (if not done)
   - Include PCs as covariates in analysis

3. Evaluate PRS performance:
   - Calculate variance explained (R²)
   - Generate percentile ranks
   - Assess predictive accuracy

================================================================================
REPORT

    echo "Report saved to: $report"
    cat $report
    echo ""
}

# ============================================================================
# FUNCTION: CLEANUP FILES
# ============================================================================

cleanup_files() {
    echo "[$(date +%H:%M:%S)] Cleaning up files based on retention settings..."
    echo ""
    
    local files_removed=0
    local categories_cleaned=()
    
    # Helper function to safely increment counter (avoids 'set -e' issues)
    increment_counter() {
        files_removed=$((files_removed + 1))
    }
    
    # ========================================
    # 1. QC'd data files
    # ========================================
    if [[ "$KEEP_QC_FILES" == false ]]; then
        echo "Removing QC'd data files..."
        
        # Base QC file
        if [[ -f "${PROJECT_NAME}.QC" ]]; then
            rm -f "${PROJECT_NAME}.QC"
            increment_counter
            echo "  ✓ Removed: ${PROJECT_NAME}.QC"
        fi
        
        # Target QC genotype files
        if [[ -f "${PROJECT_NAME}.target.qc.bed" ]]; then
            rm -f "${PROJECT_NAME}.target.qc.bed" \
                  "${PROJECT_NAME}.target.qc.bim" \
                  "${PROJECT_NAME}.target.qc.fam"
            files_removed=$((files_removed + 3))
            echo "  ✓ Removed: ${PROJECT_NAME}.target.qc.{bed,bim,fam}"
        fi
        
        categories_cleaned+=("QC'd data files")
    else
        echo "Keeping QC'd data files (KEEP_QC_FILES=true)"
    fi
    echo ""
    
    # ========================================
    # 2. Clumping files
    # ========================================
    if [[ "$KEEP_CLUMPING_FILES" == false ]]; then
        echo "Removing clumping files..."
        
        # All clumping outputs
        local clump_files=(
            "${PROJECT_NAME}.clumped.clumps"
            "${PROJECT_NAME}.clumped.clumped"
            "${PROJECT_NAME}.clumped.valid.snp"
            "${PROJECT_NAME}.clumped.clumps.missing_id"
        )
        
        for file in "${clump_files[@]}"; do
            if [[ -f "$file" ]]; then
                rm -f "$file"
                increment_counter
                echo "  ✓ Removed: $(basename $file)"
            fi
        done
        
        categories_cleaned+=("clumping files")
    else
        echo "Keeping clumping files (KEEP_CLUMPING_FILES=true)"
    fi
    echo ""
    
    # ========================================
    # 3. Temporary processing files
    # ========================================
    if [[ "$KEEP_TEMP_FILES" == false ]]; then
        echo "Removing temporary processing files..."
        
        local temp_files=(
            "qc_base.R"
            "${PROJECT_NAME}.score_file.txt"
            "${PROJECT_NAME}.SNP.pvalue"
            "${PROJECT_NAME}.range_list"
            "${PROJECT_NAME}.target.sex_discordance.txt"
            "${PROJECT_NAME}.target.sexcheck.sexcheck"
            "${PROJECT_NAME}.target.king.king.cutoff.in.id"
            "${PROJECT_NAME}.target.king.king.cutoff.out.id"
        )
        
        for file in "${temp_files[@]}"; do
            if [[ -f "$file" ]]; then
                rm -f "$file"
                increment_counter
                echo "  ✓ Removed: $(basename $file)"
            fi
        done
        
        categories_cleaned+=("temp files")
    else
        echo "Keeping temporary files (KEEP_TEMP_FILES=true)"
    fi
    echo ""
    
    # ========================================
    # 4. Log files
    # ========================================
    if [[ "$KEEP_LOGS" == false ]]; then
        echo "Removing log files..."
        
        # Find and remove all .log files EXCEPT the main pipeline log
        local main_log=$(basename "$LOG_FILE")
        
        while IFS= read -r -d '' logfile; do
            local logname=$(basename "$logfile")
            if [[ "$logname" != "$main_log" ]]; then
                rm -f "$logfile"
                increment_counter
                echo "  ✓ Removed: $logname"
            fi
        done < <(find . -maxdepth 1 -name "*.log" -print0)
        
        categories_cleaned+=("log files")
    else
        echo "Keeping log files (KEEP_LOGS=true)"
    fi
    echo ""
    
    # ========================================
    # 5. QC report
    # ========================================
    if [[ "$GENERATE_QC_REPORT" == false ]]; then
        echo "Removing QC report..."
        
        if [[ -f "${PROJECT_NAME}_QC_REPORT.txt" ]]; then
            rm -f "${PROJECT_NAME}_QC_REPORT.txt"
            increment_counter
            echo "  ✓ Removed: ${PROJECT_NAME}_QC_REPORT.txt"
        fi
        
        categories_cleaned+=("QC report")
    fi
    echo ""
    
    # ========================================
    # Summary
    # ========================================
    echo "========================================"
    echo "CLEANUP SUMMARY"
    echo "========================================"
    echo "Total files removed: $files_removed"
    
    if [[ ${#categories_cleaned[@]} -gt 0 ]]; then
        echo "Categories cleaned:"
        for category in "${categories_cleaned[@]}"; do
            echo "  - $category"
        done
    fi
    echo ""
    
    # Show what remains
    echo "Files retained in output directory:"
    echo ""
    
    # Always kept: .sscore files
    if ls ${PROJECT_NAME}.prs.*.sscore 1> /dev/null 2>&1; then
        echo "✓ PRS results (.sscore files):"
        ls -1 ${PROJECT_NAME}.prs.*.sscore | sed 's/^/    /'
    fi
    echo ""
    
    # Conditionally kept items
    if [[ "$KEEP_QC_FILES" == true ]]; then
        echo "✓ QC'd data files:"
        [[ -f "${PROJECT_NAME}.QC" ]] && echo "    ${PROJECT_NAME}.QC"
        [[ -f "${PROJECT_NAME}.target.qc.bed" ]] && echo "    ${PROJECT_NAME}.target.qc.{bed,bim,fam}"
        echo ""
    fi
    
    if [[ "$KEEP_CLUMPING_FILES" == true ]]; then
        echo "✓ Clumping files:"
        ls -1 ${PROJECT_NAME}.clumped.* 2>/dev/null | sed 's/^/    /' || echo "    (none found)"
        echo ""
    fi
    
    if [[ "$KEEP_LOGS" == true ]]; then
        echo "✓ Log files:"
        ls -1 *.log 2>/dev/null | sed 's/^/    /' || echo "    (none found)"
        echo ""
    fi
    
    if [[ "$GENERATE_QC_REPORT" == true ]]; then
        echo "✓ QC Report:"
        [[ -f "${PROJECT_NAME}_QC_REPORT.txt" ]] && echo "    ${PROJECT_NAME}_QC_REPORT.txt"
        echo ""
    fi
    
    echo "========================================"
    echo ""
}

# ============================================================================
# MAIN PIPELINE
# ============================================================================

main() {
    echo "Starting PRS calculation pipeline..."
    echo ""
    
    if [[ "$RUN_VALIDATION" == true ]]; then
        validate_inputs
    fi
    
    qc_base_data
    qc_target_data
    perform_clumping
    calculate_prs
    
    if [[ "$GENERATE_QC_REPORT" == true ]]; then
        generate_report
    fi
    
    # Clean up files based on retention settings
    cleanup_files
    
    echo "=================================="
    echo "PIPELINE COMPLETE!"
    echo "=================================="
    echo "Completed: $(date)"
    echo "Total runtime: $SECONDS seconds"
    echo ""
    echo "Output directory: $OUTPUT_DIR"
    if [[ "$KEEP_LOGS" == true ]]; then
        echo "Main log file: $LOG_FILE"
    fi
    echo ""
    echo "Your PRS scores are in: ${PROJECT_NAME}.prs.*.sscore"
    echo ""
    
    # Mark successful completion
    success_handler
}

# Run main pipeline
main

################################################################################
# END OF SCRIPT
################################################################################