#!/bin/bash

################################################################################
# COMPLETE PRS ANALYSIS PIPELINE
#
# This script runs the complete polygenic risk score analysis pipeline:
# 1. Calculate PRS scores from GWAS data (1_calculate_prs.sh)
# 2. Prepare and merge data (2_prepare_data.R) 
# 3. Analyze PRS with incremental R² (3_analyze_prs.R)
#
# IMPORTANT: This pipeline calculates INCREMENTAL R² = R²_(PRS+PCs) - R²_(PCs only)
# This shows the true PRS contribution independent of population structure.
#
# Usage: ./run_prs_pipeline.sh [phenotype_file] [output_prefix]
#
# Author: Mantou (with Claude Code fixes for incremental R²)
# Date: 2025-12-22
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Default parameters
PHENOTYPE_FILE=""
OUTPUT_PREFIX="prs_analysis"
TRAIT_TYPE="binomial"  # Default to binomial (case-control)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Parse command line arguments
if [[ $# -eq 0 ]]; then
    cat << 'USAGE'
================================================================================
PRS ANALYSIS PIPELINE - COMPLETE WORKFLOW
================================================================================

This pipeline calculates polygenic risk scores and analyzes them using 
INCREMENTAL R² to show true PRS contribution separate from population structure.

USAGE:
  ./run_prs_pipeline.sh [phenotype_file] [output_prefix] [trait_type]

ARGUMENTS:
  phenotype_file    Path to phenotype/covariate file (optional)
  output_prefix     Prefix for output files (default: prs_analysis)
  trait_type        'binomial' for case-control or 'continuous' for quantitative (default: binomial)

PIPELINE STEPS:
  1. Calculate PRS from GWAS summary statistics
  2. Prepare and merge phenotype data  
  3. Analyze PRS with incremental R² calculation

CONFIGURATION FILES:
  - prs_config.sh:           PLINK parameters and file paths
  - data_prep_config.json:   Data merging and filtering settings

EXAMPLES:
  # Case-control analysis (default)
  ./run_prs_pipeline.sh /path/to/phenotypes.csv my_scz_analysis binomial
  
  # Continuous trait analysis  
  ./run_prs_pipeline.sh /path/to/phenotypes.csv my_height_analysis continuous

OUTPUT:
  - PRS scores: {output_prefix}.prs.*.sscore
  - Analysis results: {output_prefix}_results.csv
  - Summary report: {output_prefix}_summary.txt

IMPORTANT NOTES:
  - Results show INCREMENTAL R² (true PRS contribution)
  - PC-only R² represents population structure baseline
  - Use incremental R² for reporting PRS performance
  - Edit prs_config.sh and data_prep_config.json before running

================================================================================
USAGE
    exit 0
fi

# Parse arguments
if [[ $# -ge 1 ]]; then
    PHENOTYPE_FILE="$1"
fi

if [[ $# -ge 2 ]]; then
    OUTPUT_PREFIX="$2"
fi

if [[ $# -ge 3 ]]; then
    TRAIT_TYPE="$3"
    
    # Validate trait type
    if [[ "$TRAIT_TYPE" != "binomial" && "$TRAIT_TYPE" != "continuous" ]]; then
        echo "ERROR: trait_type must be 'binomial' or 'continuous', got: $TRAIT_TYPE"
        exit 1
    fi
fi

# Setup logging
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${OUTPUT_PREFIX}_pipeline_${TIMESTAMP}.log"

echo "================================================================================"
echo "PRS ANALYSIS PIPELINE WITH INCREMENTAL R²"
echo "================================================================================"
echo "Started: $(date)"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Trait type: $TRAIT_TYPE"
echo "Log file: $LOG_FILE"
echo "Script directory: $SCRIPT_DIR"
echo ""

# Log all output
exec 1> >(tee -a "$LOG_FILE")
exec 2>&1

# Function to check if previous step completed successfully
check_step() {
    local step_name="$1"
    local expected_files="$2"
    
    echo "[$(date +%H:%M:%S)] Checking completion of: $step_name"
    
    for file in $expected_files; do
        if [[ ! -f "$file" ]]; then
            echo "ERROR: Missing expected output: $file"
            echo "Step '$step_name' may have failed."
            exit 1
        fi
    done
    
    echo "✓ $step_name completed successfully"
    echo ""
}

# Change to script directory
cd "$SCRIPT_DIR"

# ============================================================================
# STEP 1: CALCULATE PRS SCORES
# ============================================================================

echo "================================================================================"
echo "STEP 1: CALCULATING PRS SCORES"
echo "================================================================================"
echo ""

if [[ ! -f "prs_config.sh" ]]; then
    echo "ERROR: prs_config.sh not found. Please create configuration file."
    echo "See prs_config.sh template in the repository."
    exit 1
fi

echo "Running PRS calculation with PLINK..."
echo "Command: ./1_calculate_prs.sh"
echo ""

./1_calculate_prs.sh

# Source config to get output location for verification
source prs_config.sh

# Check that PRS calculation completed
echo "Checking PRS calculation outputs..."
if ! ls "${OUTPUT_DIR}/${PROJECT_NAME}.prs."*.sscore 1> /dev/null 2>&1; then
    echo "ERROR: No .sscore files found in ${OUTPUT_DIR}"
    echo "PRS calculation may have failed. Check logs above."
    exit 1
fi

PRS_FILES=$(ls "${OUTPUT_DIR}/${PROJECT_NAME}.prs."*.sscore)
echo "✓ Found PRS score files:"
for file in $PRS_FILES; do
    echo "  - $(basename $file)"
done
echo ""

# ============================================================================
# STEP 2: PREPARE DATA FOR ANALYSIS
# ============================================================================

echo "================================================================================"
echo "STEP 2: PREPARING DATA FOR ANALYSIS"
echo "================================================================================"
echo ""

if [[ ! -f "data_prep_config.json" ]]; then
    echo "ERROR: data_prep_config.json not found. Please create configuration file."
    echo "See data_prep_config_examples.json for templates."
    exit 1
fi

# Determine phenotype file
if [[ -z "$PHENOTYPE_FILE" ]]; then
    echo "No phenotype file specified. Data preparation will use PRS data only."
    echo "WARNING: You'll need phenotype data for association analysis."
    PREPARED_DATA="${OUTPUT_PREFIX}_prepared_data.csv"
    echo "Command: Rscript 2_prepare_data.R data_prep_config.json $PREPARED_DATA"
    echo ""
    
    Rscript 2_prepare_data.R data_prep_config.json "$PREPARED_DATA"
else
    echo "Using phenotype file: $PHENOTYPE_FILE"
    PREPARED_DATA="${OUTPUT_PREFIX}_prepared_data.csv"
    echo "Command: Rscript 2_prepare_data.R data_prep_config.json $PREPARED_DATA $PHENOTYPE_FILE"
    echo ""
    
    Rscript 2_prepare_data.R data_prep_config.json "$PREPARED_DATA" "$PHENOTYPE_FILE"
fi

# Check data preparation completed
check_step "Data Preparation" "$PREPARED_DATA"

# Show data summary
echo "Data preparation summary:"
if [[ -f "$PREPARED_DATA" ]]; then
    if command -v csvstat >/dev/null 2>&1; then
        echo "Dataset dimensions:"
        csvstat --count "$PREPARED_DATA"
    else
        echo "Prepared data file: $PREPARED_DATA"
        echo "Size: $(wc -l < "$PREPARED_DATA") rows"
    fi
fi
echo ""

# ============================================================================
# STEP 3: ANALYZE PRS WITH INCREMENTAL R²
# ============================================================================

echo "================================================================================"
echo "STEP 3: PRS ANALYSIS WITH INCREMENTAL R²"
echo "================================================================================"
echo ""

ANALYSIS_OUTPUT_DIR="${OUTPUT_PREFIX}_results"
mkdir -p "$ANALYSIS_OUTPUT_DIR"

echo "Running PRS analysis with incremental R² calculation..."
echo "This will show true PRS contribution separate from population structure."
echo ""
echo "Command: Rscript 3_analyze_prs.R $PREPARED_DATA $ANALYSIS_OUTPUT_DIR $OUTPUT_PREFIX $TRAIT_TYPE"
echo ""

Rscript 3_analyze_prs.R "$PREPARED_DATA" "$ANALYSIS_OUTPUT_DIR" "$OUTPUT_PREFIX" "$TRAIT_TYPE"

# Check analysis completed
EXPECTED_RESULTS=(
    "${ANALYSIS_OUTPUT_DIR}/${OUTPUT_PREFIX}_results.csv"
    "${ANALYSIS_OUTPUT_DIR}/${OUTPUT_PREFIX}_summary.txt"
)

for file in "${EXPECTED_RESULTS[@]}"; do
    if [[ -f "$file" ]]; then
        echo "✓ Created: $file"
    else
        echo "⚠ Missing: $file"
    fi
done
echo ""

# ============================================================================
# PIPELINE COMPLETION SUMMARY
# ============================================================================

echo "================================================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY!"
echo "================================================================================"
echo ""
echo "Completed: $(date)"
echo "Total runtime: $SECONDS seconds"
echo ""

echo "OUTPUT FILES:"
echo "============="
echo ""

echo "1. PRS Scores (from PLINK):"
for file in $PRS_FILES; do
    echo "   $(basename $file)"
done
echo ""

echo "2. Prepared Data:"
echo "   $(basename $PREPARED_DATA)"
echo ""

echo "3. Analysis Results:"
for file in "${EXPECTED_RESULTS[@]}"; do
    if [[ -f "$file" ]]; then
        echo "   $(basename $file)"
    fi
done
echo ""

echo "KEY RESULTS INTERPRETATION:"
echo "=========================="
echo ""
echo "• PC-only R²: Variance explained by population structure alone"
echo "• Full R²: Total variance explained by PRS + population structure"  
echo "• Incremental PRS R²: TRUE PRS contribution (Full R² - PC-only R²)"
echo ""
echo "⭐ IMPORTANT: Use INCREMENTAL R² for reporting PRS performance"
echo "   This metric remains stable across different numbers of PCs"
echo "   and shows the true genetic signal separate from population structure."
echo ""

echo "NEXT STEPS:"
echo "=========="
echo "1. Review summary report: ${ANALYSIS_OUTPUT_DIR}/${OUTPUT_PREFIX}_summary.txt"
echo "2. Check detailed results: ${ANALYSIS_OUTPUT_DIR}/${OUTPUT_PREFIX}_results.csv"
echo "3. Use incremental R² values for publication/reporting"
echo "4. Validate results in independent dataset if available"
echo ""

echo "Pipeline log saved to: $LOG_FILE"
echo "================================================================================"