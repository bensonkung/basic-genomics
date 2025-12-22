#!/bin/bash

################################################################################
# PRS CALCULATION CONFIGURATION FILE
# 
# This file contains all parameters for polygenic risk score calculation
# Modify these variables to adapt the pipeline to your data
#
# Author: Mantou
# Date: 2025-12-17
################################################################################

# ============================================================================
# INPUT DATA PATHS
# ============================================================================

# Target genotype data (the samples you want to score)
# PLINK binary format: .bed/.bim/.fam or PLINK2 format: .pgen/.pvar/.psam
TARGET_GENO="/home/groups/laramied/data/GWAS_data/PRIVATE/HBCC/COMBINED_FILES_HBCC/HBCC_imputed_clean_4datasets_HY_01_15_2021"  # Without extension
TARGET_FORMAT="bed"  # Options: "bed" (PLINK1.9) or "pgen" (PLINK2)

# ============================================================================
# SAMPLE SUBSET / POPULATION SELECTION
# ============================================================================

# Keep only specific samples (e.g., specific population/ancestry)
# Leave empty to use all samples
KEEP_SAMPLES=""  # Path to file with FID IID to keep (one per line)

# Remove specific samples
# Leave empty to keep all samples
REMOVE_SAMPLES=""  # Path to file with FID IID to remove (one per line)

# GWAS summary statistics file
# Should contain columns: SNP ID, effect allele, effect size (BETA or OR), p-value
# NOTE: Comment lines starting with '#' will be automatically skipped
GWAS_SUMSTATS="/home/groups/laramied/data/GWAS_results/PUBLIC/SCHIZOPHRENIA/PGC_SCZ_2022/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.no_heading"
GWAS_HAS_HEADER=true  # Set to false if no header (after comment lines)

# Column numbers in GWAS file (1-indexed)
# Check your file and adjust accordingly
SNP_COL=2        # SNP ID column
A1_COL=4         # Effect allele column  
BETA_COL=9      # Effect size column (BETA or log(OR))
PVAL_COL=11       # P-value column

# Does the GWAS file contain OR instead of BETA?
# If true, script will log-transform OR to BETA
USE_OR=false

# ============================================================================
# OUTPUT DIRECTORY
# ============================================================================

OUTPUT_DIR="/scratch/users/benson97/SCZ_PRS"
PROJECT_NAME="PGC_SCZ_2022_standard"

# ============================================================================
# QUALITY CONTROL THRESHOLDS - BASE DATA (GWAS)
# ============================================================================

# Minor allele frequency threshold
# Variants below this MAF will be excluded
BASE_MAF=0.01

# INFO score threshold (for imputed data)
# Variants below this will be excluded (set to 0 if not using imputed data)
BASE_INFO=0.8

# Remove ambiguous SNPs (A/T, G/C)?
# These can cause strand alignment issues
REMOVE_AMBIGUOUS=true

# Remove duplicated SNPs?
REMOVE_DUPLICATES=true

# ============================================================================
# QUALITY CONTROL THRESHOLDS - TARGET DATA
# ============================================================================

# Genotype missingness per SNP (e.g., 0.05 = remove SNPs with >5% missing)
TARGET_GENO_MISSING=0.05

# Individual missingness (e.g., 0.05 = remove samples with >5% missing)
TARGET_IND_MISSING=0.05

# Minor allele frequency in target
TARGET_MAF=0.01

# Hardy-Weinberg equilibrium p-value threshold
TARGET_HWE=1e-6

# Keep variants that fail HWE due to too few heterozygotes?
# (Recommended when population stratification is present)
HWE_KEEP_FEWHET=true

# ============================================================================
# LINKAGE DISEQUILIBRIUM (LD) CLUMPING PARAMETERS
# ============================================================================

# Clumping is used to remove SNPs in LD, keeping only independent signals
# This reduces redundancy in the PRS

# P-value threshold for index SNPs
CLUMP_P1=1.0

# P-value threshold for clumped SNPs
CLUMP_P2=1.0

# LD threshold (r²) - SNPs with r² above this with index SNP are removed
CLUMP_R2=0.1

# Physical distance threshold (kb) for clumping window
CLUMP_KB=250

# Reference panel for LD calculation (leave empty to use target data)
# If using different ancestry, provide appropriate reference
LD_REFERENCE="/oak/stanford/groups/laramied/1000Genomes/ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/EUR_1000_Genomes_biallelic_SNPs_maf001"

# ============================================================================
# PRS P-VALUE THRESHOLDS
# ============================================================================

# Define p-value thresholds for including SNPs in PRS
# Format: "threshold_name min_p max_p"
# Multiple thresholds allow finding the "best-fit" PRS

PRS_THRESHOLDS=(
    "5e-8 0 5e-8"      # Genome-wide significant
    "1e-6 0 1e-6"      
    "1e-4 0 1e-4"      
    "0.001 0 0.001"    
    "0.01 0 0.01"      
    "0.05 0 0.05"      
    "0.1 0 0.1"        
    "0.5 0 0.5"        
    "1.0 0 1.0"        # All SNPs
)

# ============================================================================
# SCORE CALCULATION OPTIONS
# ============================================================================

# Scoring method
# "average" - divide by number of non-missing alleles (PLINK default)
# "sum" - simple sum (multiply effect size by allele count)
SCORE_METHOD="average"  # Options: "average" or "sum"

# How to handle missing genotypes?
# "mean_impute" - use population MAF
# "no_mean_impute" - exclude from calculation
MISSING_GENOTYPE="mean_impute"

# ============================================================================
# COMPUTATIONAL RESOURCES
# ============================================================================

# Number of threads for parallel processing
THREADS=$(nproc)  # Auto-detect available cores, or set manually

# Memory limit (MB) - adjust based on your system
MEMORY=8000

# ============================================================================
# ADDITIONAL FILTERS (OPTIONAL)
# ============================================================================

# Exclude specific genomic regions (e.g., MHC region, high-LD regions)
# Format: chr:start-end
EXCLUDE_REGIONS=(
)

# Exclude specific SNPs (provide file path or leave empty)
# File should contain one SNP ID per line
EXCLUDE_SNPS=""

# Include only specific SNPs (provide file path or leave empty)
INCLUDE_SNPS=""

# ============================================================================
# ADVANCED OPTIONS
# ============================================================================

# Perform sex check?
SEX_CHECK=false

# Check for relatedness? (remove one of each pair with KING kinship > threshold)
RELATEDNESS_CHECK=true
RELATEDNESS_THRESHOLD=0.177   # 1st degree relatives

# Calculate principal components for population stratification?
CALCULATE_PCS=false
NUM_PCS=10

# ============================================================================
# VALIDATION OPTIONS
# ============================================================================

# Run comprehensive validation checks before main analysis?
RUN_VALIDATION=true

# Generate detailed QC report?
GENERATE_QC_REPORT=true

# ============================================================================
# FILE RETENTION CONTROL
# ============================================================================

# Keep intermediate genotype files (.bed/.bim/.fam)?
# These are: .target.geno.*, .target.mind.*, .target.maf.*
KEEP_INTERMEDIATE_GENOTYPES=false

# Keep QC'd data files?
# These are: ${PROJECT_NAME}.QC (base data), ${PROJECT_NAME}.target.qc.* (target genotypes)
KEEP_QC_FILES=false

# Keep clumping files?
# These are: ${PROJECT_NAME}.clumped.*, ${PROJECT_NAME}.clumped.valid.snp
KEEP_CLUMPING_FILES=false

# Keep temporary processing files?
# These are: qc_base.R, score_file.txt, SNP.pvalue, range_list, sex_discordance.txt
KEEP_TEMP_FILES=false

# Keep all log files (.log)?
# Recommended to keep for troubleshooting
KEEP_LOGS=false

# Summary: If ALL retention flags are false and GENERATE_QC_REPORT=false,
# only .sscore files will remain (the actual PRS results)

################################################################################
# END OF CONFIGURATION
################################################################################
