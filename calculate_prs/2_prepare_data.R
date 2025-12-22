#!/usr/bin/env Rscript

# Data preparation script for PRS analysis
# Handles all data formatting, merging, and filtering
#
# PURPOSE: This is a configurable data preparation script that applies
# user-defined filtering steps and produces a standardized output format.
#
# EXPECTED OUTPUT FORMAT (standardized for 3_analyze_prs.R):
# =================================================================
# REQUIRED COLUMNS (exact names):
#   - IID:          Sample identifier (character/string)
#   - case_control: Binary outcome (0 = control, 1 = case) 
#   - prs_score:    PRS value (numeric, any scale/normalization)
#   - PC1:          First principal component (numeric)
#   - PC2:          Second principal component (numeric)  
#   - ...
#   - PC20:         Tenth principal component (numeric)
#
# NOTE: PC columns can have any prefix (PC1, PC_EUR1, PC_AFR1, etc.)
#       The script will auto-detect and standardize to PC1-PC20
#
# METADATA COLUMNS (consistent across all rows):
#   - n_samples:    Total number of samples (integer)
#   - n_cases:      Number of cases (integer) 
#   - n_controls:   Number of controls (integer)
#
# FORMAT REQUIREMENTS:
#   - CSV file with header
#   - No missing values in required columns
#   - IID must be unique
#   - case_control must be 0 or 1 (no NAs)
#   - PC1-PC10 must be numeric (no NAs)
#   - prs_score can be any numeric scale (raw/normalized/standardized)
#
# NOTES FOR ADAPTATION:
#   - Modify the JSON config file to customize filtering steps
#   - Add/remove filter steps as needed for your dataset
#   - The final standardization step should remain unchanged
#   - The output format is a CONTRACT - 3_analyze_prs.R depends on it
#
# CONFIGURATION:
#   - Use data_prep_config.json to define:
#     * Merge column preferences
#     * Population filtering settings  
#     * Case/control definition criteria
#     * Principal component requirements
#     * PRS score column preferences
# =================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(jsonlite)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript 2_prepare_data.R <config_file> <pheno_file> <output_file>\n")
  cat("  config_file: Path to JSON config file (specifies PRS file locations)\n")
  cat("  pheno_file:  Path to phenotype CSV file\n") 
  cat("  output_file: Output path for prepared data\n")
  cat("\nExample: Rscript 2_prepare_data.R data_prep_config.json pheno.csv out.csv\n")
  cat("\nNote: PRS files are specified in the config file under 'prs_files' section\n")
  quit(status = 1)
}

config_file <- args[1]
pheno_file <- args[2]
output_file <- args[3]

# Load configuration
cat("Loading configuration from:", config_file, "\n")
if (!file.exists(config_file)) {
  stop("ERROR: Config file not found: ", config_file)
}
config <- fromJSON(config_file, simplifyVector = FALSE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
merge_cols <- unlist(config$merge_settings$columns_to_try)

cat("=================================================================\n")
cat("DATA PREPARATION FOR PRS ANALYSIS\n")
cat("=================================================================\n\n")

cat("Input files:\n")
cat("  Config file:", config_file, "\n")
cat("  Phenotype file:", pheno_file, "\n")
cat("  Output file:", output_file, "\n")
cat("  Merge columns to try:", paste(merge_cols, collapse = ", "), "\n")

# Show PRS file configuration
prs_config <- config$prs_files
cat("  PRS files config:\n")
cat("    - Type:", prs_config$input_type, "\n")
cat("    - Path:", prs_config$path, "\n")
if (prs_config$input_type == "directory") {
  cat("    - Pattern:", prs_config$file_pattern, "\n")
}

# Show enabled filter steps
enabled_steps <- Filter(function(x) x$enabled == TRUE, config$filter_steps)
enabled_names <- sapply(enabled_steps, function(x) x$name)
cat("  Enabled filter steps:", paste(enabled_names, collapse = ", "), "\n\n")

# Load and merge multiple PRS data files
cat("Loading PRS data...\n")
prs_config <- config$prs_files

if (prs_config$input_type == "directory") {
  # Find all PRS files in directory
  prs_dir <- prs_config$path
  if (!dir.exists(prs_dir)) {
    stop("ERROR: PRS directory not found: ", prs_dir)
  }
  
  # Get list of files matching pattern
  pattern <- gsub("\\*", ".*", prs_config$file_pattern)  # Convert glob to regex
  all_files <- list.files(prs_dir, full.names = TRUE)
  prs_files <- all_files[grepl(pattern, basename(all_files))]
  
  if (length(prs_files) == 0) {
    stop("ERROR: No PRS files found in ", prs_dir, " matching pattern ", prs_config$file_pattern)
  }
  
  cat("  - Found", length(prs_files), "PRS files\n")
  
  # Load and merge all PRS files
  prs_data_list <- list()
  score_col <- prs_config$score_column
  id_col <- prs_config$id_column
  
  for (i in seq_along(prs_files)) {
    file_path <- prs_files[i]
    file_name <- tools::file_path_sans_ext(basename(file_path))
    cat("    Loading:", file_name, "\n")
    
    # Load file
    temp_data <- fread(file_path, header = TRUE)
    
    # Check required columns exist
    if (!id_col %in% names(temp_data)) {
      stop("ERROR: ID column '", id_col, "' not found in ", file_path)
    }
    if (!score_col %in% names(temp_data)) {
      stop("ERROR: Score column '", score_col, "' not found in ", file_path)
    }
    
    # Extract ID and score columns, rename score column
    score_name <- paste0("prs_", gsub("[^a-zA-Z0-9_]", "_", file_name))
    temp_clean <- temp_data %>%
      select(!!id_col, !!score_col) %>%
      rename(
        IID = !!id_col,
        !!score_name := !!score_col
      )
    
    cat("      - Samples:", nrow(temp_clean), "\n")
    cat("      - Score column:", score_name, "\n")
    
    prs_data_list[[i]] <- temp_clean
  }
  
  # Merge all PRS data by IID
  cat("  - Merging PRS files...\n")
  prs_data <- prs_data_list[[1]]
  if (length(prs_data_list) > 1) {
    for (i in 2:length(prs_data_list)) {
      prs_data <- full_join(prs_data, prs_data_list[[i]], by = "IID")
    }
  }
  
  cat("  - Final PRS dataset:", nrow(prs_data), "samples x", ncol(prs_data)-1, "PRS scores\n")
  prs_score_cols <- names(prs_data)[names(prs_data) != "IID"]
  cat("  - PRS columns:", paste(prs_score_cols, collapse = ", "), "\n")
  
} else if (prs_config$input_type == "file") {
  # Single file mode
  prs_file <- prs_config$path
  if (!file.exists(prs_file)) {
    stop("ERROR: PRS file not found: ", prs_file)
  }
  
  prs_data <- fread(prs_file, header = TRUE)
  cat("  - Loaded", nrow(prs_data), "samples\n")
  cat("  - Columns:", paste(names(prs_data), collapse = ", "), "\n")
  prs_score_cols <- c(prs_config$score_column)
  
} else {
  stop("ERROR: Unsupported input_type: ", prs_config$input_type, ". Use 'directory' or 'file'.")
}

# Load phenotype data
cat("\nLoading phenotype data...\n")
if (!file.exists(pheno_file)) {
  stop("ERROR: Phenotype file not found: ", pheno_file)
}

pheno_data <- fread(pheno_file, header = TRUE)
cat("  - Loaded", nrow(pheno_data), "samples\n")

# Step 1: Merge datasets by flexible column matching
cat("\nStep 1: Merging datasets...\n")

# Try each merge column in order until we find a successful merge
merged_data <- NULL
successful_merge_col <- NULL

for (merge_col in merge_cols) {
  cat("  Trying merge column:", merge_col, "\n")
  
  # Check if column exists in both datasets
  pheno_has_col <- merge_col %in% names(pheno_data)
  prs_has_col <- merge_col %in% names(prs_data)
  
  cat("    - In phenotype file:", ifelse(pheno_has_col, "✓", "✗"), "\n")
  cat("    - In PRS file:", ifelse(prs_has_col, "✓", "✗"), "\n")
  
  if (pheno_has_col && prs_has_col) {
    # Attempt merge
    temp_merged <- inner_join(pheno_data, prs_data, by = merge_col)
    cat("    - Merge result:", nrow(temp_merged), "samples\n")
    
    if (nrow(temp_merged) > 0) {
      merged_data <- temp_merged
      successful_merge_col <- merge_col
      cat("  ✓ SUCCESS! Using merge column:", merge_col, "\n")
      break
    } else {
      cat("    ✗ No overlapping samples\n")
    }
  } else {
    cat("    ✗ Column missing from one or both files\n")
  }
}

if (is.null(merged_data) || nrow(merged_data) == 0) {
  cat("\nAvailable columns:\n")
  cat("  Phenotype file:", paste(names(pheno_data), collapse = ", "), "\n")
  cat("  PRS file:", paste(names(prs_data), collapse = ", "), "\n")
  stop("ERROR: No successful merge found with any of the specified columns: ", paste(merge_cols, collapse = ", "))
}

# Ensure we have an IID column for downstream analysis
if (!"IID" %in% names(merged_data)) {
  if (successful_merge_col != "IID") {
    cat("  - Renaming", successful_merge_col, "to IID for standardization\n")
    merged_data <- merged_data %>%
      rename(IID = !!successful_merge_col)
  }
}

cat("  - Final merged dataset:", nrow(merged_data), "samples\n")

# Step 2: Apply configurable filter steps
analysis_data <- merged_data
step_num <- 2

for (filter_step in config$filter_steps) {
  if (!filter_step$enabled) {
    next
  }
  
  cat("\nStep", step_num, ":", filter_step$name, "...\n")
  
  if (filter_step$name == "population_filter") {
    pop_settings <- filter_step$settings
    col_name <- pop_settings$column_name
    allowed_pops <- pop_settings$allowed_populations
    
    if (col_name %in% names(analysis_data)) {
      cat("  - Column:", col_name, "\n")
      cat("  - Allowed populations:", paste(allowed_pops, collapse = ", "), "\n")
      
      before_count <- nrow(analysis_data)
      analysis_data <- analysis_data %>%
        filter(!!sym(col_name) %in% allowed_pops)
      after_count <- nrow(analysis_data)
      
      cat("  - Before:", before_count, "samples\n")
      cat("  - After:", after_count, "samples\n")
      
      if (after_count == 0 && pop_settings$required) {
        stop("ERROR: No samples found for populations: ", paste(allowed_pops, collapse = ", "))
      }
    } else if (pop_settings$required) {
      stop("ERROR: Required population column not found: ", col_name)
    } else {
      cat("  WARNING: Population column '", col_name, "' not found. Skipping filter.\n")
    }
  }
  
  step_num <- step_num + 1
}

# Continue with remaining filter steps
for (filter_step in config$filter_steps) {
  if (!filter_step$enabled || filter_step$name == "population_filter") {
    next
  }
  
  cat("\nStep", step_num, ":", filter_step$name, "...\n")
  
  if (filter_step$name == "case_control_definition") {
    cc_settings <- filter_step$settings
    control_criteria <- cc_settings$control_criteria
    case_criteria <- cc_settings$case_criteria
    
    # Check if this is simple format (no conditions array) or complex format
    is_simple_format <- is.null(case_criteria$conditions)
    
    if (is_simple_format) {
      # Simple format: direct column = value mapping
      all_cols <- c(names(control_criteria), names(case_criteria))
      missing_cols <- all_cols[!all_cols %in% names(analysis_data)]
      
      if (length(missing_cols) > 0) {
        stop("ERROR: Missing required columns for case/control definition: ", 
             paste(missing_cols, collapse = ", "))
      }
      
      # Build conditions for simple format
      control_conditions <- paste(
        names(control_criteria), 
        "==", 
        unlist(control_criteria), 
        collapse = " & "
      )
      case_conditions <- paste(
        names(case_criteria), 
        "==", 
        unlist(case_criteria), 
        collapse = " & "
      )
      
      cat("  - Controls:", control_conditions, "\n")
      cat("  - Cases:", case_conditions, "\n")
      
      # Apply simple case/control definition
      analysis_data <- analysis_data %>%
        rowwise() %>%
        mutate(
          case_control = case_when(
            all(c_across(all_of(names(control_criteria))) == unlist(control_criteria)) ~ 0,
            all(c_across(all_of(names(case_criteria))) == unlist(case_criteria)) ~ 1,
            TRUE ~ NA_real_
          )
        ) %>%
        ungroup() %>%
        filter(!is.na(case_control))
        
    } else {
      # Complex format with OR conditions (original code)
      all_cols <- c(names(control_criteria), 
                    unlist(lapply(case_criteria$conditions, names)))
      missing_cols <- all_cols[!all_cols %in% names(analysis_data)]
      
      if (length(missing_cols) > 0) {
        stop("ERROR: Missing required columns for case/control definition: ", 
             paste(missing_cols, collapse = ", "))
      }
      
      # Build control condition
      control_conditions <- paste(
        names(control_criteria), 
        "==", 
        unlist(control_criteria), 
        collapse = " & "
      )
      cat("  - Controls:", control_conditions, "\n")
      
      # Build case condition  
      case_conditions <- sapply(case_criteria$conditions, function(cond) {
        paste(names(cond), "==", unlist(cond), collapse = " & ")
      })
      case_condition_str <- paste(case_conditions, collapse = paste0(" ", case_criteria$logic, " "))
      cat("  - Cases:", case_condition_str, "\n")
      
      # Apply complex case/control definition
      analysis_data <- analysis_data %>%
        rowwise() %>%
        mutate(
          case_control = case_when(
            # Control condition
            all(c_across(all_of(names(control_criteria))) == unlist(control_criteria)) ~ 0,
            # Case condition(s)
            any(sapply(case_criteria$conditions, function(cond) {
              all(c_across(all_of(names(cond))) == unlist(cond))
            })) ~ 1,
            TRUE ~ NA_real_
          )
        ) %>%
        ungroup() %>%
        filter(!is.na(case_control))
    }
  }
  
  step_num <- step_num + 1
}

cat("\nFiltering summary:\n")
cat("  - Cases:", sum(analysis_data$case_control == 1), "\n")
cat("  - Controls:", sum(analysis_data$case_control == 0), "\n")
cat("  - Total for analysis:", nrow(analysis_data), "\n")

# Validate principal components using config settings
cat("\nStep", step_num, ": Validating principal components...\n")
pc_config <- config$principal_components
min_required <- pc_config$min_required
pc_pattern <- pc_config$pattern

# Look for PC columns with configured pattern
all_pc_cols <- grep(pc_pattern, names(analysis_data), value = TRUE)
cat("  - Pattern used:", pc_pattern, "\n")
cat("  - Found", length(all_pc_cols), "total PC columns\n")

if (length(all_pc_cols) > 0) {
  # Extract numbers and sort to get PC1, PC2, etc. in order
  pc_numbers <- as.numeric(gsub(".*?([0-9]+)$", "\\1", all_pc_cols))
  pc_df <- data.frame(
    col_name = all_pc_cols,
    pc_num = pc_numbers
  ) %>%
    arrange(pc_num)
  
  # Take first min_required PCs
  pc_cols <- pc_df$col_name[1:min(min_required, nrow(pc_df))]
  cat("  - Using first", length(pc_cols), ":", paste(pc_cols, collapse = ", "), "\n")
} else {
  pc_cols <- character(0)
}

if (length(pc_cols) < min_required) {
  warning("WARNING: Only found ", length(pc_cols), " PC columns, need ", min_required, ". Will use what's available.")
  if (length(pc_cols) == 0) {
    stop("ERROR: No PC columns found matching pattern: ", pc_pattern)
  }
}
step_num <- step_num + 1

# Identify and standardize PRS score columns
cat("\nStep", step_num, ": Processing PRS score columns...\n")

# Find all PRS score columns in the merged dataset
all_prs_cols <- grep("^prs_", names(analysis_data), value = TRUE)
if (length(all_prs_cols) == 0) {
  # Fallback to single PRS column for backward compatibility
  prs_config <- config$prs_files
  if (prs_config$score_column %in% names(analysis_data)) {
    all_prs_cols <- prs_config$score_column
    # Rename to standard format
    analysis_data <- analysis_data %>%
      rename(prs_score1 = !!prs_config$score_column)
    all_prs_cols <- "prs_score1"
  } else {
    stop("ERROR: No PRS score columns found")
  }
}

cat("  - Found", length(all_prs_cols), "PRS score columns\n")
cat("  - Columns:", paste(all_prs_cols, collapse = ", "), "\n")

# For multi-PRS analysis, keep all PRS columns with their unique names
# No need for a generic 'prs_score' column when we have multiple specific ones
if (length(all_prs_cols) == 1) {
  # Single PRS case: create prs_score for backward compatibility
  first_prs_col <- all_prs_cols[1]
  analysis_data <- analysis_data %>%
    mutate(prs_score = !!sym(first_prs_col))
  cat("  - Single PRS: Using '", first_prs_col, "' as 'prs_score' for analysis\n")
} else {
  # Multi-PRS case: use only the specific named columns
  cat("  - Multi-PRS mode: Using", length(all_prs_cols), "specifically named PRS columns\n")
  cat("  - PRS columns:", paste(all_prs_cols, collapse = ", "), "\n")
}
step_num <- step_num + 1

# Final step: Create standardized dataset (DO NOT MODIFY)
# =================================================================
# This section creates the standardized output format that 
# 3_analyze_prs.R expects. Modify steps 1-5 above for your data,
# but keep this section unchanged to maintain the contract.
# =================================================================
cat("\nStep", step_num, ": Creating standardized output...\n")

# Standardize PC column names to PC1, PC2, ..., PC10
cat("  - Standardizing PC column names...\n")
for (i in 1:length(pc_cols)) {
  old_name <- pc_cols[i]
  new_name <- paste0("PC", i)
  if (old_name != new_name) {
    analysis_data <- analysis_data %>%
      rename(!!new_name := !!old_name)
  }
}

# Select and order columns for analysis
base_cols <- c(
  "IID", 
  "case_control", 
  paste0("PC", 1:length(pc_cols))
)

# Add PRS columns based on analysis mode
all_prs_cols_in_data <- grep("^prs_", names(analysis_data), value = TRUE)

if (length(all_prs_cols) == 1) {
  # Single PRS: include both the original named column and prs_score
  base_cols <- c(base_cols, "prs_score")
  standard_cols <- c(base_cols, all_prs_cols_in_data[all_prs_cols_in_data != "prs_score"])
  cat("  - Single PRS mode: including prs_score column\n")
} else {
  # Multi-PRS: include only the specifically named columns
  standard_cols <- c(base_cols, all_prs_cols_in_data)
  cat("  - Multi-PRS mode: including", length(all_prs_cols_in_data), "named PRS columns\n")
}

cat("  - PRS columns in output:", paste(all_prs_cols_in_data, collapse = ", "), "\n")

# Check essential columns exist
essential_cols <- c("IID", "case_control", paste0("PC", 1:length(pc_cols)))
missing_essential <- essential_cols[!essential_cols %in% names(analysis_data)]
if (length(missing_essential) > 0) {
  stop("ERROR: Missing essential columns: ", paste(missing_essential, collapse = ", "))
}

final_data <- analysis_data %>%
  select(all_of(standard_cols)) %>%
  mutate(
    n_samples = nrow(.),
    n_cases = sum(case_control == 1),
    n_controls = sum(case_control == 0)
  )

cat("  - Final dataset shape: ", nrow(final_data), "samples x", ncol(final_data), "columns\n")

# Save prepared data
cat("\nSaving prepared data...\n")
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("  Created directory:", output_dir, "\n")
}

fwrite(final_data, output_file)
cat("  ✓ Saved to:", output_file, "\n")

# Summary
cat("\n=================================================================\n")
cat("DATA PREPARATION COMPLETE\n")
cat("=================================================================\n")
cat("Final dataset:\n")
cat("  - Samples:", nrow(final_data), "\n")
cat("  - Cases:", final_data$n_cases[1], "\n") 
cat("  - Controls:", final_data$n_controls[1], "\n")
cat("  - PRS column: prs_score (standardized)\n")
cat("  - Principal components: PC1-", length(pc_cols), "\n", sep="")
cat("  - Ready for analysis!\n")
