#!/usr/bin/env Rscript

# Clean PRS statistical analysis script with INCREMENTAL R²
# Expects standardized input format from 2_prepare_data.R
# 
# IMPORTANT: This script calculates incremental R² = R²_(PRS+PCs) - R²_(PCs only)
# This shows the true PRS contribution independent of population structure

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(broom)
})

# Check if fmsb is installed for Nagelkerke's R2
if (!requireNamespace("fmsb", quietly = TRUE)) {
  cat("Note: fmsb not installed. Will calculate Nagelkerke's R2 manually.\n")
  use_fmsb <- FALSE
} else {
  library(fmsb)
  use_fmsb <- TRUE
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript 3_analyze_prs.R <input_file> <output_dir> [analysis_name] [trait_type]\n")
  cat("  input_file:     Prepared data from prepare_data.R\n")
  cat("  output_dir:     Directory for results\n")
  cat("  analysis_name:  Name prefix for output files (default: prs_analysis)\n")
  cat("  trait_type:     'binomial' for case-control or 'continuous' for quantitative (default: binomial)\n")
  cat("\n")
  cat("Examples:\n")
  cat("  Rscript 3_analyze_prs.R data.csv results/ scz_analysis binomial\n")
  cat("  Rscript 3_analyze_prs.R data.csv results/ height_analysis continuous\n")
  quit(status = 1)
}

input_file <- args[1]
output_dir <- args[2]
analysis_name <- ifelse(length(args) >= 3, args[3], "prs_analysis")
trait_type <- ifelse(length(args) >= 4, args[4], "binomial")

# Validate trait type
if (!trait_type %in% c("binomial", "continuous")) {
  stop("ERROR: trait_type must be 'binomial' or 'continuous', got: ", trait_type)
}

cat("=================================================================\n")
cat("PRS STATISTICAL ANALYSIS\n")
cat("=================================================================\n\n")

cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Analysis name:", analysis_name, "\n")
cat("Trait type:", trait_type, "\n\n")

# Function to validate trait type against data
validate_trait_type <- function(data, trait_type) {
  if (!"case_control" %in% names(data)) {
    stop("ERROR: 'case_control' column not found in data")
  }
  
  outcome_values <- unique(data$case_control)
  outcome_values <- outcome_values[!is.na(outcome_values)]
  
  if (trait_type == "binomial") {
    # For binomial, expect only 0 and 1
    if (!all(outcome_values %in% c(0, 1))) {
      cat("WARNING: Binomial trait specified but case_control contains non-binary values\n")
      cat("Found values:", paste(sort(outcome_values), collapse = ", "), "\n")
      cat("Expected: 0, 1\n")
      cat("Consider using trait_type='continuous' for quantitative traits\n\n")
      
      # Allow user to continue but warn
      if (length(outcome_values) > 2) {
        stop("ERROR: Too many unique values for binomial trait. Use 'continuous' trait type.")
      }
    }
    
    n_cases <- sum(data$case_control == 1, na.rm = TRUE)
    n_controls <- sum(data$case_control == 0, na.rm = TRUE)
    cat("Binomial trait validation:\n")
    cat("  Cases (1):", n_cases, "\n")
    cat("  Controls (0):", n_controls, "\n")
    
    if (n_cases < 10 || n_controls < 10) {
      cat("  WARNING: Very small group sizes may lead to unstable results\n")
    }
    
  } else if (trait_type == "continuous") {
    # For continuous, expect more than 2 unique values
    if (length(outcome_values) <= 2 && all(outcome_values %in% c(0, 1))) {
      cat("WARNING: Continuous trait specified but case_control appears binary (0, 1)\n")
      cat("Consider using trait_type='binomial' for case-control traits\n\n")
    }
    
    cat("Continuous trait validation:\n")
    cat("  Range:", round(min(data$case_control, na.rm = TRUE), 3), "to", 
        round(max(data$case_control, na.rm = TRUE), 3), "\n")
    cat("  Mean:", round(mean(data$case_control, na.rm = TRUE), 3), "\n")
    cat("  SD:", round(sd(data$case_control, na.rm = TRUE), 3), "\n")
    cat("  Unique values:", length(outcome_values), "\n")
  }
  
  cat("✓ Trait type validation complete\n\n")
}

# Function to calculate R² based on trait type
calculate_r_squared <- function(model, trait_type) {
  if (trait_type == "binomial") {
    # Use Nagelkerke R² for logistic regression
    if (use_fmsb) {
      return(NagelkerkeR2(model)$R2)
    } else {
      return(nagelkerke_r2_manual(model))
    }
  } else {
    # Use standard R² for linear regression
    return(summary(model)$r.squared)
  }
}

# Function to calculate Nagelkerke's R2 manually with validation
nagelkerke_r2_manual <- function(model) {
  null_ll <- logLik(update(model, . ~ 1))
  full_ll <- logLik(model)
  n <- length(model$fitted.values)
  
  # Calculate Cox-Snell R²
  cox_snell <- 1 - exp((2/n) * (null_ll - full_ll))
  
  # Calculate maximum possible R² for this sample size
  max_r2 <- 1 - exp((2/n) * null_ll)
  
  # Nagelkerke R² = Cox-Snell R² / Max R²
  nagelkerke <- cox_snell / max_r2
  
  # Validation checks
  if (nagelkerke < 0 || nagelkerke > 1) {
    warning("Nagelkerke R² out of bounds [0,1]: ", nagelkerke)
  }
  
  if (is.na(nagelkerke) || is.infinite(nagelkerke)) {
    warning("Invalid Nagelkerke R² calculation")
    return(0)
  }
  
  return(as.numeric(nagelkerke))
}

# Function to run PC-only regression (no PRS)
run_pc_only_regression <- function(data, pc_count, analysis_label = "", pc_col_prefix = "PC", all_pc_cols, trait_type = "binomial") {
  # Build PC variable names based on available columns
  if (pc_col_prefix == "PC_EUR") {
    pc_vars <- paste0("PC_EUR", 1:pc_count)
  } else {
    pc_vars <- paste0("PC", 1:pc_count)
  }
  
  # Check if PC columns exist, if not use available ones
  available_pc_vars <- pc_vars[pc_vars %in% names(data)]
  if (length(available_pc_vars) < pc_count) {
    # Try alternative PC naming
    alt_pc_vars <- all_pc_cols[1:min(pc_count, length(all_pc_cols))]
    if (length(alt_pc_vars) >= pc_count) {
      pc_vars <- alt_pc_vars
    } else {
      warning("Only ", length(available_pc_vars), " PC columns available, using all available")
      pc_vars <- available_pc_vars
    }
  }
  
  formula_str <- paste0("case_control ~ ", paste(pc_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  cat("  Running PC-only model with", length(pc_vars), "PCs...\n")
  cat("  Formula:", formula_str, "\n")
  
  tryCatch({
    # Use appropriate GLM family based on trait type
    if (trait_type == "binomial") {
      model <- glm(formula_obj, data = data, family = binomial())
    } else {
      model <- glm(formula_obj, data = data, family = gaussian())
    }
    
    # Calculate appropriate R²
    r2 <- calculate_r_squared(model, trait_type)
    
    results <- data.frame(
      prs_column = "PC_only",
      analysis = ifelse(analysis_label == "", paste0("PC_only_", pc_count, "PCs"), analysis_label),
      pc_count = pc_count,
      beta = NA,  # No PRS coefficient
      se = NA,
      z_value = NA,
      p_value = NA,  # No PRS p-value
      OR = NA,
      OR_lower_95CI = NA,
      OR_upper_95CI = NA,
      nagelkerke_r2 = r2,
      n_samples = nrow(data),
      n_cases = sum(data$case_control == 1),
      n_controls = sum(data$case_control == 0)
    )
    
    cat("  ✓ PC-only model successful\n")
    r2_label <- if (trait_type == "binomial") "Nagelkerke R²" else "R²"
    cat("   ", r2_label, ":", round(results$nagelkerke_r2, 4), "\n\n")
    
    # Store full model for potential additional analyses
    results$model <- list(model)
    
    return(results)
  }, error = function(e) {
    cat("  ✗ PC-only model failed:", as.character(e), "\n\n")
    return(data.frame(
      prs_column = "PC_only",
      analysis = ifelse(analysis_label == "", paste0("PC_only_", pc_count, "PCs"), analysis_label),
      pc_count = pc_count,
      beta = NA,
      se = NA,
      z_value = NA,
      p_value = NA,
      OR = NA,
      OR_lower_95CI = NA,
      OR_upper_95CI = NA,
      nagelkerke_r2 = NA,
      n_samples = nrow(data),
      n_cases = sum(data$case_control == 1),
      n_controls = sum(data$case_control == 0),
      error = as.character(e),
      model = list(NULL)
    ))
  })
}

# Function to run regression and extract results with INCREMENTAL R²
run_regression <- function(data, pc_count, prs_col = "prs_score", analysis_label = "", pc_col_prefix = "PC", all_pc_cols, trait_type = "binomial") {
  # Build PC variable names based on available columns
  if (pc_col_prefix == "PC_EUR") {
    pc_vars <- paste0("PC_EUR", 1:pc_count)
  } else {
    pc_vars <- paste0("PC", 1:pc_count)
  }
  
  # Check if PC columns exist, if not use available ones
  available_pc_vars <- pc_vars[pc_vars %in% names(data)]
  if (length(available_pc_vars) < pc_count) {
    # Try alternative PC naming
    alt_pc_vars <- all_pc_cols[1:min(pc_count, length(all_pc_cols))]
    if (length(alt_pc_vars) >= pc_count) {
      pc_vars <- alt_pc_vars
    } else {
      warning("Only ", length(available_pc_vars), " PC columns available, using all available")
      pc_vars <- available_pc_vars
    }
  }
  
  cat("  Running", analysis_label, "model with", length(pc_vars), "PCs for", prs_col, "...\n")
  
  tryCatch({
    # Model 1: PC-only baseline (for incremental R² calculation)
    pc_formula <- as.formula(paste0("case_control ~ ", paste(pc_vars, collapse = " + ")))
    if (trait_type == "binomial") {
      model_pc_only <- glm(pc_formula, data = data, family = binomial())
    } else {
      model_pc_only <- glm(pc_formula, data = data, family = gaussian())
    }
    
    r2_pc_only <- calculate_r_squared(model_pc_only, trait_type)
    
    # Model 2: Full model (PRS + PCs)
    full_formula <- as.formula(paste0("case_control ~ ", prs_col, " + ", paste(pc_vars, collapse = " + ")))
    if (trait_type == "binomial") {
      model_full <- glm(full_formula, data = data, family = binomial())
    } else {
      model_full <- glm(full_formula, data = data, family = gaussian())
    }
    
    r2_full <- calculate_r_squared(model_full, trait_type)
    
    # Calculate INCREMENTAL R² (this is the key fix!)
    r2_incremental_prs <- r2_full - r2_pc_only
    
    # Extract PRS coefficient (use the actual PRS column name)
    prs_coef <- coef(summary(model_full))[prs_col, ]
    
    # Handle different test statistic names and calculate OR for binomial traits only
    if (trait_type == "binomial") {
      # Logistic regression uses z-value and Pr(>|z|)
      test_stat <- prs_coef["z value"]
      p_value <- prs_coef["Pr(>|z|)"]
      or <- exp(prs_coef["Estimate"])
      or_lower <- exp(prs_coef["Estimate"] - 1.96 * prs_coef["Std. Error"])
      or_upper <- exp(prs_coef["Estimate"] + 1.96 * prs_coef["Std. Error"])
    } else {
      # Linear regression uses t-value and Pr(>|t|)
      test_stat <- prs_coef["t value"]
      p_value <- prs_coef["Pr(>|t|)"]
      # For continuous traits, OR doesn't make sense
      or <- NA
      or_lower <- NA
      or_upper <- NA
    }
    
    results <- data.frame(
      prs_column = prs_col,
      analysis = ifelse(analysis_label == "", paste0(pc_count, "_PCs"), analysis_label),
      pc_count = pc_count,
      beta = prs_coef["Estimate"],
      se = prs_coef["Std. Error"],
      z_value = test_stat,  # This will be z for binomial, t for continuous
      p_value = p_value,
      OR = or,
      OR_lower_95CI = or_lower,
      OR_upper_95CI = or_upper,
      r2_pc_only = r2_pc_only,
      r2_full = r2_full,
      r2_incremental_prs = r2_incremental_prs,  # This is the key metric!
      nagelkerke_r2 = r2_incremental_prs,  # Report incremental R² as main metric
      n_samples = nrow(data),
      n_cases = sum(data$case_control == 1),
      n_controls = sum(data$case_control == 0)
    )
    
    cat("  ✓ Model successful\n")
    cat("    Beta:", round(results$beta, 6), "\n")
    cat("    P-value:", format(results$p_value, scientific = TRUE), "\n")
    
    # Only show OR for binomial traits
    if (trait_type == "binomial") {
      cat("    OR:", round(results$OR, 3), "(95% CI:", round(results$OR_lower_95CI, 3), "-", round(results$OR_upper_95CI, 3), ")\n")
    }
    
    r2_label <- if (trait_type == "binomial") "Nagelkerke R²" else "R²"
    cat("    PC-only", r2_label, ":", round(r2_pc_only, 4), "(population structure)\n")
    cat("    Full model", r2_label, ":", round(r2_full, 4), "(total variance)\n")
    cat("    Incremental PRS", r2_label, ":", round(r2_incremental_prs, 4), "(PRS contribution)\n")
    
    # Diagnostic warnings for incremental R²
    if (r2_incremental_prs < 0) {
      cat("    WARNING: Negative incremental R² suggests model issues\n")
    } else if (r2_incremental_prs > 0.3) {
      cat("    NOTE: High PRS R² (>30%) - check for data issues\n")
    }
    cat("\n")
    
    # Store full model for potential additional analyses
    results$model <- list(model_full)
    
    return(results)
  }, error = function(e) {
    cat("  ✗ Model failed:", as.character(e), "\n\n")
    return(data.frame(
      prs_column = prs_col,
      analysis = ifelse(analysis_label == "", paste0(pc_count, "_PCs"), analysis_label),
      pc_count = pc_count,
      beta = NA,
      se = NA,
      z_value = NA,
      p_value = NA,
      OR = NA,
      OR_lower_95CI = NA,
      OR_upper_95CI = NA,
      r2_pc_only = NA,
      r2_full = NA,
      r2_incremental_prs = NA,
      nagelkerke_r2 = NA,
      n_samples = nrow(data),
      n_cases = sum(data$case_control == 1),
      n_controls = sum(data$case_control == 0),
      error = as.character(e),
      model = list(NULL)
    ))
  })
}

# Load prepared data
cat("Loading prepared data...\n")
if (!file.exists(input_file)) {
  stop("ERROR: Input file not found: ", input_file)
}

analysis_data <- fread(input_file, header = TRUE)
cat("  - Loaded", nrow(analysis_data), "samples\n")
cat("  - Cases:", sum(analysis_data$case_control == 1), "\n")
cat("  - Controls:", sum(analysis_data$case_control == 0), "\n")

# Check for potential data issues that could inflate R²
case_prop <- mean(analysis_data$case_control)
cat("  - Case proportion:", round(case_prop, 3), "\n")

if (case_prop < 0.01 || case_prop > 0.99) {
  cat("  WARNING: Extreme case/control imbalance detected!\n")
  cat("           This can lead to inflated R² values.\n")
}

# Find PC columns (flexible naming)
pc_patterns <- c("^PC\\d+$", "^PC_EUR\\d+$", "^PC_AFR\\d+$", "^PC1$")
all_pc_cols <- character(0)
for (pattern in pc_patterns) {
  pc_matches <- grep(pattern, names(analysis_data), value = TRUE)
  if (length(pc_matches) > 0) {
    all_pc_cols <- pc_matches
    break
  }
}

# Check for population stratification in PC space
if (length(all_pc_cols) >= 2) {
  pc1_cases <- analysis_data[analysis_data$case_control == 1, all_pc_cols[1], with=FALSE]
  pc1_controls <- analysis_data[analysis_data$case_control == 0, all_pc_cols[1], with=FALSE]
  pc1_t_test <- t.test(pc1_cases[[1]], pc1_controls[[1]])
  
  if (pc1_t_test$p.value < 1e-10) {
    cat("  WARNING: Strong PC1 difference between cases/controls (p < 1e-10)\n")
    cat("           This suggests population stratification.\n")
  }
}

# Find all PRS columns
prs_cols <- grep("^prs_", names(analysis_data), value = TRUE)

# If no prs_* columns found, fall back to prs_score for single PRS analysis
if (length(prs_cols) == 0 && "prs_score" %in% names(analysis_data)) {
  prs_cols <- "prs_score"
  cat("  - Single PRS mode: using 'prs_score' column\n")
} else if (length(prs_cols) > 0) {
  # Multi-PRS mode: use only the specifically named columns, skip generic prs_score
  prs_cols <- prs_cols[prs_cols != "prs_score"]
  cat("  - Multi-PRS mode: using", length(prs_cols), "specifically named PRS columns\n")
} else {
  stop("ERROR: No PRS columns found")
}

cat("  - Found", length(prs_cols), "PRS columns:", paste(prs_cols, collapse = ", "), "\n")

if (length(all_pc_cols) < 3) {
  stop("ERROR: Need at least 3 PC columns, found: ", length(all_pc_cols))
}

cat("  - Found", length(all_pc_cols), "PC columns:", paste(all_pc_cols[1:min(10, length(all_pc_cols))], collapse = ", "), "\n")

# Validate essential columns
essential_cols <- c("IID", "case_control")
missing_essential <- essential_cols[!essential_cols %in% names(analysis_data)]

if (length(missing_essential) > 0) {
  stop("ERROR: Missing essential columns: ", paste(missing_essential, collapse = ", "))
}

cat("  ✓ All essential columns present\n\n")

# Validate trait type against data
validate_trait_type(analysis_data, trait_type)

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
}

# Run analyses
cat("=================================================================\n")
if (trait_type == "binomial") {
  cat("RUNNING LOGISTIC REGRESSION MODELS FOR", length(prs_cols), "PRS COLUMNS\n")
} else {
  cat("RUNNING LINEAR REGRESSION MODELS FOR", length(prs_cols), "PRS COLUMNS\n")
}
cat("=================================================================\n\n")

# Determine PC prefix from available columns
pc_prefix <- "PC"
if (any(grepl("^PC_EUR", all_pc_cols))) {
  pc_prefix <- "PC_EUR"
}

# Run analyses for each PRS column
all_results_list <- list()
result_counter <- 1

# First, run PC-only baseline models
cat("Running PC-only baseline models...\n")
cat("=====================================\n")

max_pcs <- min(10, length(all_pc_cols))

for (pc_count in c(3, max_pcs)) {
  if (pc_count <= length(all_pc_cols)) {
    label <- paste0("PC_only_", pc_count, "PC")
    result <- run_pc_only_regression(analysis_data, pc_count, label, pc_prefix, all_pc_cols, trait_type)
    
    all_results_list[[result_counter]] <- result %>% select(-model)
    result_counter <- result_counter + 1
  }
}

cat("\n")

# Now run PRS analyses
for (prs_col in prs_cols) {
  cat("Analyzing PRS column:", prs_col, "\n")
  cat("-----------------------------------\n")
  
  # Check if the PRS column has any variation
  if (all(is.na(analysis_data[[prs_col]])) || var(analysis_data[[prs_col]], na.rm = TRUE) == 0) {
    cat("  WARNING: PRS column", prs_col, "has no variation, skipping\n\n")
    next
  }
  
  # Standard analyses: 3 PCs and 10 PCs for each PRS
  max_pcs <- min(20, length(all_pc_cols))
  
  for (pc_count in 3:max_pcs) {
    if (pc_count <= length(all_pc_cols)) {
      label <- paste0(prs_col, "_", pc_count, "PC")
      result <- run_regression(analysis_data, pc_count, prs_col, label, pc_prefix, all_pc_cols, trait_type)
      
      all_results_list[[result_counter]] <- result %>% select(-model)
      result_counter <- result_counter + 1
    }
  }
  cat("\n")
}

# Combine all results
if (length(all_results_list) > 0) {
  all_results <- bind_rows(all_results_list)
} else {
  stop("ERROR: No valid PRS analyses completed")
}

# Save detailed results
cat("=================================================================\n")
cat("SAVING RESULTS\n")
cat("=================================================================\n\n")

results_file <- file.path(output_dir, paste0(analysis_name, "_results.csv"))
fwrite(all_results, results_file)
cat("✓ Detailed results saved to:", results_file, "\n")

# Create summary report
summary_file <- file.path(output_dir, paste0(analysis_name, "_summary.txt"))
sink(summary_file)

cat("=================================================================\n")
cat("PRS ANALYSIS SUMMARY REPORT\n")
cat("=================================================================\n\n")

cat("Analysis:", analysis_name, "\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("Input file:", input_file, "\n\n")

cat("DATASET CHARACTERISTICS\n")
cat("=================================================================\n")
cat("Total samples:", nrow(analysis_data), "\n")
cat("Cases:", sum(analysis_data$case_control == 1), "\n")
cat("Controls:", sum(analysis_data$case_control == 0), "\n")
cat("Case proportion:", round(mean(analysis_data$case_control), 3), "\n")
cat("PRS columns analyzed:", length(prs_cols), "\n")
cat("PRS columns:", paste(prs_cols, collapse = ", "), "\n\n")

# Best performing PRS comparison (using 10 PC models for fair comparison, excluding PC-only)
max_pcs <- min(10, length(all_pc_cols))
best_results <- all_results %>%
  filter(pc_count == max_pcs, !is.na(p_value), prs_column != "PC_only") %>%
  arrange(p_value)

# Get PC-only baseline results for comparison
pc_only_results <- all_results %>%
  filter(prs_column == "PC_only") %>%
  arrange(pc_count)

# Show PC-only baseline first
if (nrow(pc_only_results) > 0) {
  cat("PC-ONLY BASELINE MODELS\n")
  cat("=================================================================\n")
  
  for (i in 1:nrow(pc_only_results)) {
    result <- pc_only_results[i, ]
    cat(sprintf("%d PCs: R² = %.2f%% (population structure baseline)\n", 
               result$pc_count, result$nagelkerke_r2 * 100))
  }
  cat("\n")
}

if (nrow(best_results) > 0) {
  cat("BEST PERFORMING PRS (", max_pcs, " PC models)\n")
  cat("=================================================================\n")
  
  # Show baseline for comparison
  baseline_r2 <- pc_only_results[pc_only_results$pc_count == max_pcs, "nagelkerke_r2"]
  if (length(baseline_r2) > 0) {
    cat("Baseline (PC-only):", sprintf("R² = %.2f%%\n", baseline_r2 * 100))
    cat("\n")
  }
  
  for (i in 1:min(5, nrow(best_results))) {
    result <- best_results[i, ]
    cat(sprintf("%d. %s\n", i, result$prs_column))
    cat(sprintf("   P-value: %s\n", format(result$p_value, scientific = TRUE)))
    cat(sprintf("   OR: %.3f (95%% CI: %.3f-%.3f)\n", 
               result$OR, result$OR_lower_95CI, result$OR_upper_95CI))
    
    # Use incremental R² if available, otherwise fall back to total R²
    if (!is.na(result$r2_incremental_prs)) {
      cat(sprintf("   PC-only R²: %.2f%%, Full R²: %.2f%%, Incremental PRS R²: %.2f%%\n", 
                 result$r2_pc_only * 100, result$r2_full * 100, result$r2_incremental_prs * 100))
    } else {
      # Fallback for older data format
      additional_r2 <- result$nagelkerke_r2 - baseline_r2
      cat(sprintf("   Total R²: %.2f%% (Additional R² from PRS: %.2f%%)\n", 
                 result$nagelkerke_r2 * 100, additional_r2 * 100))
    }
    cat(sprintf("   Beta: %.4f (SE: %.4f)\n\n", result$beta, result$se))
  }
}

cat("DETAILED RESULTS BY PRS\n")
cat("=================================================================\n\n")

for (prs_col in prs_cols) {
  prs_results <- all_results %>% filter(prs_column == prs_col)
  if (nrow(prs_results) > 0) {
    cat("PRS Column:", prs_col, "\n")
    cat("-----------------------------------\n")
    
    for (i in 1:nrow(prs_results)) {
      result <- prs_results[i, ]
      pc_label <- paste0(result$pc_count, " PCs")
      
      if (!is.na(result$p_value)) {
        sig_status <- ifelse(result$p_value < 0.05, "SIGNIFICANT", "Not significant")
        cat("  ", pc_label, ":", sig_status, "(p =", format(result$p_value, scientific = TRUE), ")\n")
        cat("    Beta:", round(result$beta, 6), "(SE =", round(result$se, 6), ")\n")
        cat("    OR:", round(result$OR, 3), "(95% CI:", round(result$OR_lower_95CI, 3), "-", round(result$OR_upper_95CI, 3), ")\n")
        if (!is.na(result$r2_pc_only)) {
          cat("    PC-only R²:", round(result$r2_pc_only * 100, 2), "%,")
          cat(" Full R²:", round(result$r2_full * 100, 2), "%,")
          cat(" Incremental PRS R²:", round(result$r2_incremental_prs * 100, 2), "%\n")
        } else {
          cat("    Incremental PRS R²:", round(result$nagelkerke_r2 * 100, 2), "%\n")
        }
      } else {
        cat("  ", pc_label, ": FAILED\n")
        if ("error" %in% names(result)) {
          cat("    Error:", result$error, "\n")
        }
      }
    }
    cat("\n")
  }
}

cat("INTERPRETATION NOTES\n")
cat("=================================================================\n")
cat("- Beta: Log-odds change per unit increase in PRS\n")
cat("- OR: Odds ratio (exp(beta))\n")
cat("- PC-only R²: Variance explained by population structure alone\n")
cat("- Full R²: Total variance explained by PRS + population structure\n")
cat("- Incremental PRS R²: True PRS contribution = Full R² - PC-only R²\n")
cat("- IMPORTANT: Use incremental R² for reporting PRS performance\n")
cat("- Incremental R² should remain stable across different PC counts\n")
cat("- Use", max_pcs, "PC model for publication to control population structure\n\n")

sink()
cat("✓ Summary report saved to:", summary_file, "\n")

# Create a comparison file showing best PRS per threshold
comparison_file <- file.path(output_dir, paste0(analysis_name, "_prs_comparison.csv"))
if (nrow(best_results) > 0) {
  prs_comparison <- best_results %>%
    select(prs_column, p_value, OR, OR_lower_95CI, OR_upper_95CI, nagelkerke_r2, beta, se) %>%
    mutate(rank = row_number()) %>%
    select(rank, everything())
  
  fwrite(prs_comparison, comparison_file)
  cat("✓ PRS comparison saved to:", comparison_file, "\n\n")
}

# Display final results  
cat("=================================================================\n")
cat("FINAL RESULTS - MULTI-PRS ANALYSIS\n")
cat("=================================================================\n\n")

cat("Total PRS columns analyzed:", length(unique(all_results$prs_column)), "\n")
cat("Total models run:", nrow(all_results), "\n\n")

# Show top 3 PRS results
if (nrow(best_results) > 0) {
  cat("TOP PERFORMING PRS (", max_pcs, " PC models):\n")
  for (i in 1:min(3, nrow(best_results))) {
    result <- best_results[i, ]
    # Use incremental R² for display
    display_r2 <- if (!is.na(result$r2_incremental_prs)) result$r2_incremental_prs else result$nagelkerke_r2
    cat(sprintf("  %d. %s: p=%s, OR=%.3f, Incremental PRS R²=%.1f%%\n", 
               i, result$prs_column, format(result$p_value, scientific = TRUE), 
               result$OR, display_r2 * 100))
  }
  cat("\n")
}

cat("Detailed results table:\n")
print(all_results)

cat("\n=================================================================\n")
cat("MULTI-PRS ANALYSIS COMPLETE!\n")
cat("=================================================================\n\n")

cat("Output files:\n")
cat("  - Detailed results:", results_file, "\n")
cat("  - Summary report:", summary_file, "\n")
if (exists("comparison_file")) {
  cat("  - PRS comparison:", comparison_file, "\n")
}