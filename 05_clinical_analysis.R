# HER2 Integration Project - Session 5: Clinical Validation & Survival Analysis (FIXED)
# Author: Research Team
# Date: January 2025
# Objective: Clinical validation using REAL survival data and Bayesian scores

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(gridExtra)
  library(jsonlite)
  library(RColorBrewer)
  library(tidyr)
})

# Ensure dplyr functions are available
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("dplyr package is required but not available")
}

# Test dplyr functionality
test_df <- data.frame(a = 1:3, b = 4:6)
test_result <- test_df %>% dplyr::select(a)
if (nrow(test_result) != 3) {
  stop("dplyr is not functioning correctly")
}

# Session metadata initialization
session_metadata <- list(
  session_id = 5,
  date = Sys.Date(),
  objective = "Clinical validation and survival analysis using real TCGA data",
  start_time = Sys.time(),
  data_source = "modeling_results.rds with real Bayesian scores"
)

cat("=== HER2 Integration Project - Session 5 (FIXED) ===\n")
cat("Objective: Clinical validation using REAL modeling results\n")
cat("Start time:", as.character(session_metadata$start_time), "\n\n")

# Create output directories
dirs_to_create <- c("figures/clinical", "tables/clinical", "results/clinical")
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# =============================================================================
# 1. LOAD AND MERGE REAL DATA FROM MULTIPLE SOURCES
# =============================================================================

cat("1. Loading and merging real data from multiple sources...\n")

# Load Bayesian scores (your main results!)
if (!file.exists("data/results/latent_her2_scores.rds")) {
  stop("ERROR: latent_her2_scores.rds not found. Please run Sessions 3-4 first.")
}

latent_scores <- readRDS("data/results/latent_her2_scores.rds")
cat("✓ Bayesian scores loaded\n")
cat(sprintf("  - %d patients with latent HER2 scores\n", nrow(latent_scores)))

# Load clinical data for survival information
if (!file.exists("data/processed/tcga_clinical.rds")) {
  stop("ERROR: tcga_clinical.rds not found. Please run Session 1 first.")
}

clinical_data <- readRDS("data/processed/tcga_clinical.rds")
cat("✓ Clinical data loaded\n")
cat(sprintf("  - %d patients in clinical dataset\n", nrow(clinical_data)))

# Check clinical data columns for survival
clinical_cols <- names(clinical_data)
survival_cols <- clinical_cols[grepl("survival|time|event|death|status|vital", clinical_cols, ignore.case = TRUE)]
cat("  - Survival-related columns:", paste(survival_cols, collapse = ", "), "\n")

# Load harmonized platform data  
if (!file.exists("data/processed/harmonized_dataset.rds")) {
  stop("ERROR: harmonized_dataset.rds not found. Please run Session 2 first.")
}

harmonized_data <- readRDS("data/processed/harmonized_dataset.rds")
cat("✓ Harmonized platform data loaded\n")
cat(sprintf("  - %d patients in harmonized dataset\n", nrow(harmonized_data)))

# =============================================================================
# 2. MERGE ALL DATA SOURCES
# =============================================================================

cat("\n2. Merging data sources for survival analysis...\n")

# Start with Bayesian scores as the base - use explicit dplyr namespace
survival_data <- latent_scores %>%
  dplyr::rename(
    zeta_continuous = latent_her2_score,
    her2_traditional = her2_ihc_status
  ) %>%
  # Add platform data
  dplyr::left_join(
    harmonized_data %>% dplyr::select(patient_id, dplyr::any_of(c("age", "her2_status", "her2_binary", 
                                                                  "erbb2_rna_log2", "erbb2_cnv", "erbb2_rna", "erbb2_cnv_log2"))),
    by = "patient_id"
  ) %>%
  # Add clinical survival data
  dplyr::left_join(
    clinical_data %>% dplyr::select(patient_id, dplyr::any_of(survival_cols)),
    by = "patient_id"
  )

cat(sprintf("✓ Data merged: %d patients\n", nrow(survival_data)))

# Create survival variables from available clinical data
cat("Checking for survival data in merged dataset...\n")
merged_cols <- names(survival_data)
survival_related <- merged_cols[grepl("survival|time|event|death|status|vital", merged_cols, ignore.case = TRUE)]
cat("Available survival-related columns:", paste(survival_related, collapse = ", "), "\n")

# Check specifically for TCGA standard survival columns
tcga_survival_cols <- c("days_to_death", "days_to_last_follow_up", "vital_status", "survival_time", "event")
existing_survival_cols <- tcga_survival_cols[tcga_survival_cols %in% merged_cols]
cat("TCGA survival columns found:", paste(existing_survival_cols, collapse = ", "), "\n")

# Initialize survival variables as NA
survival_data$os_time_years <- NA_real_
survival_data$os_event <- NA_real_

# Only create survival variables if we have the necessary columns
if ("survival_time" %in% merged_cols) {
  cat("Using survival_time column\n")
  survival_data$os_time_years <- ifelse(!is.na(survival_data$survival_time), 
                                        pmax(as.numeric(survival_data$survival_time), 0.01), 
                                        NA_real_)
}

if ("days_to_death" %in% merged_cols && "days_to_last_follow_up" %in% merged_cols) {
  cat("Using days_to_death and days_to_last_follow_up columns\n")
  survival_data$os_time_years <- ifelse(
    !is.na(survival_data$days_to_death) | !is.na(survival_data$days_to_last_follow_up),
    pmax(ifelse(!is.na(survival_data$days_to_death), 
                as.numeric(survival_data$days_to_death), 
                as.numeric(survival_data$days_to_last_follow_up)) / 365.25, 0.01),
    survival_data$os_time_years  # Keep existing value if both are NA
  )
}

if ("days_to_death" %in% merged_cols && !"days_to_last_follow_up" %in% merged_cols) {
  cat("Using days_to_death column only\n")
  survival_data$os_time_years <- ifelse(!is.na(survival_data$days_to_death),
                                        pmax(as.numeric(survival_data$days_to_death) / 365.25, 0.01),
                                        survival_data$os_time_years)
}

if ("days_to_last_follow_up" %in% merged_cols && !"days_to_death" %in% merged_cols) {
  cat("Using days_to_last_follow_up column only\n")
  survival_data$os_time_years <- ifelse(!is.na(survival_data$days_to_last_follow_up),
                                        pmax(as.numeric(survival_data$days_to_last_follow_up) / 365.25, 0.01),
                                        survival_data$os_time_years)
}

# Create event indicators
if ("event" %in% merged_cols) {
  cat("Using event column\n")
  survival_data$os_event <- ifelse(!is.na(survival_data$event), 
                                   as.numeric(survival_data$event), 
                                   NA_real_)
}

if ("vital_status" %in% merged_cols) {
  cat("Using vital_status column\n")
  survival_data$os_event <- ifelse(!is.na(survival_data$vital_status), 
                                   as.numeric(survival_data$vital_status == "Dead"), 
                                   survival_data$os_event)
}

if ("days_to_death" %in% merged_cols) {
  cat("Using days_to_death for event indicator\n")
  survival_data$os_event <- ifelse(!is.na(survival_data$days_to_death), 
                                   1,  # Death occurred
                                   survival_data$os_event)
}

if ("days_to_last_follow_up" %in% merged_cols && !"days_to_death" %in% merged_cols) {
  cat("Using days_to_last_follow_up for censoring\n")
  survival_data$os_event <- ifelse(!is.na(survival_data$days_to_last_follow_up) & is.na(survival_data$os_event), 
                                   0,  # Censored
                                   survival_data$os_event)
}

# Check how much real survival data we have
patients_with_survival <- sum(!is.na(survival_data$os_time_years) & !is.na(survival_data$os_event))
cat(sprintf("✓ Real TCGA survival data available for %d patients\n", patients_with_survival))

if (patients_with_survival < 50) {
  cat("ERROR: Insufficient real survival data for analysis\n")
  cat("Please ensure TCGA clinical data contains survival information\n")
  cat("Available survival columns in clinical data:\n")
  print(survival_cols)
  stop("Cannot proceed without real survival data - no simulation allowed")
}

# Final data cleaning - use more permissive filtering
survival_data <- survival_data %>%
  filter(
    !is.na(os_time_years),           # Must have survival time
    !is.na(os_event),                # Must have event status
    os_time_years > 0,               # Positive survival time
    os_time_years <= 20,             # Reasonable upper bound
    !is.na(patient_id)               # Must have patient ID
  )

cat("✓ Survival data prepared using real modeling results\n")
cat(sprintf("   - Final sample size: %d patients\n", nrow(survival_data)))
cat(sprintf("   - Median follow-up: %.1f years\n", median(survival_data$os_time_years, na.rm = TRUE)))
cat(sprintf("   - Events (deaths): %d (%.1f%%)\n", 
            sum(survival_data$os_event, na.rm = TRUE), 
            mean(survival_data$os_event, na.rm = TRUE) * 100))

# Check data availability for different methods
cat("\nData availability by method:\n")
cat(sprintf("   - Traditional HER2: %d patients\n", sum(!is.na(survival_data$her2_traditional) & survival_data$her2_traditional != "Unknown")))
cat(sprintf("   - RNA quartiles: %d patients\n", sum(!is.na(survival_data$rna_quartile))))
cat(sprintf("   - CNV categories: %d patients\n", sum(!is.na(survival_data$cnv_category))))
cat(sprintf("   - Bayesian scores: %d patients\n", sum(!is.na(survival_data$zeta_continuous))))

# =============================================================================
# 3. SURVIVAL ANALYSIS BY DIFFERENT HER2 METHODS
# =============================================================================

cat("\n3. Performing survival analysis by different HER2 methods...\n")

# Function to perform comprehensive survival analysis
perform_survival_analysis <- function(data, grouping_var, method_name) {
  
  # Filter for complete cases for this specific analysis
  analysis_data <- data %>%
    filter(!is.na(.data[[grouping_var]]), .data[[grouping_var]] != "Unknown")
  
  if (nrow(analysis_data) < 10) {
    cat(sprintf("   Warning: Insufficient data for %s (%d patients)\n", method_name, nrow(analysis_data)))
    return(list(
      method = method_name,
      n_total = nrow(analysis_data),
      n_events = 0,
      error = "Insufficient data"
    ))
  }
  
  # Create survival object
  surv_obj <- Surv(analysis_data$os_time_years, analysis_data$os_event)
  
  # Fit survival curves
  formula_str <- paste0("surv_obj ~ ", grouping_var)
  fit <- survfit(as.formula(formula_str), data = analysis_data)
  
  # Log-rank test
  logrank_test <- survdiff(as.formula(formula_str), data = analysis_data)
  logrank_p <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
  
  # Cox proportional hazards model
  tryCatch({
    cox_formula <- as.formula(paste0("Surv(os_time_years, os_event) ~ ", grouping_var))
    cox_model <- coxph(cox_formula, data = analysis_data)
    c_index <- concordance(cox_model)$concordance
    
    # Extract hazard ratios
    cox_summary <- summary(cox_model)
    hr_table <- data.frame(
      variable = rownames(cox_summary$coefficients),
      hr = cox_summary$coefficients[, "exp(coef)"],
      hr_lower = cox_summary$conf.int[, "lower .95"],
      hr_upper = cox_summary$conf.int[, "upper .95"], 
      p_value = cox_summary$coefficients[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("   Warning: Cox model failed for %s: %s\n", method_name, e$message))
    cox_model <- NULL
    c_index <- NA
    hr_table <- data.frame()
  })
  
  # 5-year survival rates
  surv_5yr <- tryCatch({
    summary(fit, times = 5)
  }, error = function(e) NULL)
  
  surv_rates <- if (!is.null(surv_5yr) && length(surv_5yr$time) > 0) {
    data.frame(
      group = surv_5yr$strata,
      survival_5yr = surv_5yr$surv,
      survival_5yr_lower = surv_5yr$lower,
      survival_5yr_upper = surv_5yr$upper,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }
  
  return(list(
    method = method_name,
    fit = fit,
    logrank_p = logrank_p,
    cox_model = cox_model,
    c_index = c_index,
    hr_table = hr_table,
    surv_rates = surv_rates,
    n_events = sum(analysis_data$os_event),
    n_total = nrow(analysis_data),
    analysis_data = analysis_data
  ))
}

# Define survival analysis methods using real data
survival_methods <- list(
  list(var = "her2_traditional", name = "Traditional HER2"),
  list(var = "rna_quartile_label", name = "RNA Quartiles"),
  list(var = "cnv_category_label", name = "CNV Categories"), 
  list(var = "zeta_tertile_label", name = "Integrated Bayesian Score")
)

# Perform survival analysis for different methods
os_results <- list()
for (method in survival_methods) {
  cat(sprintf("   Analyzing %s...\n", method$name))
  os_results[[method$name]] <- perform_survival_analysis(
    survival_data, method$var, method$name
  )
  
  # Print summary if successful
  if (!is.null(os_results[[method$name]]$fit)) {
    result <- os_results[[method$name]]
    cat(sprintf("     ✓ %d patients, %d events, log-rank p=%.3f\n", 
                result$n_total, result$n_events, result$logrank_p))
  }
}

# =============================================================================
# 4. CREATE KAPLAN-MEIER PLOTS
# =============================================================================

cat("\n4. Creating Kaplan-Meier survival plots...\n")

# Function to create publication-quality KM plot
create_km_plot <- function(analysis_result, colors = NULL) {
  
  if (is.null(analysis_result$fit) || is.null(analysis_result$analysis_data)) {
    return(ggplot() + theme_void() + labs(title = paste("No data available for", analysis_result$method)))
  }
  
  # Default colors if not provided
  if (is.null(colors)) {
    n_groups <- length(unique(analysis_result$analysis_data[[names(analysis_result$analysis_data)[length(names(analysis_result$analysis_data))]]]))
    colors <- RColorBrewer::brewer.pal(min(n_groups, 8), "Set2")
  }
  
  # Create the plot
  km_plot <- ggsurvplot(
    analysis_result$fit,
    data = analysis_result$analysis_data,
    pval = TRUE,
    pval.coord = c(0.02, 0.02),
    conf.int = TRUE,
    conf.int.alpha = 0.1,
    risk.table = TRUE,
    risk.table.height = 0.3,
    ncensor.plot = FALSE,
    surv.median.line = "hv",
    ggtheme = theme_bw(base_size = 12),
    palette = colors,
    title = paste("Overall Survival -", analysis_result$method),
    subtitle = sprintf("Log-rank p = %.3f, C-index = %.3f", 
                       analysis_result$logrank_p, 
                       ifelse(is.na(analysis_result$c_index), 0, analysis_result$c_index)),
    xlab = "Time (years)",
    ylab = "Overall Survival Probability",
    legend.title = analysis_result$method,
    font.main = c(14, "bold"),
    font.submain = c(12, "plain"),
    font.x = c(12, "bold"),
    font.y = c(12, "bold"),
    font.legend = c(10, "plain")
  )
  
  return(km_plot)
}

# Create and save KM plots for each method
km_plots <- list()

for (method_name in names(os_results)) {
  if (!is.null(os_results[[method_name]]$fit)) {
    cat(sprintf("   Creating KM plot for %s...\n", method_name))
    
    # Create plot
    km_plots[[method_name]] <- create_km_plot(os_results[[method_name]])
    
    # Save plot
    method_clean <- gsub("[^A-Za-z0-9]", "_", method_name)
    filename <- sprintf("figures/clinical/km_survival_%s.png", method_clean)
    
    ggsave(
      filename = filename,
      plot = print(km_plots[[method_name]]),
      width = 12, height = 10, dpi = 300
    )
    
    cat(sprintf("     ✓ Saved: %s\n", filename))
  } else {
    cat(sprintf("   ⚠ Skipping %s - insufficient data\n", method_name))
  }
}

# =============================================================================
# 5. COMPARE MODEL PERFORMANCE
# =============================================================================

cat("\n5. Comparing model performance...\n")

# Create performance comparison table
performance_comparison <- data.frame()

for (method_name in names(os_results)) {
  result <- os_results[[method_name]]
  
  if (!is.null(result$fit)) {
    performance_comparison <- rbind(performance_comparison, data.frame(
      Method = result$method,
      N_Patients = result$n_total,
      N_Events = result$n_events,
      Event_Rate = round(result$n_events / result$n_total * 100, 1),
      LogRank_P = round(result$logrank_p, 4),
      C_Index = round(ifelse(is.na(result$c_index), 0, result$c_index), 3),
      Significant = ifelse(result$logrank_p < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
}

# Sort by C-index (descending)
if (nrow(performance_comparison) > 0) {
  performance_comparison <- performance_comparison[order(-performance_comparison$C_Index), ]
  
  cat("✓ Performance comparison completed\n")
  cat("\nModel Performance Summary:\n")
  print(performance_comparison)
  
  # Save performance comparison
  write.csv(performance_comparison, "tables/clinical/performance_comparison.csv", row.names = FALSE)
  
} else {
  cat("⚠ No successful analyses to compare\n")
}

# =============================================================================
# 6. RISK STRATIFICATION ANALYSIS
# =============================================================================

cat("\n6. Performing risk stratification analysis...\n")

# Create comprehensive risk analysis using Bayesian scores
if (sum(!is.na(survival_data$zeta_continuous)) > 20) {
  
  # Continuous risk analysis
  risk_data <- survival_data %>%
    filter(!is.na(zeta_continuous)) %>%
    mutate(
      zeta_binary = ifelse(zeta_continuous > median(zeta_continuous, na.rm = TRUE), "High", "Low"),
      zeta_quartile = ntile(zeta_continuous, 4),
      zeta_quartile_label = paste0("Q", zeta_quartile)
    )
  
  # Binary risk stratification
  binary_analysis <- perform_survival_analysis(risk_data, "zeta_binary", "Binary Risk (Median Split)")
  
  if (!is.null(binary_analysis$fit)) {
    # Create risk stratification plot
    risk_km_plot <- create_km_plot(binary_analysis, colors = c("#2E86AB", "#E63946"))
    
    # Save risk plot
    ggsave(
      filename = "figures/clinical/risk_stratification_binary.png",
      plot = print(risk_km_plot),
      width = 12, height = 10, dpi = 300
    )
    
    cat("   ✓ Binary risk stratification completed\n")
    cat(sprintf("     - High risk: %d patients\n", sum(risk_data$zeta_binary == "High")))
    cat(sprintf("     - Low risk: %d patients\n", sum(risk_data$zeta_binary == "Low")))
    cat(sprintf("     - Log-rank p: %.4f\n", binary_analysis$logrank_p))
  }
  
  # Quartile risk analysis
  quartile_analysis <- perform_survival_analysis(risk_data, "zeta_quartile_label", "Quartile Risk Stratification")
  
  if (!is.null(quartile_analysis$fit)) {
    # Create quartile plot
    quartile_km_plot <- create_km_plot(quartile_analysis)
    
    # Save quartile plot
    ggsave(
      filename = "figures/clinical/risk_stratification_quartiles.png",
      plot = print(quartile_km_plot),
      width = 12, height = 10, dpi = 300
    )
    
    cat("   ✓ Quartile risk stratification completed\n")
    cat(sprintf("     - Log-rank p: %.4f\n", quartile_analysis$logrank_p))
  }
  
} else {
  cat("   ⚠ Insufficient Bayesian score data for risk stratification\n")
}

# =============================================================================
# 7. CLINICAL VALIDATION SUMMARY
# =============================================================================

cat("\n7. Creating clinical validation summary...\n")

# Compile all results
clinical_validation_results <- list(
  study_summary = list(
    total_patients = nrow(survival_data),
    median_followup = median(survival_data$os_time_years, na.rm = TRUE),
    total_events = sum(survival_data$os_event, na.rm = TRUE),
    event_rate = mean(survival_data$os_event, na.rm = TRUE)
  ),
  survival_results = os_results,
  performance_comparison = performance_comparison,
  data_source = "Real TCGA modeling results with Bayesian scores"
)

# Save complete results
saveRDS(clinical_validation_results, "results/clinical/clinical_validation_results.rds")

# Update session metadata
session_metadata$end_time <- Sys.time()
session_metadata$duration_minutes <- as.numeric(difftime(
  session_metadata$end_time, 
  session_metadata$start_time, 
  units = "mins"
))

session_metadata$results_summary <- list(
  patients_analyzed = nrow(survival_data),
  methods_compared = length(os_results),
  successful_analyses = sum(sapply(os_results, function(x) !is.null(x$fit))),
  best_method = if(nrow(performance_comparison) > 0) performance_comparison$Method[1] else "None",
  significant_results = if(nrow(performance_comparison) > 0) sum(performance_comparison$Significant == "Yes") else 0
)

# Save session metadata
write_json(session_metadata, "metadata/session5_metadata.json", pretty = TRUE)

# =============================================================================
# 8. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("=== SESSION 5: CLINICAL VALIDATION COMPLETED ===\n")
cat("✓ Used REAL modeling results with Bayesian scores\n")
cat("✓ Survival analysis with multiple HER2 assessment methods\n")
cat("✓ Kaplan-Meier curves with statistical testing\n")
cat("✓ Performance comparison across methods\n")
cat("✓ Risk stratification analysis\n")

cat(sprintf("\nSTUDY SUMMARY:\n"))
cat(sprintf("- Patients analyzed: %d\n", nrow(survival_data)))
cat(sprintf("- Median follow-up: %.1f years\n", median(survival_data$os_time_years, na.rm = TRUE)))
cat(sprintf("- Total events: %d (%.1f%%)\n", 
            sum(survival_data$os_event, na.rm = TRUE),
            mean(survival_data$os_event, na.rm = TRUE) * 100))

if (nrow(performance_comparison) > 0) {
  cat(sprintf("\nBEST PERFORMING METHOD:\n"))
  best_method <- performance_comparison[1, ]
  cat(sprintf("- Method: %s\n", best_method$Method))
  cat(sprintf("- C-index: %.3f\n", best_method$C_Index))
  cat(sprintf("- Log-rank p: %.4f\n", best_method$LogRank_P))
  cat(sprintf("- Significant: %s\n", best_method$Significant))
}

cat("\nFILES CREATED:\n")
cat("==============\n")

# List created files
figures <- list.files("figures/clinical", pattern = "\\.png$", full.names = FALSE)
if (length(figures) > 0) {
  cat("📊 Survival Plots:\n")
  for (fig in figures) {
    cat(sprintf("  - figures/clinical/%s\n", fig))
  }
}

tables <- list.files("tables/clinical", pattern = "\\.csv$", full.names = FALSE)
if (length(tables) > 0) {
  cat("\n📋 Tables:\n")
  for (tbl in tables) {
    cat(sprintf("  - tables/clinical/%s\n", tbl))
  }
}

cat("\n🎯 SESSION 5 COMPLETED SUCCESSFULLY!\n")
cat("Ready for Session 6: Publication Figures & Tables\n")

# Clean up environment
rm(list = ls()[!ls() %in% c("session_metadata")])
gc()

cat("\n=== Session 5 Complete ===\n")
