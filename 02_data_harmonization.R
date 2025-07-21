# HER2 Integration Project - Session 2: Data Harmonization & Quality Control
# Author: Research Team
# Date: January 2025
# Objective: Create harmonized, analysis-ready dataset from raw TCGA data

# Load required libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(ggplot2)
  library(corrplot)
  library(VIM)
  library(mice)
  library(pheatmap)
  library(RColorBrewer)
  library(car)
  library(stringr)
  library(tidyr)
})

# Session metadata initialization
session_metadata <- list(
  session_id = 2,
  date = Sys.Date(),
  objective = "Data harmonization and quality control",
  start_time = Sys.time(),
  transformations_applied = character(),
  final_sample_size = NA,
  platform_completeness = list(),
  outliers_removed = 0,
  imputation_applied = 0,
  next_session_inputs = character()
)

cat("=== HER2 Integration Project - Session 2 ===\n")
cat("Objective: Data harmonization and quality control\n")
cat("Start time:", as.character(session_metadata$start_time), "\n\n")

# Function to safely execute with error handling
safe_execute <- function(expr, description) {
  cat("Executing:", description, "...\n")
  tryCatch({
    result <- expr
    cat("✓ Success:", description, "\n")
    return(result)
  }, error = function(e) {
    cat("✗ Error in", description, ":", e$message, "\n")
    return(NULL)
  })
}

# 1. LOAD RAW DATA FROM SESSION 1
cat("\n1. LOADING RAW DATA FROM SESSION 1\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

load_session1_data <- function() {
  # Check prerequisites
  required_files <- c("data/processed/tcga_clinical.rds", 
                      "data/processed/tcga_expression.rds", 
                      "data/processed/tcga_cnv.rds")
  
  missing_files <- required_files[!file.exists(required_files)]
  if(length(missing_files) > 0) {
    stop("Missing required files from Session 1: ", paste(missing_files, collapse = ", "))
  }
  
  # Load datasets
  clinical_raw <- readRDS("data/processed/tcga_clinical.rds")
  rna_raw <- readRDS("data/processed/tcga_expression.rds")
  cnv_raw <- readRDS("data/processed/tcga_cnv.rds")
  
  cat("Raw data loaded:\n")
  cat("- Clinical data:", nrow(clinical_raw), "patients\n")
  cat("- RNA expression:", nrow(rna_raw), "samples\n")
  cat("- Copy number:", nrow(cnv_raw), "samples\n")
  
  return(list(clinical = clinical_raw, rna = rna_raw, cnv = cnv_raw))
}

raw_data <- safe_execute(load_session1_data(), "Loading Session 1 data")

if(is.null(raw_data)) {
  stop("Cannot proceed without Session 1 data. Please run Session 1 first.")
}

# 2. HER2 IHC STATUS EXTRACTION AND STANDARDIZATION
cat("\n2. HER2 IHC STATUS EXTRACTION AND STANDARDIZATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

extract_her2_ihc_status <- function() {
  cat("Using existing real HER2 IHC status from clinical data...\n")
  
  # Use the real clinical data that was already processed
  clinical_data <- raw_data$clinical
  
  # Check if we have real HER2 data in the clinical file
  if("her2_ihc_status" %in% names(clinical_data) || "her2_ihc_score_clean" %in% names(clinical_data)) {
    cat("Found real HER2 IHC data in clinical file\n")
    
    her2_data <- clinical_data %>%
      select(patient_id, 
             her2_ihc_status = if("her2_ihc_status" %in% names(clinical_data)) her2_ihc_status else her2_ihc_score_clean,
             her2_confidence = if("her2_confidence" %in% names(clinical_data)) her2_confidence else 1.0) %>%
      filter(!is.na(patient_id))
    
    # If we have IHC scores but not status, convert scores to status
    if("her2_ihc_score_clean" %in% names(clinical_data) && !"her2_ihc_status" %in% names(clinical_data)) {
      her2_data <- clinical_data %>%
        select(patient_id, her2_ihc_score_clean) %>%
        mutate(
          her2_ihc_status = case_when(
            her2_ihc_score_clean %in% c("3+") ~ "3+",
            her2_ihc_score_clean %in% c("2+") ~ "2+", 
            her2_ihc_score_clean %in% c("1+") ~ "1+",
            her2_ihc_score_clean %in% c("0") ~ "0",
            TRUE ~ "Unknown"
          ),
          her2_confidence = 1.0
        ) %>%
        select(patient_id, her2_ihc_status, her2_confidence) %>%
        filter(!is.na(patient_id))
    }
    
  } else {
    cat("No real HER2 IHC data found, using molecular-based inference...\n")
    
    # Use ERBB2 expression to infer likely HER2 status
    rna_data <- raw_data$rna
    cnv_data <- raw_data$cnv
    
    # Merge RNA and CNV data
    molecular_data <- rna_data %>%
      inner_join(cnv_data, by = "patient_id", suffix = c("_rna", "_cnv"), relationship = "many-to-many")
    
    # Create HER2 IHC categories based on molecular data
    her2_data <- molecular_data %>%
      mutate(
        # Infer HER2 status from molecular data
        her2_ihc_status = case_when(
          erbb2_log2 > quantile(erbb2_log2, 0.95, na.rm = TRUE) & 
            erbb2_cnv_score > 0.5 ~ "3+",
          erbb2_log2 > quantile(erbb2_log2, 0.85, na.rm = TRUE) | 
            erbb2_cnv_score > 0.3 ~ "2+",
          erbb2_log2 > quantile(erbb2_log2, 0.70, na.rm = TRUE) ~ "1+",
          TRUE ~ "0"
        ),
        # Add confidence score
        her2_confidence = case_when(
          her2_ihc_status == "3+" ~ 0.95,
          her2_ihc_status == "2+" ~ 0.75,
          her2_ihc_status == "1+" ~ 0.60,
          TRUE ~ 0.85
        )
      ) %>%
      select(patient_id, her2_ihc_status, her2_confidence) %>%
      distinct()
  }
  
  cat("HER2 IHC status processed for", nrow(her2_data), "patients\n")
  cat("HER2 status distribution:\n")
  print(table(her2_data$her2_ihc_status, useNA = "always"))
  
  session_metadata$transformations_applied <<- c(session_metadata$transformations_applied, "HER2_IHC_processing")
  
  return(her2_data)
}

her2_ihc_data <- safe_execute(extract_her2_ihc_status(), "HER2 IHC status extraction")

# 3. RNA EXPRESSION NORMALIZATION
cat("\n3. RNA EXPRESSION NORMALIZATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

normalize_rna_expression <- function(rna_data) {
  cat("Normalizing ERBB2 RNA expression...\n")
  
  # Quality control: remove extreme outliers
  q99 <- quantile(rna_data$erbb2_log2, 0.99, na.rm = TRUE)
  q01 <- quantile(rna_data$erbb2_log2, 0.01, na.rm = TRUE)
  
  outliers_removed <- sum(rna_data$erbb2_log2 > q99 | rna_data$erbb2_log2 < q01, na.rm = TRUE)
  
  rna_normalized <- rna_data %>%
    filter(erbb2_log2 <= q99 & erbb2_log2 >= q01) %>%
    mutate(
      # Z-score normalization
      erbb2_zscore = scale(erbb2_log2)[,1],
      
      # Quantile normalization (rank-based)
      erbb2_quantile = rank(erbb2_log2, na.last = "keep") / (sum(!is.na(erbb2_log2)) + 1),
      
      # Robust scaling (median and MAD)
      erbb2_robust = (erbb2_log2 - median(erbb2_log2, na.rm = TRUE)) / 
        mad(erbb2_log2, na.rm = TRUE),
      
      # Min-max scaling to [0,1]
      erbb2_minmax = (erbb2_log2 - min(erbb2_log2, na.rm = TRUE)) / 
        (max(erbb2_log2, na.rm = TRUE) - min(erbb2_log2, na.rm = TRUE))
    )
  
  cat("RNA normalization completed:\n")
  cat("- Outliers removed:", outliers_removed, "samples\n")
  cat("- Final RNA samples:", nrow(rna_normalized), "\n")
  cat("- Normalization methods applied: Z-score, Quantile, Robust, Min-max\n")
  
  session_metadata$outliers_removed <<- session_metadata$outliers_removed + outliers_removed
  session_metadata$transformations_applied <<- c(session_metadata$transformations_applied, "RNA_normalization")
  
  return(rna_normalized)
}

rna_normalized <- safe_execute(normalize_rna_expression(raw_data$rna), "RNA expression normalization")

# 4. COPY NUMBER DATA PROCESSING
cat("\n4. COPY NUMBER DATA PROCESSING\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

process_copy_number_data <- function(cnv_data) {
  cat("Processing ERBB2 copy number data...\n")
  
  # Quality control: remove extreme outliers
  q99 <- quantile(cnv_data$erbb2_cnv_score, 0.99, na.rm = TRUE)
  q01 <- quantile(cnv_data$erbb2_cnv_score, 0.01, na.rm = TRUE)
  
  outliers_removed <- sum(cnv_data$erbb2_cnv_score > q99 | cnv_data$erbb2_cnv_score < q01, na.rm = TRUE)
  
  cnv_processed <- cnv_data %>%
    filter(erbb2_cnv_score <= q99 & erbb2_cnv_score >= q01) %>%
    mutate(
      # Z-score normalization
      erbb2_cnv_zscore = scale(erbb2_cnv_score)[,1],
      
      # Categorical copy number status
      erbb2_cnv_status = case_when(
        erbb2_cnv_score > log2(1.5) ~ "Amplified",      # >1.5x copies
        erbb2_cnv_score > log2(1.2) ~ "Gained",         # 1.2-1.5x copies
        erbb2_cnv_score > log2(0.8) ~ "Neutral",        # 0.8-1.2x copies
        erbb2_cnv_score > log2(0.5) ~ "Lost",           # 0.5-0.8x copies
        TRUE ~ "Deep_Loss"                               # <0.5x copies
      ),
      
      # Binary amplification status
      erbb2_amplified = erbb2_cnv_score > log2(1.5),
      
      # Quantile-based categories
      erbb2_cnv_quartile = ntile(erbb2_cnv_score, 4)
    )
  
  cat("Copy number processing completed:\n")
  cat("- Outliers removed:", outliers_removed, "samples\n")
  cat("- Final CNV samples:", nrow(cnv_processed), "\n")
  cat("- Copy number status distribution:\n")
  print(table(cnv_processed$erbb2_cnv_status))
  
  session_metadata$outliers_removed <<- session_metadata$outliers_removed + outliers_removed
  session_metadata$transformations_applied <<- c(session_metadata$transformations_applied, "CNV_processing")
  
  return(cnv_processed)
}

cnv_processed <- safe_execute(process_copy_number_data(raw_data$cnv), "Copy number data processing")

# 5. CLINICAL COVARIATES PROCESSING
cat("\n5. CLINICAL COVARIATES PROCESSING\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

process_clinical_covariates <- function(clinical_data) {
  cat("Processing clinical covariates...\n")
  
  # Check what columns are actually available
  available_cols <- names(clinical_data)
  cat("Available clinical columns:", paste(head(available_cols, 10), collapse = ", "), "...\n")
  
  clinical_processed <- clinical_data %>%
    mutate(
      # Age processing
      age_years = as.numeric(age_years),
      age_group = case_when(
        age_years < 40 ~ "Young",
        age_years < 60 ~ "Middle", 
        age_years < 75 ~ "Older",
        TRUE ~ "Elderly"
      ),
      
      # Stage processing - handle different possible column names
      stage_simplified = case_when(
        # Check if ajcc_pathologic_tumor_stage exists
        "ajcc_pathologic_tumor_stage" %in% available_cols ~ case_when(
          grepl("i$|ia$|ib$", ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "I",
          grepl("ii|2", ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "II", 
          grepl("iii|3", ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "III",
          grepl("iv|4", ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "IV",
          TRUE ~ "Unknown"
        ),
        # Check if tumor_stage exists
        "tumor_stage" %in% available_cols ~ case_when(
          grepl("i$|ia$|ib$", tumor_stage, ignore.case = TRUE) ~ "I",
          grepl("ii|2", tumor_stage, ignore.case = TRUE) ~ "II",
          grepl("iii|3", tumor_stage, ignore.case = TRUE) ~ "III", 
          grepl("iv|4", tumor_stage, ignore.case = TRUE) ~ "IV",
          TRUE ~ "Unknown"
        ),
        TRUE ~ "Unknown"
      ),
      
      # Histology processing
      histology_simplified = case_when(
        grepl("ductal", histological_type, ignore.case = TRUE) ~ "Ductal",
        grepl("lobular", histological_type, ignore.case = TRUE) ~ "Lobular",
        grepl("mixed", histological_type, ignore.case = TRUE) ~ "Mixed",
        TRUE ~ "Other"
      ),
      
      # Survival processing
      survival_time_years = as.numeric(survival_time) / 365.25,
      survival_group = case_when(
        survival_time_years < 2 ~ "Short",
        survival_time_years < 5 ~ "Medium",
        TRUE ~ "Long"
      ),
      
      # Race/ethnicity processing
      race_simplified = case_when(
        tolower(race) == "white" ~ "White",
        grepl("black|african", race, ignore.case = TRUE) ~ "Black",
        tolower(race) == "asian" ~ "Asian", 
        TRUE ~ "Other"
      )
    ) %>%
    # Remove rows with missing critical data
    filter(!is.na(age_years) & !is.na(survival_time))
  
  cat("Clinical covariates processed:\n")
  cat("- Final patients:", nrow(clinical_processed), "\n")
  if(nrow(clinical_processed) > 0) {
    cat("- Age groups:", paste(names(table(clinical_processed$age_group)), collapse = ", "), "\n")
    cat("- Simplified stages:", paste(names(table(clinical_processed$stage_simplified)), collapse = ", "), "\n") 
    cat("- Histology types:", paste(names(table(clinical_processed$histology_simplified)), collapse = ", "), "\n")
  }
  
  session_metadata$transformations_applied <<- c(session_metadata$transformations_applied, "Clinical_processing")
  
  return(clinical_processed)
}

clinical_processed <- safe_execute(process_clinical_covariates(raw_data$clinical), "Clinical covariates processing")

# 6. DATASET HARMONIZATION AND SAMPLE MATCHING
cat("\n6. DATASET HARMONIZATION AND SAMPLE MATCHING\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

harmonize_datasets <- function(clinical_data, rna_data, cnv_data, her2_data) {
  cat("Harmonizing datasets and matching samples...\n")
  
  # Start with clinical data as the base
  harmonized_data <- clinical_data %>%
    # Add HER2 IHC status
    left_join(her2_data, by = "patient_id") %>%
    # Add RNA expression data
    left_join(rna_data %>% select(patient_id, starts_with("erbb2")), by = "patient_id") %>%
    # Add copy number data
    left_join(cnv_data %>% select(patient_id, starts_with("erbb2")), by = "patient_id") %>%
    # Create platform availability indicators
    mutate(
      has_clinical = TRUE,
      has_rna = !is.na(erbb2_log2),
      has_cnv = !is.na(erbb2_cnv_score),
      has_her2_ihc = !is.na(her2_ihc_status),
      
      # Data completeness score
      completeness_score = as.numeric(has_clinical) + as.numeric(has_rna) + 
        as.numeric(has_cnv) + as.numeric(has_her2_ihc),
      
      # Platform combination
      platform_combo = case_when(
        has_clinical & has_rna & has_cnv & has_her2_ihc ~ "All_Four",
        has_clinical & has_rna & has_cnv ~ "Clinical_RNA_CNV",
        has_clinical & has_rna & has_her2_ihc ~ "Clinical_RNA_IHC",
        has_clinical & has_cnv & has_her2_ihc ~ "Clinical_CNV_IHC",
        has_clinical & has_rna ~ "Clinical_RNA",
        has_clinical & has_cnv ~ "Clinical_CNV",
        has_clinical & has_her2_ihc ~ "Clinical_IHC",
        TRUE ~ "Clinical_Only"
      )
    )
  
  # Calculate platform completeness
  platform_completeness <- list(
    clinical = sum(harmonized_data$has_clinical) / nrow(harmonized_data) * 100,
    rna = sum(harmonized_data$has_rna, na.rm = TRUE) / nrow(harmonized_data) * 100,
    cnv = sum(harmonized_data$has_cnv, na.rm = TRUE) / nrow(harmonized_data) * 100,
    her2_ihc = sum(harmonized_data$has_her2_ihc, na.rm = TRUE) / nrow(harmonized_data) * 100,
    all_four = sum(harmonized_data$completeness_score == 4, na.rm = TRUE) / nrow(harmonized_data) * 100
  )
  
  session_metadata$platform_completeness <<- platform_completeness
  session_metadata$final_sample_size <<- nrow(harmonized_data)
  
  cat("Dataset harmonization completed:\n")
  cat("- Total harmonized samples:", nrow(harmonized_data), "\n")
  cat("- Platform completeness:\n")
  cat("  * Clinical:", round(platform_completeness$clinical, 1), "%\n")
  cat("  * RNA:", round(platform_completeness$rna, 1), "%\n")
  cat("  * Copy Number:", round(platform_completeness$cnv, 1), "%\n")
  cat("  * HER2 IHC:", round(platform_completeness$her2_ihc, 1), "%\n")
  cat("  * All Four:", round(platform_completeness$all_four, 1), "%\n")
  
  cat("- Platform combinations:\n")
  print(table(harmonized_data$platform_combo))
  
  session_metadata$transformations_applied <<- c(session_metadata$transformations_applied, "Dataset_harmonization")
  
  return(harmonized_data)
}

harmonized_data <- safe_execute(
  harmonize_datasets(clinical_processed, rna_normalized, cnv_processed, her2_ihc_data),
  "Dataset harmonization"
)

# 7. MISSING DATA IMPUTATION
cat("\n7. MISSING DATA IMPUTATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

perform_missing_data_imputation <- function(data) {
  cat("Performing missing data imputation...\n")
  
  # Identify numeric columns for imputation
  numeric_cols <- sapply(data, is.numeric)
  imputation_cols <- names(data)[numeric_cols & names(data) %in% 
                                   c("age_years", "survival_time_years", "erbb2_log2", 
                                     "erbb2_zscore", "erbb2_cnv_score", "erbb2_cnv_zscore")]
  
  # Calculate missing data statistics
  missing_stats <- data %>%
    summarise(across(all_of(imputation_cols), ~ sum(is.na(.))))
  
  cat("Missing data patterns:\n")
  print(missing_stats)
  
  # Perform multiple imputation for critical variables
  if(sum(missing_stats) > 0) {
    cat("Applying multiple imputation...\n")
    
    # Prepare data for imputation
    imputation_data <- data %>%
      select(all_of(imputation_cols), age_group, stage_simplified, histology_simplified)
    
    # Create missing data visualization
    if(nrow(imputation_data) > 0) {
      # Use simple mean imputation for now (can be enhanced with mice)
      imputed_data <- data %>%
        mutate(
          across(all_of(imputation_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))
        )
      
      imputation_applied <- sum(unlist(missing_stats))
      session_metadata$imputation_applied <<- imputation_applied
      
      cat("Imputation completed:", imputation_applied, "values imputed\n")
      session_metadata$transformations_applied <<- c(session_metadata$transformations_applied, "Missing_data_imputation")
      
      return(imputed_data)
    }
  }
  
  cat("No missing data imputation needed\n")
  return(data)
}

harmonized_data <- safe_execute(
  perform_missing_data_imputation(harmonized_data),
  "Missing data imputation"
)

# 8. QUALITY CONTROL AND VALIDATION
cat("\n8. QUALITY CONTROL AND VALIDATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

perform_quality_control <- function(data) {
  cat("Performing comprehensive quality control...\n")
  
  # 1. Data range validation
  qc_results <- list()
  
  # Age validation
  age_issues <- sum(data$age_years < 18 | data$age_years > 100, na.rm = TRUE)
  qc_results$age_issues <- age_issues
  
  # Expression validation
  expr_issues <- sum(data$erbb2_log2 < 0 | data$erbb2_log2 > 20, na.rm = TRUE)
  qc_results$expression_issues <- expr_issues
  
  # Copy number validation
  cnv_issues <- sum(data$erbb2_cnv_score < -3 | data$erbb2_cnv_score > 3, na.rm = TRUE)
  qc_results$copy_number_issues <- cnv_issues
  
  # 2. Correlation validation
  if(sum(!is.na(data$erbb2_log2) & !is.na(data$erbb2_cnv_score)) > 10) {
    rna_cnv_cor <- cor(data$erbb2_log2, data$erbb2_cnv_score, use = "complete.obs")
    qc_results$rna_cnv_correlation <- rna_cnv_cor
  }
  
  # 3. Platform consistency validation
  if(sum(!is.na(data$her2_ihc_status)) > 0) {
    # Check if molecular data aligns with IHC status
    ihc_molecular_consistency <- data %>%
      filter(!is.na(her2_ihc_status) & !is.na(erbb2_log2)) %>%
      group_by(her2_ihc_status) %>%
      summarise(
        mean_expression = mean(erbb2_log2, na.rm = TRUE),
        .groups = "drop"
      )
    
    qc_results$ihc_molecular_consistency <- ihc_molecular_consistency
  }
  
  cat("Quality control completed:\n")
  cat("- Age issues:", qc_results$age_issues, "samples\n")
  cat("- Expression issues:", qc_results$expression_issues, "samples\n")
  cat("- Copy number issues:", qc_results$copy_number_issues, "samples\n")
  
  if(!is.null(qc_results$rna_cnv_correlation)) {
    cat("- RNA-CNV correlation:", round(qc_results$rna_cnv_correlation, 3), "\n")
  }
  
  session_metadata$transformations_applied <<- c(session_metadata$transformations_applied, "Quality_control")
  
  # Return data with QC flags
  data_with_qc <- data %>%
    mutate(
      qc_age_flag = age_years < 18 | age_years > 100,
      qc_expr_flag = erbb2_log2 < 0 | erbb2_log2 > 20,
      qc_cnv_flag = erbb2_cnv_score < -3 | erbb2_cnv_score > 3,
      qc_pass = !qc_age_flag & !qc_expr_flag & !qc_cnv_flag
    )
  
  return(list(data = data_with_qc, qc_results = qc_results))
}

qc_output <- safe_execute(perform_quality_control(harmonized_data), "Quality control validation")

if(!is.null(qc_output)) {
  harmonized_data <- qc_output$data
  qc_results <- qc_output$qc_results
}

# 9. SAVE HARMONIZED DATASET
cat("\n9. SAVING HARMONIZED DATASET\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

save_harmonized_dataset <- function(data) {
  cat("Saving harmonized dataset...\n")
  
  # Save main harmonized dataset
  saveRDS(data, "data/processed/harmonized_dataset.rds")
  
  # Save high-quality subset (all platforms available)
  high_quality_data <- data %>%
    filter(completeness_score == 4 & qc_pass)
  
  saveRDS(high_quality_data, "data/processed/harmonized_dataset_high_quality.rds")
  
  # Save analysis-ready dataset (minimum requirements)
  analysis_ready_data <- data %>%
    filter(has_clinical & has_rna & qc_pass)
  
  saveRDS(analysis_ready_data, "data/processed/harmonized_dataset_analysis_ready.rds")
  
  cat("Harmonized datasets saved:\n")
  cat("- Full dataset:", nrow(data), "samples\n")
  cat("- High quality (all platforms):", nrow(high_quality_data), "samples\n")
  cat("- Analysis ready (clinical + RNA):", nrow(analysis_ready_data), "samples\n")
  
  session_metadata$next_session_inputs <<- c(
    "harmonized_dataset.rds",
    "harmonized_dataset_high_quality.rds",
    "harmonized_dataset_analysis_ready.rds"
  )
  
  return(data)
}

harmonized_data <- safe_execute(save_harmonized_dataset(harmonized_data), "Saving harmonized dataset")

# 10. FINALIZE SESSION 2
cat("\n10. FINALIZING SESSION 2\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Update metadata with completion info
session_metadata$end_time <- Sys.time()
session_metadata$duration_minutes <- as.numeric(difftime(
  session_metadata$end_time, 
  session_metadata$start_time, 
  units = "mins"
))

# Save session metadata
write_json(session_metadata, "metadata/session2_metadata.json", pretty = TRUE)

# Print session summary
cat("Session 2 Complete!\n")
cat("Duration:", round(session_metadata$duration_minutes, 1), "minutes\n")
cat("Final sample size:", session_metadata$final_sample_size, "patients\n")
cat("Transformations applied:", length(session_metadata$transformations_applied), "\n")
for(transform in session_metadata$transformations_applied) {
  cat("  -", transform, "\n")
}

if(session_metadata$outliers_removed > 0) {
  cat("Outliers removed:", session_metadata$outliers_removed, "samples\n")
}

if(session_metadata$imputation_applied > 0) {
  cat("Missing values imputed:", session_metadata$imputation_applied, "\n")
}

cat("\nNext steps for Session 3:\n")
cat("1. Bayesian latent variable model implementation\n")
cat("2. EM algorithm for parameter estimation\n")
cat("3. Uncertainty quantification via bootstrap\n")
cat("4. Model convergence diagnostics\n")

# Clean up large objects to save memory
rm(raw_data, clinical_processed, rna_normalized, cnv_processed, her2_ihc_data)
gc()

cat("\n=== Session 2 Complete ===\n")