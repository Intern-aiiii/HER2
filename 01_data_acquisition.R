# HER2 Integration Project - Session 1: Complete TCGA Data Acquisition
# Author: Research Team
# Date: January 2025
# Objective: Complete HER2 biomarker harmonization data acquisition
# Target: Journal of Molecular Diagnostics

# PROJECT CONTEXT:
# - Developing Bayesian latent variable model for HER2 IHC, RNA, and CNV integration
# - Need complete TCGA-BRCA dataset with multi-platform HER2 data
# - Requires comprehensive quality assessment and metadata tracking

# Load required libraries with progress tracking
cat("=== HER2 BIOMARKER HARMONIZATION PROJECT ===\n")
cat("Session 1: Complete Data Acquisition\n")
cat("Target: Journal of Molecular Diagnostics\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(stringr)
  library(survival)
  library(tidyr)
  library(ggplot2)  # For quality plots
  library(DT)       # For HTML tables
})

# Load TCGAbiolinks
tryCatch({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  tcga_available <- TRUE
  cat("✓ TCGAbiolinks available for data downloads\n")
}, error = function(e) {
  tcga_available <- FALSE
  cat("⚠️  TCGAbiolinks not available - local files only\n")
})

# Create comprehensive project structure
create_project_structure <- function() {
  dirs <- c("data/raw/TCGA-BRCA", "data/processed", "data/results", 
            "metadata", "results/figures", "results/tables", "results/reports",
            "quality_control", "logs")
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

create_project_structure()

# Comprehensive session metadata initialization
session_metadata <- list(
  session_id = 1,
  date = Sys.Date(),
  objective = "Complete HER2 biomarker data acquisition and initial processing",
  project_context = "Bayesian latent variable model for HER2 harmonization",
  target_journal = "Journal of Molecular Diagnostics",
  start_time = Sys.time(),
  datasets_created = character(),
  sample_counts = list(
    total_cases = 0,
    her2_positive = 0,
    her2_negative = 0,
    complete_data = 0,
    clinical_only = 0,
    with_rna = 0,
    with_cnv = 0,
    all_platforms = 0
  ),
  quality_metrics = list(
    missing_data_percent = list(),
    outliers_detected = 0,
    data_completeness = list(),
    platform_concordance = list()
  ),
  data_sources = list(
    clinical = "TBD",
    survival = "TBD", 
    rna = "TBD",
    cnv = "TBD"
  ),
  next_session_inputs = character(),
  errors_encountered = character(),
  warnings_generated = character()
)

cat("Start time:", as.character(session_metadata$start_time), "\n\n")

# =============================================================================
# 1. DATA SOURCE IDENTIFICATION AND VALIDATION
# =============================================================================

cat("1. Identifying and validating data sources...\n")

# Your actual files (confirmed from previous check)
local_files <- list(
  clinical = "data/raw/TCGA-BRCA/clinical_clean.tsv",
  survival = "data/raw/TCGA-BRCA/survival_clean.tsv"
)

# Validate local files
validate_local_file <- function(filepath, file_type) {
  if(!file.exists(filepath)) {
    return(list(status = "missing", size_mb = 0, message = "File not found"))
  }
  
  tryCatch({
    # Test read first few rows
    test_data <- read_tsv(filepath, n_max = 5, show_col_types = FALSE)
    file_info <- file.info(filepath)
    
    return(list(
      status = "valid",
      size_mb = round(file_info$size / (1024^2), 2),
      rows_sample = nrow(test_data),
      cols = ncol(test_data),
      column_names = names(test_data),
      message = "File validated successfully"
    ))
  }, error = function(e) {
    return(list(
      status = "error",
      size_mb = 0,
      message = paste("Validation failed:", e$message)
    ))
  })
}

# Validate all local files
cat("Validating local data files:\n")
file_validation <- list()
for(name in names(local_files)) {
  validation <- validate_local_file(local_files[[name]], name)
  file_validation[[name]] <- validation
  
  if(validation$status == "valid") {
    cat("  ✓", name, ":", validation$size_mb, "MB,", validation$cols, "columns\n")
    session_metadata$data_sources[[name]] <- "Local file"
  } else {
    cat("  ✗", name, ":", validation$message, "\n")
  }
}

# =============================================================================
# 2. COMPREHENSIVE CLINICAL DATA ACQUISITION
# =============================================================================

cat("\n2. Comprehensive clinical data acquisition with HER2 focus...\n")

acquire_clinical_data_comprehensive <- function() {
  cat("Loading clinical data with HER2-specific processing...\n")
  
  # First try the clinical supplement file with biomarker data
  clinical_supplement_path <- "data/GDCdata/TCGA-BRCA/Clinical/Clinical_Supplement/8162d394-8b64-4da2-9f5b-d164c54b9608/nationwidechildrens.org_clinical_patient_brca.txt"
  
  if(file.exists(clinical_supplement_path)) {
    cat("✓ Found clinical supplement file with biomarker data\n")
    # Read the header row to get column names
    header_line <- readLines(clinical_supplement_path, n = 1)
    col_names <- unlist(strsplit(header_line, "\t"))
    
    # Read the data skipping the first 2 rows (header and CDE_ID row)
    clinical_raw <- read_tsv(clinical_supplement_path, skip = 2, col_names = col_names, show_col_types = FALSE)
    cat("✓ Loaded clinical supplement data:", nrow(clinical_raw), "rows,", ncol(clinical_raw), "columns\n")
    use_supplement_file <- TRUE
  } else {
    # Fall back to original clinical file
    cat("Clinical supplement file not found, using basic clinical file...\n")
    
    if(file_validation$clinical$status != "valid") {
      stop("Clinical data file validation failed - cannot proceed")
    }
    
    clinical_raw <- read_tsv(local_files$clinical, show_col_types = FALSE)
    cat("✓ Loaded clinical data:", nrow(clinical_raw), "rows,", ncol(clinical_raw), "columns\n")
    use_supplement_file <- FALSE
  }
  
  # Display all available columns for HER2 identification
  cat("Available clinical columns:\n")
  for(i in 1:length(names(clinical_raw))) {
    cat(sprintf("  %2d: %s\n", i, names(clinical_raw)[i]))
  }
  
  # Enhanced column mapping for HER2 project
  safe_coalesce <- function(data, col_options, default_value = NA) {
    for(col in col_options) {
      if(col %in% names(data)) {
        cat("    → Using column:", col, "\n")
        return(data[[col]])
      }
    }
    cat("    → No column found from:", paste(col_options, collapse = ", "), "\n")
    return(rep(default_value, nrow(data)))
  }
  
  cat("\nProcessing clinical variables with HER2-specific focus:\n")
  
  if(use_supplement_file) {
    # Process the clinical supplement file with actual biomarker columns
    clinical_processed <- clinical_raw %>%
      mutate(
        # Core identifiers
        patient_id = str_sub(bcr_patient_barcode, 1, 12),
        
        # Demographics
        age = as.numeric(age_at_diagnosis),
        gender = toupper(gender),
        race = race,
        ethnicity = ethnicity,
        
        # HER2 IHC Status (using actual column names from supplement file)
        her2_ihc_status = her2_status_by_ihc,
        her2_ihc_score = her2_ihc_score,
        her2_ihc_percent = her2_ihc_percent_positive,
        
        # HER2 FISH Status
        her2_fish_status = her2_fish_status,
        her2_copy_number = as.numeric(her2_copy_number),
        her2_cent17_ratio = as.numeric(her2_cent17_ratio),
        
        # Hormone receptor status
        er_status = er_status_by_ihc,
        er_percent = er_status_ihc_Percent_Positive,
        pr_status = pr_status_by_ihc,
        pr_percent = pr_status_ihc_percent_positive,
        
        # Tumor characteristics
        histological_type = histological_type,
        histological_grade = NA_character_,  # Not available in this file
        pathologic_stage = ajcc_pathologic_tumor_stage,
        pathologic_t = ajcc_tumor_pathologic_pt,
        pathologic_n = ajcc_nodes_pathologic_pn,
        pathologic_m = ajcc_metastasis_pathologic_pm,
        
        # Survival variables
        vital_status = vital_status,
        days_to_death = as.numeric(death_days_to),
        days_to_last_follow_up = as.numeric(last_contact_days_to),
        
        # Additional fields from basic file
        molecular_subtype = NA_character_
      )
  } else {
    # Process the basic clinical file (original logic)
    clinical_processed <- clinical_raw %>%
      mutate(
        # Core identifiers
        patient_id = safe_coalesce(., c("submitter_id", "id", "bcr_patient_barcode", "patient_id", "Patient_ID")),
        
        # Demographics
        age = safe_coalesce(., c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis", "age")),
        gender = safe_coalesce(., c("gender", "sex")),
        race = safe_coalesce(., c("race", "race_list")),
        ethnicity = safe_coalesce(., c("ethnicity", "hispanic_or_latino_ethnicity")),
        
        # HER2 IHC Status (CRITICAL for project)
        her2_ihc_status = safe_coalesce(., c("her2_ihc_status", "her2_status", "her2_ihc", "erbb2_ihc", 
                                             "her2_immunohistochemistry", "her2_protein_expression")),
        her2_ihc_score = NA_character_,
        her2_ihc_percent = NA_character_,
        
        # Additional HER2-related markers
        her2_fish_status = safe_coalesce(., c("her2_fish_status", "her2_fish", "erbb2_fish")),
        her2_copy_number = safe_coalesce(., c("her2_copy_number", "erbb2_copy_number")),
        her2_cent17_ratio = NA_real_,
        
        # Hormone receptor status
        er_status = safe_coalesce(., c("er_status", "estrogen_receptor_status", "er_ihc_status")),
        er_percent = NA_character_,
        pr_status = safe_coalesce(., c("pr_status", "progesterone_receptor_status", "pr_ihc_status")),
        pr_percent = NA_character_,
        
        # Tumor characteristics
        histological_type = safe_coalesce(., c("histological_type", "histologic_diagnosis", "histology")),
        histological_grade = safe_coalesce(., c("histological_grade", "grade", "tumor_grade", "nottingham_grade")),
        pathologic_stage = safe_coalesce(., c("pathologic_stage", "ajcc_pathologic_stage", "stage")),
        pathologic_t = safe_coalesce(., c("pathologic_t", "tumor_stage", "ajcc_pathologic_t")),
        pathologic_n = safe_coalesce(., c("pathologic_n", "node_stage", "ajcc_pathologic_n")),
        pathologic_m = safe_coalesce(., c("pathologic_m", "metastasis_stage", "ajcc_pathologic_m")),
        
        # Molecular subtype
        molecular_subtype = safe_coalesce(., c("molecular_subtype", "subtype", "breast_carcinoma_subtype")),
        
        # Survival variables
        vital_status = safe_coalesce(., c("vital_status", "overall_survival_status")),
        days_to_death = safe_coalesce(., c("days_to_death", "overall_survival")),
        days_to_last_follow_up = safe_coalesce(., c("days_to_last_follow_up", "days_to_last_followup"))
      )
  }
  
  # Common processing for both file types
  clinical_processed <- clinical_processed %>%
    
    # Data processing and standardization
    mutate(
      # Clean patient ID
      patient_id = str_sub(as.character(patient_id), 1, 12),
      
      # Standardize age
      age_years = case_when(
        !is.na(age) & age > 365 ~ age / 365.25,  # Convert days to years
        !is.na(age) ~ age,                       # Already in years
        TRUE ~ NA_real_
      ),
      
      # Standardize HER2 IHC status (CRITICAL)
      her2_ihc_standardized = case_when(
        # Use her2_ihc_score if available (from supplement file)
        !is.na(her2_ihc_score) & her2_ihc_score == "0" ~ "0",
        !is.na(her2_ihc_score) & her2_ihc_score == "1+" ~ "1+",
        !is.na(her2_ihc_score) & her2_ihc_score == "2+" ~ "2+",
        !is.na(her2_ihc_score) & her2_ihc_score == "3+" ~ "3+",
        # Otherwise use her2_ihc_status
        !is.na(her2_ihc_status) & str_detect(her2_ihc_status, regex("Negative", ignore_case = TRUE)) ~ "0",
        !is.na(her2_ihc_status) & str_detect(her2_ihc_status, regex("Equivocal|Indeterminate", ignore_case = TRUE)) ~ "2+",
        !is.na(her2_ihc_status) & str_detect(her2_ihc_status, regex("Positive", ignore_case = TRUE)) ~ "3+",
        !is.na(her2_ihc_status) & str_detect(her2_ihc_status, regex("0\\+|0", ignore_case = TRUE)) ~ "0",
        !is.na(her2_ihc_status) & str_detect(her2_ihc_status, regex("1\\+|1", ignore_case = TRUE)) ~ "1+",
        !is.na(her2_ihc_status) & str_detect(her2_ihc_status, regex("2\\+|2", ignore_case = TRUE)) ~ "2+", 
        !is.na(her2_ihc_status) & str_detect(her2_ihc_status, regex("3\\+|3", ignore_case = TRUE)) ~ "3+",
        TRUE ~ NA_character_
      ),
      
      # Binary HER2 classification
      her2_positive = case_when(
        her2_ihc_standardized %in% c("3+") ~ TRUE,
        her2_ihc_standardized %in% c("0", "1+") ~ FALSE,
        her2_ihc_standardized == "2+" & !is.na(her2_fish_status) & 
          str_detect(her2_fish_status, regex("Positive|Amplified", ignore_case = TRUE)) ~ TRUE,
        her2_ihc_standardized == "2+" & !is.na(her2_fish_status) & 
          str_detect(her2_fish_status, regex("Negative|Not Amplified", ignore_case = TRUE)) ~ FALSE,
        !is.na(her2_fish_status) & str_detect(her2_fish_status, regex("Positive|Amplified", ignore_case = TRUE)) ~ TRUE,
        !is.na(her2_fish_status) & str_detect(her2_fish_status, regex("Negative|Not Amplified", ignore_case = TRUE)) ~ FALSE,
        !is.na(her2_cent17_ratio) & her2_cent17_ratio >= 2.0 ~ TRUE,
        !is.na(her2_cent17_ratio) & her2_cent17_ratio < 2.0 ~ FALSE,
        TRUE ~ NA
      ),
      
      # Standardize hormone receptor status
      er_positive = case_when(
        !is.na(er_status) & str_detect(er_status, regex("Positive|POS", ignore_case = TRUE)) ~ TRUE,
        !is.na(er_status) & str_detect(er_status, regex("Negative|NEG", ignore_case = TRUE)) ~ FALSE,
        TRUE ~ NA
      ),
      
      pr_positive = case_when(
        !is.na(pr_status) & str_detect(pr_status, regex("Positive|POS", ignore_case = TRUE)) ~ TRUE,
        !is.na(pr_status) & str_detect(pr_status, regex("Negative|NEG", ignore_case = TRUE)) ~ FALSE,
        TRUE ~ NA
      ),
      
      # Molecular subtype classification
      molecular_subtype_standardized = case_when(
        # HER2-enriched
        her2_positive == TRUE & (er_positive == FALSE | is.na(er_positive)) & 
          (pr_positive == FALSE | is.na(pr_positive)) ~ "HER2-enriched",
        
        # Luminal HER2
        her2_positive == TRUE & (er_positive == TRUE | pr_positive == TRUE) ~ "Luminal-HER2",
        
        # Luminal A/B (HER2 negative, hormone receptor positive)
        her2_positive == FALSE & (er_positive == TRUE | pr_positive == TRUE) ~ "Luminal",
        
        # Triple negative
        her2_positive == FALSE & er_positive == FALSE & pr_positive == FALSE ~ "Triple-negative",
        
        # Use original classification if available
        !is.na(molecular_subtype) ~ molecular_subtype,
        TRUE ~ "Unknown"
      ),
      
      # Standardize histological grade
      histological_grade_standardized = case_when(
        !is.na(histological_grade) & str_detect(histological_grade, regex("1|I|Low|Well", ignore_case = TRUE)) ~ "Grade 1",
        !is.na(histological_grade) & str_detect(histological_grade, regex("2|II|Intermediate|Moderate", ignore_case = TRUE)) ~ "Grade 2", 
        !is.na(histological_grade) & str_detect(histological_grade, regex("3|III|High|Poor", ignore_case = TRUE)) ~ "Grade 3",
        !is.na(histological_grade) ~ histological_grade,
        TRUE ~ NA_character_
      ),
      
      # Standardize pathologic stage
      pathologic_stage_clean = case_when(
        !is.na(pathologic_stage) & str_detect(pathologic_stage, regex("Stage I[^IV]|Stage 1", ignore_case = TRUE)) ~ "I",
        !is.na(pathologic_stage) & str_detect(pathologic_stage, regex("Stage II|Stage 2", ignore_case = TRUE)) ~ "II",
        !is.na(pathologic_stage) & str_detect(pathologic_stage, regex("Stage III|Stage 3", ignore_case = TRUE)) ~ "III", 
        !is.na(pathologic_stage) & str_detect(pathologic_stage, regex("Stage IV|Stage 4", ignore_case = TRUE)) ~ "IV",
        !is.na(pathologic_stage) ~ pathologic_stage,
        TRUE ~ NA_character_
      ),
      
      # Create survival variables
      survival_days = case_when(
        vital_status == "Dead" | vital_status == "1" ~ as.numeric(days_to_death),
        !is.na(days_to_last_follow_up) ~ as.numeric(days_to_last_follow_up),
        TRUE ~ NA_real_
      ),
      
      survival_event = case_when(
        vital_status == "Dead" | vital_status == "1" ~ 1,
        vital_status == "Alive" | vital_status == "0" ~ 0,
        TRUE ~ NA_real_
      ),
      
      survival_years = survival_days / 365.25
    ) %>%
    
    filter(!is.na(patient_id)) %>%
    distinct(patient_id, .keep_all = TRUE)
  
  # Quality assessment
  cat("\nClinical data quality assessment:\n")
  cat(sprintf("  - Total patients: %d\n", nrow(clinical_processed)))
  cat(sprintf("  - Age available: %d (%.1f%%)\n", 
              sum(!is.na(clinical_processed$age_years)),
              100 * sum(!is.na(clinical_processed$age_years)) / nrow(clinical_processed)))
  cat(sprintf("  - HER2 IHC status available: %d (%.1f%%)\n",
              sum(!is.na(clinical_processed$her2_ihc_standardized)),
              100 * sum(!is.na(clinical_processed$her2_ihc_standardized)) / nrow(clinical_processed)))
  cat(sprintf("  - HER2 positive cases: %d (%.1f%%)\n",
              sum(clinical_processed$her2_positive == TRUE, na.rm = TRUE),
              100 * sum(clinical_processed$her2_positive == TRUE, na.rm = TRUE) / nrow(clinical_processed)))
  cat(sprintf("  - Survival data available: %d (%.1f%%)\n",
              sum(!is.na(clinical_processed$survival_days)),
              100 * sum(!is.na(clinical_processed$survival_days)) / nrow(clinical_processed)))
  
  return(clinical_processed)
}

clinical_data <- acquire_clinical_data_comprehensive()

# =============================================================================
# 3. SURVIVAL DATA INTEGRATION
# =============================================================================

cat("\n3. Integrating survival data with enhanced processing...\n")

integrate_survival_data <- function(survival_file_path, clinical_data) {
  if(file_validation$survival$status != "valid") {
    cat("Survival file not valid - using clinical survival data only\n")
    return(clinical_data)
  }
  
  survival_raw <- read_tsv(survival_file_path, show_col_types = FALSE)
  cat("✓ Loaded survival file:", nrow(survival_raw), "rows\n")
  cat("Survival columns:", paste(names(survival_raw), collapse = ", "), "\n")
  
  # Process survival data (based on your file structure: sample, OS.time, OS, _PATIENT)
  survival_processed <- survival_raw %>%
    mutate(
      patient_id = case_when(
        "sample" %in% names(.) ~ str_sub(as.character(sample), 1, 12),
        "_PATIENT" %in% names(.) ~ str_sub(as.character(`_PATIENT`), 1, 12),
        TRUE ~ NA_character_
      ),
      
      # Overall survival time
      os_time_external = case_when(
        "OS.time" %in% names(.) & !is.na(OS.time) & OS.time > 10 ~ OS.time / 365.25,
        "OS.time" %in% names(.) & !is.na(OS.time) ~ OS.time,
        TRUE ~ NA_real_
      ),
      
      # Overall survival event
      os_event_external = case_when(
        "OS" %in% names(.) & OS == 1 ~ 1,
        "OS" %in% names(.) & OS == 0 ~ 0,
        TRUE ~ NA_real_
      ),
      
      survival_days_external = os_time_external * 365.25
    ) %>%
    filter(!is.na(patient_id)) %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    select(patient_id, survival_days_external, os_event_external)
  
  # Merge with clinical data
  final_data <- clinical_data %>%
    left_join(survival_processed, by = "patient_id") %>%
    mutate(
      # Use external survival data preferentially
      survival_days_final = coalesce(survival_days_external, survival_days),
      survival_event_final = coalesce(os_event_external, survival_event),
      
      # Update main survival columns
      survival_days = survival_days_final,
      survival_event = survival_event_final,
      survival_years = survival_days / 365.25
    ) %>%
    select(-survival_days_external, -os_event_external, -survival_days_final, -survival_event_final)
  
  cat("Survival integration complete\n")
  session_metadata$data_sources$survival <- "Local file integrated"
  
  return(final_data)
}

clinical_final <- integrate_survival_data(local_files$survival, clinical_data) %>%
  # Final filtering for valid cases
  filter(
    !is.na(survival_days) & !is.na(survival_event),
    survival_days > 0 & survival_days < 365.25 * 20  # Reasonable follow-up times
  )

# Update sample counts
session_metadata$sample_counts$total_cases <- nrow(clinical_final)
session_metadata$sample_counts$her2_positive <- sum(clinical_final$her2_positive == TRUE, na.rm = TRUE)
session_metadata$sample_counts$her2_negative <- sum(clinical_final$her2_positive == FALSE, na.rm = TRUE)
session_metadata$sample_counts$clinical_only <- nrow(clinical_final)

# =============================================================================
# 4. RNA EXPRESSION DATA ACQUISITION (ERBB2 FOCUS)
# =============================================================================

cat("\n4. Acquiring RNA expression data with ERBB2 focus...\n")

acquire_rna_expression <- function() {
  if(!tcga_available) {
    cat("⚠️  TCGAbiolinks not available - skipping RNA expression\n")
    session_metadata$data_sources$rna <- "Not available"
    return(NULL)
  }
  
  cat("Downloading TCGA-BRCA RNA-seq data with HER2 gene focus...\n")
  
  tryCatch({
    # Query for RNA-seq data
    rna_query <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling", 
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal")
    )
    
    cat("Found", nrow(getResults(rna_query)), "RNA-seq files\n")
    cat("Downloading RNA expression data (this may take several minutes)...\n")
    
    GDCdownload(rna_query, method = "api", files.per.chunk = 15)
    rna_data <- GDCprepare(rna_query)
    
    cat("✓ RNA expression data downloaded successfully\n")
    session_metadata$data_sources$rna <- "Downloaded from TCGA"
    
    return(rna_data)
    
  }, error = function(e) {
    cat("❌ RNA download failed:", e$message, "\n")
    session_metadata$errors_encountered <- c(session_metadata$errors_encountered, 
                                             paste("RNA download:", e$message))
    session_metadata$data_sources$rna <- "Download failed"
    return(NULL)
  })
}

rna_raw <- acquire_rna_expression()

# Process RNA data with HER2 focus
process_rna_her2_focus <- function(rna_data) {
  if(is.null(rna_data)) {
    cat("⚠️  No RNA data to process\n")
    return(NULL)
  }
  
  cat("Processing RNA data with HER2 pathway focus...\n")
  
  # Extract components
  gene_info <- rowData(rna_data)
  sample_info <- colData(rna_data)
  expression_data <- assay(rna_data)
  
  # HER2 pathway and related genes
  target_genes <- c("ERBB2", "ESR1", "PGR", "MKI67", "PCNA", "GREB1", "CCND1", "MYC")
  gene_indices <- list()
  
  cat("Finding HER2-related genes:\n")
  for (gene in target_genes) {
    idx <- which(gene_info$gene_name == gene)
    if (length(idx) == 0) {
      idx <- grep(paste0("^", gene, "$"), gene_info$gene_name, ignore.case = TRUE)
    }
    if (length(idx) > 0) {
      gene_indices[[gene]] <- idx[1]
      cat("  ✓ Found", gene, "at index", idx[1], "\n")
    } else {
      cat("  ✗", gene, "not found\n")
    }
  }
  
  if(length(gene_indices) == 0) {
    cat("❌ No target genes found\n")
    return(NULL)
  }
  
  # Create expression dataset
  rna_processed <- data.frame(
    patient_id = substr(sample_info$submitter_id, 1, 12),
    sample_id = sample_info$submitter_id,
    sample_type = sample_info$sample_type,
    sample_submitter_id = sample_info$sample_submitter_id,
    stringsAsFactors = FALSE
  )
  
  # Add expression values for HER2 pathway genes
  for(gene in names(gene_indices)) {
    idx <- gene_indices[[gene]]
    raw_counts <- as.numeric(expression_data[idx, ])
    
    rna_processed[[paste0(gene, "_counts")]] <- raw_counts
    rna_processed[[paste0(gene, "_log2")]] <- log2(raw_counts + 1)
    rna_processed[[paste0(gene, "_log2_centered")]] <- scale(log2(raw_counts + 1))[,1]
    
    # Calculate percentiles for clinical interpretation
    rna_processed[[paste0(gene, "_percentile")]] <- rank(raw_counts) / length(raw_counts) * 100
  }
  
  # Filter for primary tumors and quality control
  rna_final <- rna_processed %>%
    filter(sample_type == "Primary Tumor") %>%
    filter(!is.na(ERBB2_counts) & ERBB2_counts > 0) %>%
    # Remove outliers (optional)
    filter(ERBB2_log2 > quantile(ERBB2_log2, 0.01, na.rm = TRUE) & 
             ERBB2_log2 < quantile(ERBB2_log2, 0.99, na.rm = TRUE)) %>%
    arrange(patient_id)
  
  cat("✓ RNA processing complete:", nrow(rna_final), "primary tumor samples\n")
  cat("  Available genes:", paste(names(gene_indices), collapse = ", "), "\n")
  
  return(rna_final)
}

rna_final <- process_rna_her2_focus(rna_raw)

# Update sample counts
if(!is.null(rna_final)) {
  session_metadata$sample_counts$with_rna <- nrow(rna_final)
}

# =============================================================================
# 5. COPY NUMBER DATA ACQUISITION (ERBB2 LOCUS FOCUS)
# =============================================================================

cat("\n5. Acquiring copy number data with ERBB2 locus focus...\n")

acquire_cnv_erbb2_focus <- function() {
  if(!tcga_available) {
    cat("⚠️  TCGAbiolinks not available - skipping copy number data\n")
    session_metadata$data_sources$cnv <- "Not available"
    return(NULL)
  }
  
  cat("Downloading TCGA-BRCA copy number data for ERBB2 analysis...\n")
  
  tryCatch({
    # Try different CNV data types
    cnv_types <- list(
      list(data.type = "Masked Copy Number Segment", desc = "Segment-level"),
      list(data.type = "Gene Level Copy Number", desc = "Gene-level")
    )
    
    for(cnv_type in cnv_types) {
      tryCatch({
        cnv_query <- GDCquery(
          project = "TCGA-BRCA",
          data.category = "Copy Number Variation",
          data.type = cnv_type$data.type,
          sample.type = "Primary Tumor"
        )
        
        if(nrow(getResults(cnv_query)) > 0) {
          cat("Found", nrow(getResults(cnv_query)), cnv_type$desc, "CNV files\n")
          cat("Downloading", cnv_type$desc, "copy number data...\n")
          
          GDCdownload(cnv_query, method = "api", files.per.chunk = 20)
          cnv_data <- GDCprepare(cnv_query)
          
          cat("✓", cnv_type$desc, "copy number data downloaded\n")
          session_metadata$data_sources$cnv <- paste("Downloaded:", cnv_type$desc)
          
          return(cnv_data)
        }
      }, error = function(e) {
        cat("Failed with", cnv_type$desc, ":", e$message, "\n")
      })
    }
    
    cat("❌ All CNV download attempts failed\n")
    session_metadata$data_sources$cnv <- "Download failed"
    return(NULL)
    
  }, error = function(e) {
    cat("❌ CNV download error:", e$message, "\n")
    session_metadata$errors_encountered <- c(session_metadata$errors_encountered,
                                             paste("CNV download:", e$message))
    return(NULL)
  })
}

cnv_raw <- acquire_cnv_erbb2_focus()

# Process CNV data focusing on ERBB2 locus
process_cnv_erbb2_locus <- function(cnv_data) {
  if(is.null(cnv_data)) {
    cat("⚠️  No CNV data to process\n")
    return(NULL)
  }
  
  cat("Processing copy number data for ERBB2 locus (chr17:37,844,167-37,886,679)...\n")
  
  if(inherits(cnv_data, c("RangedSummarizedExperiment", "SummarizedExperiment"))) {
    # Gene-level data
    cat("Processing gene-level copy number data...\n")
    
    gene_info <- rowData(cnv_data)
    sample_info <- colData(cnv_data)
    
    # Find ERBB2 gene
    erbb2_idx <- which(gene_info$gene_name == "ERBB2" | gene_info$Hugo_Symbol == "ERBB2")
    if (length(erbb2_idx) == 0) {
      erbb2_idx <- grep("ERBB2", gene_info$gene_name, ignore.case = TRUE)
    }
    
    if (length(erbb2_idx) > 0) {
      erbb2_cnv <- assay(cnv_data)[erbb2_idx[1], ]
      
      cnv_processed <- data.frame(
        patient_id = substr(sample_info$submitter_id, 1, 12),
        sample_id = sample_info$submitter_id,
        erbb2_cnv_score = as.numeric(erbb2_cnv),
        sample_type = sample_info$sample_type,
        stringsAsFactors = FALSE
      )
    } else {
      cat("❌ ERBB2 not found in gene-level data\n")
      return(NULL)
    }
    
  } else if(is.data.frame(cnv_data)) {
    # Segment-level data
    cat("Processing segment-level copy number data...\n")
    
    # First, let's check what columns we have
    cat("Available columns:", paste(names(cnv_data), collapse=", "), "\n")
    
    # Create a copy of the data to work with
    cnv_work <- cnv_data
    
    # Map column names flexibly
    available_cols <- names(cnv_work)
    
    # Find and rename Sample column
    sample_cols <- available_cols[grepl("Sample|sample", available_cols, ignore.case = TRUE)]
    if(length(sample_cols) > 0) {
      names(cnv_work)[names(cnv_work) == sample_cols[1]] <- "Sample"
    }
    
    # Find and rename Chromosome column
    chr_cols <- available_cols[grepl("Chromosome|Chr|chrom", available_cols, ignore.case = TRUE)]
    if(length(chr_cols) > 0) {
      names(cnv_work)[names(cnv_work) == chr_cols[1]] <- "Chromosome"
    }
    
    # Find and rename Start column
    start_cols <- available_cols[grepl("Start|start", available_cols, ignore.case = TRUE)]
    if(length(start_cols) > 0) {
      names(cnv_work)[names(cnv_work) == start_cols[1]] <- "Start"
    }
    
    # Find and rename End column
    end_cols <- available_cols[grepl("End|end", available_cols, ignore.case = TRUE)]
    if(length(end_cols) > 0) {
      names(cnv_work)[names(cnv_work) == end_cols[1]] <- "End"
    }
    
    # Find and rename Segment_Mean column
    seg_cols <- available_cols[grepl("Segment_Mean|Segment.Mean|Copy_Number|Log2", available_cols, ignore.case = TRUE)]
    if(length(seg_cols) > 0) {
      names(cnv_work)[names(cnv_work) == seg_cols[1]] <- "Segment_Mean"
    }
    
    # Check if we have all required columns
    required_cols <- c("Sample", "Chromosome", "Start", "End", "Segment_Mean")
    missing_cols <- setdiff(required_cols, names(cnv_work))
    
    if(length(missing_cols) > 0) {
      cat("❌ Missing required columns:", paste(missing_cols, collapse=", "), "\n")
      cat("Available columns after mapping:", paste(names(cnv_work), collapse=", "), "\n")
      return(NULL)
    }
    
    # Process ERBB2 locus segments
    cnv_processed <- cnv_work %>%
      filter(
        Chromosome == "17" | Chromosome == "chr17",
        Start <= 37886679,
        End >= 37844167
      ) %>%
      group_by(Sample) %>%
      summarise(
        erbb2_cnv_score = mean(Segment_Mean, na.rm = TRUE),
        erbb2_cnv_max = max(Segment_Mean, na.rm = TRUE),
        erbb2_cnv_min = min(Segment_Mean, na.rm = TRUE),
        n_segments = n(),
        .groups = "drop"
      ) %>%
      mutate(
        patient_id = substr(Sample, 1, 12),
        sample_id = Sample,
        sample_type = "Primary Tumor"
      )
  } else {
    cat("❌ Unrecognized CNV data format\n")
    return(NULL)
  }
  
  # Add clinical interpretation
  cnv_final <- cnv_processed %>%
    filter(!is.na(erbb2_cnv_score)) %>%
    mutate(
      # Standard thresholds for copy number interpretation
      erbb2_amplified = erbb2_cnv_score > 0.3,        # Log2 ratio > 0.3 suggests amplification
      erbb2_deleted = erbb2_cnv_score < -0.3,         # Log2 ratio < -0.3 suggests deletion
      erbb2_neutral = abs(erbb2_cnv_score) <= 0.3,    # Within normal range
      
      # Clinical categories
      erbb2_cnv_category = case_when(
        erbb2_cnv_score > 0.7 ~ "High amplification",
        erbb2_cnv_score > 0.3 ~ "Amplification", 
        erbb2_cnv_score > 0.1 ~ "Gain",
        erbb2_cnv_score > -0.1 ~ "Neutral",
        erbb2_cnv_score > -0.3 ~ "Loss",
        TRUE ~ "Deletion"
      ),
      
      # Percentile ranking
      erbb2_cnv_percentile = rank(erbb2_cnv_score) / length(erbb2_cnv_score) * 100
    ) %>%
    distinct(patient_id, .keep_all = TRUE)
  
  cat("✓ CNV processing complete:", nrow(cnv_final), "samples\n")
  cat("  - Amplified samples:", sum(cnv_final$erbb2_amplified), 
      sprintf("(%.1f%%)\n", 100 * sum(cnv_final$erbb2_amplified) / nrow(cnv_final)))
  cat("  - Deleted samples:", sum(cnv_final$erbb2_deleted),
      sprintf("(%.1f%%)\n", 100 * sum(cnv_final$erbb2_deleted) / nrow(cnv_final)))
  
  return(cnv_final)
}

cnv_final <- process_cnv_erbb2_locus(cnv_raw)

# Update sample counts
if(!is.null(cnv_final)) {
  session_metadata$sample_counts$with_cnv <- nrow(cnv_final)
}

# =============================================================================
# 6. COMPREHENSIVE QUALITY CONTROL ASSESSMENT
# =============================================================================

cat("\n6. Comprehensive quality control assessment...\n")

perform_quality_assessment <- function(clinical_data, rna_data, cnv_data) {
  cat("Performing comprehensive data quality assessment...\n")
  
  # Clinical data quality metrics
  clinical_quality <- list(
    total_patients = nrow(clinical_data),
    age_missing = sum(is.na(clinical_data$age_years)),
    age_missing_pct = 100 * sum(is.na(clinical_data$age_years)) / nrow(clinical_data),
    her2_ihc_missing = sum(is.na(clinical_data$her2_ihc_standardized)),
    her2_ihc_missing_pct = 100 * sum(is.na(clinical_data$her2_ihc_standardized)) / nrow(clinical_data),
    survival_missing = sum(is.na(clinical_data$survival_days)),
    survival_missing_pct = 100 * sum(is.na(clinical_data$survival_days)) / nrow(clinical_data)
  )
  
  # Multi-platform completeness
  common_patients <- clinical_data$patient_id
  if(!is.null(rna_data)) {
    common_patients <- intersect(common_patients, rna_data$patient_id)
  }
  if(!is.null(cnv_data)) {
    common_patients <- intersect(common_patients, cnv_data$patient_id)
  }
  
  multi_platform_quality <- list(
    clinical_patients = nrow(clinical_data),
    rna_patients = ifelse(is.null(rna_data), 0, nrow(rna_data)),
    cnv_patients = ifelse(is.null(cnv_data), 0, nrow(cnv_data)),
    complete_data_patients = length(common_patients),
    complete_data_pct = 100 * length(common_patients) / nrow(clinical_data)
  )
  
  # HER2 status distribution
  her2_distribution <- tryCatch({
    clinical_data %>%
      filter(!is.na(her2_ihc_standardized)) %>%
      group_by(her2_ihc_standardized) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(percentage = 100 * count / sum(count))
  }, error = function(e) {
    cat("Error in HER2 distribution calculation:", e$message, "\n")
    data.frame(her2_ihc_standardized = NA, count = 0, percentage = 0)
  })
  
  # Outlier detection (basic)
  outliers_detected <- 0
  if(!is.null(rna_data)) {
    outliers_detected <- outliers_detected + sum(rna_data$ERBB2_log2 > quantile(rna_data$ERBB2_log2, 0.99, na.rm = TRUE) | 
                                                   rna_data$ERBB2_log2 < quantile(rna_data$ERBB2_log2, 0.01, na.rm = TRUE), na.rm = TRUE)
  }
  
  # Store quality metrics
  session_metadata$quality_metrics$missing_data_percent <- list(
    clinical_age = clinical_quality$age_missing_pct,
    clinical_her2 = clinical_quality$her2_ihc_missing_pct,
    clinical_survival = clinical_quality$survival_missing_pct
  )
  
  session_metadata$quality_metrics$outliers_detected <- outliers_detected
  session_metadata$quality_metrics$data_completeness <- multi_platform_quality
  
  session_metadata$sample_counts$complete_data <- multi_platform_quality$complete_data_patients
  session_metadata$sample_counts$all_platforms <- length(common_patients)
  
  cat("Quality assessment complete:\n")
  cat(sprintf("  - Clinical data completeness: %.1f%% age, %.1f%% HER2, %.1f%% survival\n",
              100 - clinical_quality$age_missing_pct,
              100 - clinical_quality$her2_ihc_missing_pct,
              100 - clinical_quality$survival_missing_pct))
  cat(sprintf("  - Multi-platform complete data: %d patients (%.1f%%)\n",
              multi_platform_quality$complete_data_patients,
              multi_platform_quality$complete_data_pct))
  
  return(list(
    clinical = clinical_quality,
    multiplatform = multi_platform_quality,
    her2_distribution = her2_distribution,
    outliers = outliers_detected
  ))
}

quality_metrics <- perform_quality_assessment(clinical_final, rna_final, cnv_final)

# =============================================================================
# 7. SAVE PROCESSED DATASETS
# =============================================================================

cat("\n7. Saving processed datasets...\n")

# Save clinical data
saveRDS(clinical_final, "data/processed/tcga_clinical.rds")
session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_clinical.rds")
cat("✓ Saved tcga_clinical.rds:", nrow(clinical_final), "patients\n")

# Save RNA data if available
if(!is.null(rna_final)) {
  saveRDS(rna_final, "data/processed/tcga_expression.rds")
  session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_expression.rds")
  cat("✓ Saved tcga_expression.rds:", nrow(rna_final), "samples\n")
}

# Save CNV data if available
if(!is.null(cnv_final)) {
  saveRDS(cnv_final, "data/processed/tcga_cnv.rds")
  session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_cnv.rds")
  cat("✓ Saved tcga_cnv.rds:", nrow(cnv_final), "samples\n")
}

# =============================================================================
# 8. GENERATE COMPREHENSIVE DATA QUALITY REPORT (HTML)
# =============================================================================

cat("\n8. Generating comprehensive data quality report...\n")

# Calculate duration before using it
session_metadata$end_time <- Sys.time()
session_metadata$duration_minutes <- as.numeric(difftime(
  session_metadata$end_time, 
  session_metadata$start_time, 
  units = "mins"
))

generate_quality_report_html <- function(session_metadata, quality_metrics) {
  
  # Create HTML content
  html_content <- paste0(
    "<!DOCTYPE html>",
    "<html><head>",
    "<title>TCGA-BRCA HER2 Biomarker Data Quality Report</title>",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 40px; }",
    "h1 { color: #2c3e50; border-bottom: 3px solid #3498db; }",
    "h2 { color: #34495e; border-bottom: 1px solid #bdc3c7; }",
    "table { border-collapse: collapse; width: 100%; margin: 20px 0; }",
    "th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }",
    "th { background-color: #f2f2f2; font-weight: bold; }",
    ".success { color: #27ae60; font-weight: bold; }",
    ".warning { color: #f39c12; font-weight: bold; }",
    ".error { color: #e74c3c; font-weight: bold; }",
    ".metric { background-color: #ecf0f1; padding: 10px; margin: 10px 0; border-left: 4px solid #3498db; }",
    "</style>",
    "</head><body>",
    
    "<h1>TCGA-BRCA HER2 Biomarker Harmonization Project</h1>",
    "<h2>Session 1: Data Quality Report</h2>",
    
    "<div class='metric'>",
    "<strong>Report Generated:</strong> ", Sys.time(), "<br>",
    "<strong>Session Duration:</strong> ", round(session_metadata$duration_minutes, 1), " minutes<br>",
    "<strong>Project Objective:</strong> ", session_metadata$objective, "<br>",
    "<strong>Target Journal:</strong> ", session_metadata$target_journal,
    "</div>",
    
    "<h2>Data Source Summary</h2>",
    "<table>",
    "<tr><th>Data Type</th><th>Source</th><th>Status</th><th>Sample Count</th></tr>",
    "<tr><td>Clinical/Survival</td><td>", session_metadata$data_sources$clinical, "</td><td class='success'>✓ Available</td><td>", session_metadata$sample_counts$total_cases, "</td></tr>",
    "<tr><td>RNA Expression</td><td>", session_metadata$data_sources$rna, "</td><td class='", 
    ifelse(session_metadata$sample_counts$with_rna > 0, "success'>✓ Available", "warning'>⚠ Missing"), "</td><td>", session_metadata$sample_counts$with_rna, "</td></tr>",
    "<tr><td>Copy Number</td><td>", session_metadata$data_sources$cnv, "</td><td class='",
    ifelse(session_metadata$sample_counts$with_cnv > 0, "success'>✓ Available", "warning'>⚠ Missing"), "</td><td>", session_metadata$sample_counts$with_cnv, "</td></tr>",
    "</table>",
    
    "<h2>HER2 Biomarker Analysis Readiness</h2>",
    "<div class='metric'>",
    "<strong>Total Cases:</strong> ", session_metadata$sample_counts$total_cases, "<br>",
    "<strong>HER2 Positive:</strong> ", session_metadata$sample_counts$her2_positive, " (", 
    round(100 * session_metadata$sample_counts$her2_positive / session_metadata$sample_counts$total_cases, 1), "%)<br>",
    "<strong>Complete Multi-Platform Data:</strong> ", session_metadata$sample_counts$complete_data, " (", 
    round(100 * session_metadata$sample_counts$complete_data / session_metadata$sample_counts$total_cases, 1), "%)<br>",
    "</div>",
    
    "<h2>Data Quality Metrics</h2>",
    "<table>",
    "<tr><th>Variable</th><th>Missing Data %</th><th>Quality Status</th></tr>",
    "<tr><td>Age at Diagnosis</td><td>", round(session_metadata$quality_metrics$missing_data_percent$clinical_age, 1), "%</td><td class='",
    ifelse(session_metadata$quality_metrics$missing_data_percent$clinical_age < 10, "success'>Excellent", 
           ifelse(session_metadata$quality_metrics$missing_data_percent$clinical_age < 25, "warning'>Acceptable", "error'>Poor")), "</td></tr>",
    "<tr><td>HER2 IHC Status</td><td>", round(session_metadata$quality_metrics$missing_data_percent$clinical_her2, 1), "%</td><td class='",
    ifelse(session_metadata$quality_metrics$missing_data_percent$clinical_her2 < 10, "success'>Excellent",
           ifelse(session_metadata$quality_metrics$missing_data_percent$clinical_her2 < 25, "warning'>Acceptable", "error'>Poor")), "</td></tr>",
    "<tr><td>Survival Data</td><td>", round(session_metadata$quality_metrics$missing_data_percent$clinical_survival, 1), "%</td><td class='",
    ifelse(session_metadata$quality_metrics$missing_data_percent$clinical_survival < 10, "success'>Excellent",
           ifelse(session_metadata$quality_metrics$missing_data_percent$clinical_survival < 25, "warning'>Acceptable", "error'>Poor")), "</td></tr>",
    "</table>",
    
    "<h2>Analysis Capabilities</h2>",
    if(session_metadata$sample_counts$complete_data > 100) {
      "<div class='success'>✓ COMPLETE BAYESIAN INTEGRATION ANALYSIS READY</div>
       <ul>
       <li>✓ HER2 IHC, RNA expression, and copy number integration</li>
       <li>✓ Latent variable modeling with clinical covariates</li>
       <li>✓ Platform concordance assessment</li>
       <li>✓ Clinical outcome validation</li>
       </ul>"
    } else {
      "<div class='warning'>⚠ PARTIAL ANALYSIS CAPABILITIES</div>
       <p>Missing data platforms limit full Bayesian integration. Consider:</p>
       <ul>
       <li>Clinical associations and survival analysis (available)</li>
       <li>Single-platform molecular analysis (where available)</li>
       <li>Two-platform integration (if 2+ platforms available)</li>
       </ul>"
    },
    
    "<h2>Files Created</h2>",
    "<ul>",
    paste0("<li>", session_metadata$datasets_created, "</li>", collapse = ""),
    "</ul>",
    
    "<h2>Next Steps</h2>",
    if(session_metadata$sample_counts$complete_data > 100) {
      "<div class='success'>Ready for Session 2: Data Harmonization</div>
       <p>Proceed with complete multi-platform Bayesian modeling pipeline.</p>"
    } else {
      "<div class='warning'>Consider Data Enhancement</div>
       <p>Options: (1) Proceed with available platforms, (2) Acquire missing data types, (3) Focus on available data strengths.</p>"
    },
    
    "<hr>",
    "<p><em>Generated by HER2 Biomarker Harmonization Project - Session 1</em></p>",
    "</body></html>"
  )
  
  # Save HTML report
  writeLines(html_content, "results/reports/data_quality_report.html")
  cat("✓ HTML quality report saved: results/reports/data_quality_report.html\n")
  
  return("results/reports/data_quality_report.html")
}

html_report_path <- generate_quality_report_html(session_metadata, quality_metrics)

# =============================================================================
# 9. FINALIZE SESSION METADATA
# =============================================================================

cat("\n9. Finalizing session metadata...\n")

# Session metadata already updated in step 8

session_metadata$next_session_inputs <- session_metadata$datasets_created

# Save comprehensive metadata
write_json(session_metadata, "metadata/session1_metadata.json", pretty = TRUE)
cat("✓ Session metadata saved: metadata/session1_metadata.json\n")

# =============================================================================
# 10. FINAL SESSION SUMMARY
# =============================================================================

cat("\n=== SESSION 1 COMPLETE: HER2 BIOMARKER DATA ACQUISITION ===\n")
cat("Target: Journal of Molecular Diagnostics\n")
cat("Duration:", round(session_metadata$duration_minutes, 1), "minutes\n\n")

cat("📊 FINAL DATA SUMMARY:\n")
cat("Clinical Data:", session_metadata$sample_counts$total_cases, "patients\n")
cat("  - HER2 positive:", session_metadata$sample_counts$her2_positive, "\n")
cat("  - Deaths observed:", sum(clinical_final$survival_event == 1), "\n")
cat("  - Median follow-up:", round(median(clinical_final$survival_years), 1), "years\n")

if(session_metadata$sample_counts$with_rna > 0) {
  cat("RNA Expression:", session_metadata$sample_counts$with_rna, "samples ✓\n")
} else {
  cat("RNA Expression: Not available ❌\n")
}

if(session_metadata$sample_counts$with_cnv > 0) {
  cat("Copy Number:", session_metadata$sample_counts$with_cnv, "samples ✓\n")
} else {
  cat("Copy Number: Not available ❌\n")
}

cat("Complete Multi-Platform Data:", session_metadata$sample_counts$complete_data, "patients\n")

# Analysis readiness assessment
analysis_ready <- session_metadata$sample_counts$complete_data > 100
cat("\n🎯 ANALYSIS READINESS:\n")
if(analysis_ready) {
  cat("🎉 COMPLETE BAYESIAN INTEGRATION READY!\n")
  cat("✅ Proceed to Session 2: Data Harmonization\n")
  cat("✅ Full HER2 biomarker harmonization analysis possible\n")
} else {
  platforms_available <- sum(c(TRUE, session_metadata$sample_counts$with_rna > 0, session_metadata$sample_counts$with_cnv > 0))
  cat("⚠️  PARTIAL ANALYSIS READY (", platforms_available, " of 3 platforms)\n")
  cat("✅ Clinical associations and survival analysis ready\n")
  cat("❓ Consider acquiring missing molecular data for complete analysis\n")
}

cat("\n📁 OUTPUT FILES:\n")
for(file in session_metadata$datasets_created) {
  cat("  ✓", file, "\n")
}
cat("  ✓ session1_metadata.json\n")
cat("  ✓ data_quality_report.html\n")

cat("\n🔄 NEXT SESSION PREPARATION:\n")
cat("Input files for Session 2:", paste(session_metadata$next_session_inputs, collapse = ", "), "\n")

if(length(session_metadata$errors_encountered) > 0) {
  cat("\n⚠️  ERRORS ENCOUNTERED:\n")
  for(error in session_metadata$errors_encountered) {
    cat("  -", error, "\n")
  }
}

cat("\n✅ Session 1 complete - HER2 biomarker data ready for harmonization!\n")
cat("📖 View detailed report: results/reports/data_quality_report.html\n")
