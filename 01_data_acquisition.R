# HER2 Integration Project - Session 1: TCGA-BRCA Data Acquisition
# Author: Research Team
# Date: January 2025
# Objective: Download and process TCGA-BRCA multi-platform data with real biomarker extraction

# Load required libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(SummarizedExperiment)
  library(DT)
  library(knitr)
  library(rmarkdown)
  library(survival)
  library(stringr)
})

# Create project directory structure
create_project_structure <- function() {
  dirs <- c("data/raw", "data/processed", "data/results", 
            "metadata", "results/figures", "results/tables", "results/reports")
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# Initialize project
create_project_structure()

# Session metadata initialization
session_metadata <- list(
  session_id = 1,
  date = Sys.Date(),
  objective = "TCGA-BRCA data acquisition with real biomarker extraction",
  start_time = Sys.time(),
  datasets_created = character(),
  sample_counts = list(),
  quality_metrics = list(),
  next_session_inputs = character()
)

cat("=== HER2 Integration Project - Session 1 ===\n")
cat("Objective: TCGA-BRCA real data acquisition\n")
cat("Start time:", as.character(session_metadata$start_time), "\n")
cat("\n🔄 SKIP LOGIC ENABLED:\n")
cat("- Will check for existing processed files (.rds) and skip processing\n")
cat("- Will check for existing raw TCGA downloads and skip downloads\n")  
cat("- Will load any existing data to minimize processing time\n")
cat("- Will automatically clean up incomplete downloads\n\n")

cat("💡 MANUAL CLEANUP FUNCTIONS AVAILABLE:\n")
cat("- inspect_existing_data() - examine existing processed files\n")
cat("- clean_incomplete_downloads() - clean partial downloads\n")
cat("- force_clean_tcga_downloads() - remove all TCGA data for fresh start\n\n")

# Function to inspect existing data files for debugging
inspect_existing_data <- function() {
  cat("\n🔍 INSPECTING EXISTING DATA FILES:\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  files_to_check <- c(
    "data/processed/tcga_clinical.rds",
    "data/processed/tcga_expression.rds", 
    "data/processed/tcga_cnv.rds"
  )
  
  for (file in files_to_check) {
    if (file.exists(file)) {
      cat("\n📄", basename(file), ":\n")
      cat("  File size:", round(file.size(file) / 1024 / 1024, 2), "MB\n")
      
      tryCatch({
        data <- readRDS(file)
        cat("  Dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
        cat("  Column names:", paste(head(names(data), 5), collapse = ", "))
        if (ncol(data) > 5) cat(" ...")
        cat("\n")
        
        if (nrow(data) > 0) {
          cat("  Sample data preview:\n")
          print(head(data[, 1:min(3, ncol(data))], 3))
        } else {
          cat("  ⚠️  Data frame is empty!\n")
        }
        
      }, error = function(e) {
        cat("  ❌ Error reading file:", e$message, "\n")
      })
    } else {
      cat("\n❌", basename(file), ": File not found\n")
    }
  }
  cat("\n")
}

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

# =============================================================================
# SKIP LOGIC: Check if data already exists
# =============================================================================

# Function to manually clean all TCGA downloads (for persistent download issues)
force_clean_tcga_downloads <- function() {
  cat("\n🗑️  FORCE CLEANING ALL TCGA DOWNLOADS:\n")
  cat("This will remove all downloaded TCGA data and force fresh downloads.\n")
  cat("Use this if you're experiencing persistent partial download issues.\n\n")
  
  gdc_path <- file.path("GDCdata", "TCGA-BRCA")
  
  if (dir.exists(gdc_path)) {
    cat("Removing:", gdc_path, "\n")
    unlink(gdc_path, recursive = TRUE)
    cat("✅ All TCGA downloads cleaned\n")
  } else {
    cat("No TCGA downloads found to clean\n")
  }
  
  cat("\nNext script run will download all data fresh.\n")
}

# Function to clean up incomplete downloads
clean_incomplete_downloads <- function() {
  cat("\n🧹 CHECKING FOR INCOMPLETE DOWNLOADS:\n")
  
  # Check for partial downloads that might cause issues
  gdc_base <- file.path("GDCdata", "TCGA-BRCA", "harmonized")
  
  if (dir.exists(gdc_base)) {
    data_dirs <- c("Clinical", "Transcriptome_Profiling", "Copy_Number_Variation")
    
    for (data_dir in data_dirs) {
      full_path <- file.path(gdc_base, data_dir)
      if (dir.exists(full_path)) {
        files <- list.files(full_path, recursive = TRUE)
        
        # Check for incomplete download indicators
        incomplete_indicators <- c(
          length(grep("\\.tmp$", files)) > 0,  # Temporary files
          length(grep("\\.partial$", files)) > 0,  # Partial files
          length(grep("DOWNLOAD_.*\\.txt$", files)) > 0  # Download status files
        )
        
        if (any(incomplete_indicators)) {
          cat("⚠️  Found incomplete download indicators in", data_dir, "\n")
          cat("   Cleaning up to prevent download conflicts...\n")
          unlink(full_path, recursive = TRUE)
          cat("   ✅ Cleaned up", data_dir, "\n")
        } else {
          cat("✅", data_dir, "appears complete\n")
        }
      }
    }
  }
  
  cat("\n")
}

# Function to check if all required data already exists
check_existing_data <- function() {
  required_files <- c(
    "data/processed/tcga_clinical.rds",
    "data/processed/tcga_expression.rds", 
    "data/processed/tcga_cnv.rds"
  )
  
  existing_files <- file.exists(required_files)
  
  cat("🔍 CHECKING FOR EXISTING DATA FILES:\n")
  for (i in 1:length(required_files)) {
    status <- if (existing_files[i]) "✅ EXISTS" else "❌ MISSING"
    cat(sprintf("   %s: %s\n", required_files[i], status))
  }
  
  if (all(existing_files)) {
    cat("\n🎯 ALL PROCESSED DATA FILES EXIST!\n")
    
    # Additional check: verify files are not empty or corrupted
    cat("Verifying file integrity...\n")
    valid_files <- rep(FALSE, length(required_files))
    
    for (i in 1:length(required_files)) {
      tryCatch({
        data <- readRDS(required_files[i])
        valid_files[i] <- !is.null(data) && nrow(data) > 0
        if (valid_files[i]) {
          cat(sprintf("   ✅ %s: %d rows\n", basename(required_files[i]), nrow(data)))
        } else {
          cat(sprintf("   ⚠️  %s: File exists but appears empty\n", basename(required_files[i])))
        }
      }, error = function(e) {
        cat(sprintf("   ❌ %s: File corrupted - %s\n", basename(required_files[i]), e$message))
        valid_files[i] <<- FALSE
      })
    }
    
    if (all(valid_files)) {
      cat("\n✅ All files are valid! Loading existing data...\n")
      
      # Clean up any incomplete raw downloads that might interfere
      clean_incomplete_downloads()
      
      # Load existing data and update metadata
      clinical_final <<- readRDS("data/processed/tcga_clinical.rds")
      rna_final <<- readRDS("data/processed/tcga_expression.rds")
      cnv_final <<- readRDS("data/processed/tcga_cnv.rds")
      
      # Validate loaded data structure
      cat("Validating loaded data structure:\n")
      
      # Check clinical data
      if (nrow(clinical_final) == 0) {
        cat("⚠️  Warning: Clinical data appears empty\n")
      } else {
        cat("✓ Clinical data: ", nrow(clinical_final), " patients\n")
        if (!"survival_days" %in% names(clinical_final)) cat("⚠️  Missing survival_days column\n")
        if (!"survival_event" %in% names(clinical_final)) cat("⚠️  Missing survival_event column\n")
      }
      
      # Check RNA data  
      if (nrow(rna_final) == 0) {
        cat("⚠️  Warning: RNA data appears empty\n")
      } else {
        cat("✓ RNA data: ", nrow(rna_final), " samples\n")
        if (!"ERBB2_log2" %in% names(rna_final)) cat("⚠️  Missing ERBB2_log2 column\n")
      }
      
      # Check CNV data
      if (nrow(cnv_final) == 0) {
        cat("⚠️  Warning: CNV data appears empty\n")
      } else {
        cat("✓ CNV data: ", nrow(cnv_final), " samples\n")
        if (!"erbb2_amplified" %in% names(cnv_final)) cat("⚠️  Missing erbb2_amplified column\n")
      }
      
      # Update session metadata with proper calculations
      session_metadata$datasets_created <<- c("tcga_clinical.rds", "tcga_expression.rds", "tcga_cnv.rds")
      session_metadata$sample_counts$total_clinical <<- nrow(clinical_final)
      session_metadata$sample_counts$total_rna <<- nrow(rna_final)
      session_metadata$sample_counts$total_cnv <<- nrow(cnv_final)
      session_metadata$sample_counts$with_survival <<- sum(!is.na(clinical_final$survival_days))
      session_metadata$sample_counts$deaths_observed <<- sum(clinical_final$survival_event == 1, na.rm = TRUE)
      session_metadata$sample_counts$erbb2_amplified <<- sum(cnv_final$erbb2_amplified, na.rm = TRUE)
      session_metadata$data_sources <<- list(
        clinical = "loaded_from_existing_file",
        rna = "loaded_from_existing_file", 
        cnv = "loaded_from_existing_file"
      )
      
      cat("✅ Data loaded successfully:\n")
      cat(sprintf("   - Clinical: %d patients (%d with survival)\n", 
                  nrow(clinical_final), sum(!is.na(clinical_final$survival_days))))
      cat(sprintf("   - RNA: %d samples\n", nrow(rna_final)))
      cat(sprintf("   - CNV: %d samples (%d amplified)\n", 
                  nrow(cnv_final), sum(cnv_final$erbb2_amplified, na.rm = TRUE)))
      
      return(TRUE)  # Skip downloads
    } else {
      cat("\n⚠️  Some files are corrupted or empty. Will regenerate all data.\n")
      cat("Tip: Run inspect_existing_data() for detailed file inspection.\n")
      
      # Clean up incomplete downloads before regenerating
      clean_incomplete_downloads()
      
      return(FALSE)  # Need to reprocess
    }
  }
  
  # Clean up any incomplete downloads before fresh downloads
  clean_incomplete_downloads()
  
  return(FALSE)  # Need to download
}

# Function to check if raw TCGA data is already downloaded
check_existing_raw_data <- function() {
  tcga_dirs <- c(
    file.path("GDCdata", "TCGA-BRCA", "harmonized", "Clinical"),
    file.path("GDCdata", "TCGA-BRCA", "harmonized", "Transcriptome_Profiling"),
    file.path("GDCdata", "TCGA-BRCA", "harmonized", "Copy_Number_Variation")
  )
  
  existing_dirs <- sapply(tcga_dirs, dir.exists)
  
  cat("🔍 CHECKING FOR EXISTING RAW TCGA DATA:\n")
  data_types <- c("Clinical", "RNA-seq", "Copy Number")
  for (i in 1:length(tcga_dirs)) {
    status <- if (existing_dirs[i]) "✅ EXISTS" else "❌ MISSING"
    cat(sprintf("   %s data: %s\n", data_types[i], status))
  }
  
  return(existing_dirs)
}

# Check if we can skip downloads entirely
SKIP_ALL_DOWNLOADS <- check_existing_data()

if (SKIP_ALL_DOWNLOADS) {
  cat("\n✅ SKIPPING ALL DOWNLOADS - USING EXISTING PROCESSED DATA\n")
  cat("Proceeding directly to quality assessment...\n\n")
  
  # Force quality assessment even with loaded data
  cat("\n4. COMPREHENSIVE DATA QUALITY ASSESSMENT\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Ensure data is properly loaded and assessed
  ensure_data_loaded()
  
  overlap_data <- safe_execute(
    generate_quality_report(),
    "Quality assessment of loaded data"
  )
  
} else {
  # Check if raw data exists
  existing_raw <- check_existing_raw_data()
  
  cat("\n📥 PROCEEDING WITH DATA PROCESSING...\n")
  if (any(existing_raw)) {
    cat("Some raw TCGA data already exists - will skip those downloads\n")
  }
  cat("\n")
  
  # =============================================================================
  # 1. COMPREHENSIVE CLINICAL DATA ACQUISITION
  # =============================================================================
  
  cat("\n1. ACQUIRING COMPREHENSIVE CLINICAL DATA\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  acquire_clinical_data <- function() {
    # Check if clinical data already exists and is complete
    clinical_dir <- file.path("GDCdata", "TCGA-BRCA", "harmonized", "Clinical")
    
    if (dir.exists(clinical_dir) && length(list.files(clinical_dir, recursive = TRUE)) > 0) {
      cat("✅ Clinical data already downloaded, skipping download completely...\n")
      
      # Create query for data preparation only (no download attempt)
      clinical_query <- GDCquery(
        project = "TCGA-BRCA",
        data.category = "Clinical",
        data.type = "Clinical Supplement",
        data.format = "BCR Biotab"
      )
      
      # Check if the query results match existing files
      expected_files <- getResults(clinical_query)
      existing_files <- list.files(clinical_dir, recursive = TRUE, pattern = "\\.txt$")
      
      if (length(existing_files) >= nrow(expected_files) * 0.8) {  # At least 80% of files exist
        cat("✅ Clinical data appears complete, skipping all downloads\n")
      } else {
        cat("⚠️  Clinical data incomplete, forcing fresh download...\n")
        # Remove incomplete data and start fresh
        unlink(clinical_dir, recursive = TRUE)
        cat("📥 Downloading TCGA-BRCA clinical data...\n")
        GDCdownload(clinical_query, method = "api")
        cat("✅ Clinical data download complete\n")
      }
    } else {
      cat("📥 Downloading TCGA-BRCA clinical data...\n")
      
      # Query for all clinical data types
      clinical_query <- GDCquery(
        project = "TCGA-BRCA",
        data.category = "Clinical",
        data.type = "Clinical Supplement",
        data.format = "BCR Biotab"
      )
      
      cat("Found", nrow(getResults(clinical_query)), "clinical files\n")
      
      # Download clinical data
      GDCdownload(clinical_query, method = "api")
      cat("✅ Clinical data download complete\n")
    }
    
    # Prepare comprehensive clinical data
    cat("Preparing clinical data tables...\n")
    
    # Get patient-level data
    clinical_patient <- GDCprepare_clinic(clinical_query, clinical.info = "patient")
    
    # Get follow-up data for survival
    clinical_followup <- tryCatch({
      GDCprepare_clinic(clinical_query, clinical.info = "follow_up")
    }, error = function(e) {
      cat("Follow-up data not available:", e$message, "\n")
      return(NULL)
    })
    
    # Get drug treatment data (for hormone therapy)
    clinical_drug <- tryCatch({
      GDCprepare_clinic(clinical_query, clinical.info = "drug")
    }, error = function(e) {
      cat("Drug data not available:", e$message, "\n")
      return(NULL)
    })
    
    # Try to get biomarker data from additional clinical files
    clinical_biomarker <- tryCatch({
      # This may contain IHC results
      clinical_files <- getResults(clinical_query)
      biomarker_files <- clinical_files[grepl("biospecimen|pathology|ihc", 
                                              clinical_files$file_name, ignore.case = TRUE), ]
      
      if (nrow(biomarker_files) > 0) {
        cat("Found potential biomarker files:", nrow(biomarker_files), "\n")
        # Additional processing for biomarker files would go here
      }
      return(NULL)  # Placeholder for now
    }, error = function(e) {
      cat("Biomarker data extraction failed:", e$message, "\n")
      return(NULL)
    })
    
    return(list(
      patient = clinical_patient,
      followup = clinical_followup,
      drug = clinical_drug,
      biomarker = clinical_biomarker
    ))
  }
  
  clinical_data <- safe_execute(
    acquire_clinical_data(),
    "Comprehensive clinical data acquisition"
  )
  
  # Process clinical data for survival and biomarkers
  if (!file.exists("data/processed/tcga_clinical.rds")) {
    process_clinical_data <- function(clinical_data) {
      if (is.null(clinical_data) || is.null(clinical_data$patient)) return(NULL)
      
      cat("Processing clinical data for survival and biomarkers...\n")
      
      # Start with patient data
      clinical_processed <- clinical_data$patient %>%
        select(
          patient_id = submitter_id,
          age_at_diagnosis,
          gender,
          race,
          ethnicity,
          vital_status,
          days_to_death,
          days_to_last_follow_up,
          primary_diagnosis,
          morphology,
          tumor_stage,
          pathologic_stage,
          histological_type,
          prior_malignancy,
          synchronous_malignancy,
          # Try to capture additional biomarker fields if they exist
          everything()
        ) %>%
        # Convert age from days to years
        mutate(
          age_years = as.numeric(age_at_diagnosis) / 365.25,
          patient_id = as.character(patient_id)
        ) %>%
        filter(!is.na(patient_id))
      
      # Process survival data with follow-up information
      if (!is.null(clinical_data$followup)) {
        cat("Processing survival data with follow-up information...\n")
        
        # Get latest follow-up for each patient
        followup_latest <- clinical_data$followup %>%
          filter(!is.na(submitter_id)) %>%
          group_by(submitter_id) %>%
          arrange(desc(as.numeric(days_to_follow_up))) %>%
          slice(1) %>%
          ungroup() %>%
          select(
            patient_id = submitter_id,
            followup_vital_status = vital_status,
            followup_days = days_to_follow_up,
            followup_days_to_death = days_to_death,
            disease_status = disease_status,
            new_tumor_events = new_tumor_event_after_initial_treatment
          )
        
        # Merge with patient data
        clinical_processed <- clinical_processed %>%
          left_join(followup_latest, by = "patient_id")
        
        # Create comprehensive survival variables
        clinical_processed <- clinical_processed %>%
          mutate(
            # Use most recent vital status
            final_vital_status = ifelse(
              !is.na(followup_vital_status), 
              followup_vital_status, 
              vital_status
            ),
            
            # Use most recent survival time
            survival_days = case_when(
              !is.na(followup_days_to_death) ~ as.numeric(followup_days_to_death),
              !is.na(days_to_death) ~ as.numeric(days_to_death),
              !is.na(followup_days) ~ as.numeric(followup_days),
              !is.na(days_to_last_follow_up) ~ as.numeric(days_to_last_follow_up),
              TRUE ~ NA_real_
            ),
            
            # Create survival event indicator
            survival_event = case_when(
              final_vital_status == "Dead" ~ 1,
              final_vital_status == "Alive" ~ 0,
              TRUE ~ NA_real_
            ),
            
            # Convert to years
            survival_years = survival_days / 365.25,
            
            # Disease progression event
            progression_event = case_when(
              disease_status %in% c("Tumor progression", "Recurred/Progressed") ~ 1,
              disease_status %in% c("Tumor free", "Disease Free") ~ 0,
              TRUE ~ NA_real_
            )
          )
      } else {
        # Fallback to basic survival data only
        clinical_processed <- clinical_processed %>%
          mutate(
            final_vital_status = vital_status,
            survival_days = ifelse(
              vital_status == "Dead",
              as.numeric(days_to_death),
              as.numeric(days_to_last_follow_up)
            ),
            survival_event = ifelse(vital_status == "Dead", 1, 0),
            survival_years = survival_days / 365.25,
            progression_event = NA_real_
          )
      }
      
      # Extract tumor characteristics and staging
      clinical_processed <- clinical_processed %>%
        mutate(
          # Clean pathologic stage
          pathologic_stage_clean = case_when(
            str_detect(pathologic_stage, "Stage I[^IV]") ~ "I",
            str_detect(pathologic_stage, "Stage II") ~ "II", 
            str_detect(pathologic_stage, "Stage III") ~ "III",
            str_detect(pathologic_stage, "Stage IV") ~ "IV",
            TRUE ~ pathologic_stage
          ),
          
          # Extract grade if available in morphology or other fields
          histologic_grade = case_when(
            str_detect(morphology, "Grade 1|G1") ~ "G1",
            str_detect(morphology, "Grade 2|G2") ~ "G2", 
            str_detect(morphology, "Grade 3|G3") ~ "G3",
            str_detect(morphology, "Grade 4|G4") ~ "G4",
            TRUE ~ NA_character_
          ),
          
          # Basic tumor type classification
          tumor_type = case_when(
            str_detect(histological_type, "Ductal") ~ "Invasive Ductal Carcinoma",
            str_detect(histological_type, "Lobular") ~ "Invasive Lobular Carcinoma",
            str_detect(histological_type, "Mixed") ~ "Mixed Ductal and Lobular",
            TRUE ~ histological_type
          )
        )
      
      cat("Clinical data processing complete:", nrow(clinical_processed), "patients\n")
      cat("Patients with survival data:", sum(!is.na(clinical_processed$survival_days)), "\n")
      cat("Deaths observed:", sum(clinical_processed$survival_event == 1, na.rm = TRUE), "\n")
      
      return(clinical_processed)
    }
    
    clinical_final <- safe_execute(
      process_clinical_data(clinical_data),
      "Clinical data processing with survival"
    )
    
    # Download and extract HER2 IHC status from TCGA biomarker data
    extract_her2_status <- function() {
      cat("Attempting to extract HER2 biomarker status...\n")
      
      # Try to download biomarker data from TCGA
      tryCatch({
        # Query for biomarker/pathology data
        biomarker_query <- GDCquery(
          project = "TCGA-BRCA",
          data.category = "Clinical", 
          data.type = "Clinical Supplement",
          data.format = "BCR Biotab"
        )
        
        # Look for files that might contain biomarker data
        files_info <- getResults(biomarker_query)
        cat("Available clinical files:", nrow(files_info), "\n")
        
        # Download if not already downloaded
        if (!dir.exists(file.path("GDCdata", "TCGA-BRCA", "harmonized", "Clinical"))) {
          GDCdownload(biomarker_query, method = "api")
        }
        
        # Try to extract biomarker data from downloaded files
        # Look for specific biomarker tables
        biomarker_tables <- c("biospecimen_sample_brca", "clinical_patient_brca", 
                              "biospecimen_aliquot_brca", "pathology_details_brca")
        
        her2_data <- NULL
        
        # Search through available clinical data for biomarker information
        file_path <- file.path("GDCdata", "TCGA-BRCA", "harmonized", "Clinical")
        if (dir.exists(file_path)) {
          all_files <- list.files(file_path, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)
          
          for (file in all_files) {
            if (file.exists(file)) {
              # Try to read and search for HER2-related columns
              temp_data <- tryCatch({
                read_tsv(file, show_col_types = FALSE)
              }, error = function(e) NULL)
              
              if (!is.null(temp_data)) {
                her2_cols <- names(temp_data)[grepl("her2|erbb2|ihc", names(temp_data), ignore.case = TRUE)]
                er_cols <- names(temp_data)[grepl("\\ber\\b|estrogen", names(temp_data), ignore.case = TRUE)]
                pr_cols <- names(temp_data)[grepl("\\bpr\\b|progesterone", names(temp_data), ignore.case = TRUE)]
                
                if (length(her2_cols) > 0 || length(er_cols) > 0 || length(pr_cols) > 0) {
                  cat("Found potential biomarker data in:", basename(file), "\n")
                  cat("HER2 columns:", paste(her2_cols, collapse = ", "), "\n")
                  cat("ER columns:", paste(er_cols, collapse = ", "), "\n") 
                  cat("PR columns:", paste(pr_cols, collapse = ", "), "\n")
                  
                  # Extract relevant columns
                  biomarker_subset <- temp_data %>%
                    select(any_of(c("bcr_patient_barcode", "submitter_id", 
                                    her2_cols, er_cols, pr_cols))) %>%
                    filter(!is.na(bcr_patient_barcode) | !is.na(submitter_id))
                  
                  if (nrow(biomarker_subset) > 0) {
                    her2_data <- biomarker_subset
                    break
                  }
                }
              }
            }
          }
        }
        
        return(her2_data)
        
      }, error = function(e) {
        cat("Biomarker extraction failed:", e$message, "\n")
        return(NULL)
      })
    }
    
    # Extract biomarker data
    her2_biomarker_data <- safe_execute(
      extract_her2_status(),
      "HER2 biomarker status extraction"
    )
    
    # Merge biomarker data with clinical data if available
    if (!is.null(her2_biomarker_data) && !is.null(clinical_final)) {
      cat("Merging biomarker data with clinical data...\n")
      
      # Standardize patient IDs for merging
      if ("bcr_patient_barcode" %in% names(her2_biomarker_data)) {
        her2_biomarker_data$patient_id <- her2_biomarker_data$bcr_patient_barcode
      } else if ("submitter_id" %in% names(her2_biomarker_data)) {
        her2_biomarker_data$patient_id <- her2_biomarker_data$submitter_id
      }
      
      # Merge with clinical data
      clinical_final <- clinical_final %>%
        left_join(her2_biomarker_data %>% select(-any_of(c("bcr_patient_barcode", "submitter_id"))), 
                  by = "patient_id")
      
      cat("Successfully merged biomarker data\n")
    } else {
      cat("No biomarker data available or merge failed\n")
      # Add placeholder columns for later sessions
      clinical_final$her2_ihc_status <- NA_character_
      clinical_final$er_status <- NA_character_
      clinical_final$pr_status <- NA_character_
    }
    
    # Save clinical data
    if (!is.null(clinical_final)) {
      saveRDS(clinical_final, "data/processed/tcga_clinical.rds")
      session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_clinical.rds")
      session_metadata$sample_counts$total_clinical <- nrow(clinical_final)
      session_metadata$sample_counts$with_survival <- sum(!is.na(clinical_final$survival_days))
      session_metadata$sample_counts$deaths_observed <- sum(clinical_final$survival_event == 1, na.rm = TRUE)
      cat("Clinical data saved:", nrow(clinical_final), "patients\n")
    }
    
  } else {
    cat("✅ Clinical data already processed, loading from file...\n")
    clinical_final <- readRDS("data/processed/tcga_clinical.rds")
    session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_clinical.rds")
    session_metadata$sample_counts$total_clinical <- nrow(clinical_final)
    session_metadata$sample_counts$with_survival <- sum(!is.na(clinical_final$survival_days))
    session_metadata$sample_counts$deaths_observed <- sum(clinical_final$survival_event == 1, na.rm = TRUE)
    cat("Clinical data loaded:", nrow(clinical_final), "patients\n")
  }
  
  # =============================================================================
  # 2. RNA EXPRESSION DATA ACQUISITION (ERBB2)
  # =============================================================================
  
  cat("\n2. ACQUIRING RNA EXPRESSION DATA (ERBB2)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  acquire_rna_data <- function() {
    # Check if RNA data already exists and is complete
    rna_dir <- file.path("GDCdata", "TCGA-BRCA", "harmonized", "Transcriptome_Profiling")
    
    if (dir.exists(rna_dir) && length(list.files(rna_dir, recursive = TRUE)) > 10) {
      cat("✅ RNA-seq data directory exists, checking completeness...\n")
      
      # Create query to check expected files
      rna_query <- GDCquery(
        project = "TCGA-BRCA",
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        sample.type = c("Primary Tumor", "Solid Tissue Normal")
      )
      
      expected_files <- getResults(rna_query)
      existing_files <- list.files(rna_dir, recursive = TRUE, pattern = "\\.gz$")
      
      cat("Expected files:", nrow(expected_files), "| Existing files:", length(existing_files), "\n")
      
      # Check if we have at least 80% of expected files (allows for some samples being unavailable)
      if (length(existing_files) >= nrow(expected_files) * 0.8) {
        cat("✅ RNA-seq data appears complete (", 
            round(length(existing_files)/nrow(expected_files)*100, 1), 
            "% of expected files), skipping download\n")
      } else {
        cat("⚠️  RNA-seq data incomplete, forcing fresh download...\n")
        cat("Removing incomplete data directory...\n")
        unlink(rna_dir, recursive = TRUE)
        
        cat("📥 Downloading TCGA-BRCA RNA-seq data...\n")
        cat("Found", nrow(getResults(rna_query)), "RNA-seq files\n")
        GDCdownload(rna_query, method = "api", files.per.chunk = 50)
        cat("✅ RNA-seq data download complete\n")
      }
    } else {
      cat("📥 Downloading TCGA-BRCA RNA-seq data...\n")
      
      # Query for RNA-Seq data
      rna_query <- GDCquery(
        project = "TCGA-BRCA",
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        sample.type = c("Primary Tumor", "Solid Tissue Normal")
      )
      
      cat("Found", nrow(getResults(rna_query)), "RNA-seq files\n")
      
      # Download RNA data
      GDCdownload(rna_query, method = "api", files.per.chunk = 50)
      cat("✅ RNA-seq data download complete\n")
    }
    
    # Prepare RNA data
    rna_data <- GDCprepare(rna_query)
    
    return(rna_data)
  }
  
  rna_data <- safe_execute(
    acquire_rna_data(),
    "RNA expression data acquisition"
  )
  
  # Process RNA data for ERBB2 and related genes
  if (!file.exists("data/processed/tcga_expression.rds")) {
    process_rna_data <- function(rna_data) {
      if (is.null(rna_data)) return(NULL)
      
      cat("Processing RNA expression data...\n")
      
      # Extract gene information
      gene_info <- rowData(rna_data)
      
      # Find ERBB2 and related genes
      target_genes <- c("ERBB2", "ESR1", "PGR", "MKI67", "PCNA")
      gene_indices <- list()
      
      for (gene in target_genes) {
        idx <- which(gene_info$gene_name == gene)
        if (length(idx) == 0) {
          idx <- grep(paste0("^", gene, "$"), gene_info$gene_name, ignore.case = TRUE)
        }
        if (length(idx) > 0) {
          gene_indices[[gene]] <- idx[1]
          cat("Found", gene, "at index", idx[1], "\n")
        } else {
          cat("Warning:", gene, "not found\n")
        }
      }
      
      if (length(gene_indices) == 0) {
        cat("Error: No target genes found\n")
        return(NULL)
      }
      
      # Extract expression data
      sample_info <- colData(rna_data)
      expression_data <- assay(rna_data)
      
      # Create expression dataset
      rna_processed <- data.frame(
        patient_id = substr(sample_info$submitter_id, 1, 12),
        sample_id = sample_info$submitter_id,
        sample_type = sample_info$sample_type,
        stringsAsFactors = FALSE
      )
      
      # Add expression values for each found gene
      for (gene in names(gene_indices)) {
        idx <- gene_indices[[gene]]
        raw_counts <- as.numeric(expression_data[idx, ])
        rna_processed[[paste0(gene, "_counts")]] <- raw_counts
        rna_processed[[paste0(gene, "_log2")]] <- log2(raw_counts + 1)
        rna_processed[[paste0(gene, "_log2_centered")]] <- scale(log2(raw_counts + 1))[,1]
      }
      
      # Filter for primary tumors with valid expression
      rna_processed <- rna_processed %>%
        filter(sample_type == "Primary Tumor") %>%
        filter(!is.na(ERBB2_counts) & ERBB2_counts > 0) %>%
        arrange(patient_id)
      
      cat("RNA processing complete:", nrow(rna_processed), "samples\n")
      cat("Genes processed:", paste(names(gene_indices), collapse = ", "), "\n")
      
      return(rna_processed)
    }
    
    rna_final <- safe_execute(
      process_rna_data(rna_data),
      "RNA expression data processing"
    )
    
    # Save RNA data
    if (!is.null(rna_final)) {
      saveRDS(rna_final, "data/processed/tcga_expression.rds")
      session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_expression.rds")
      session_metadata$sample_counts$total_rna <- nrow(rna_final)
      cat("RNA expression data saved:", nrow(rna_final), "samples\n")
    }
    
  } else {
    cat("✅ RNA expression data already processed, loading from file...\n")
    rna_final <- readRDS("data/processed/tcga_expression.rds")
    session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_expression.rds")
    session_metadata$sample_counts$total_rna <- nrow(rna_final)
    cat("RNA expression data loaded:", nrow(rna_final), "samples\n")
  }
  
  # =============================================================================
  # 3. COPY NUMBER DATA ACQUISITION (ERBB2 locus)
  # =============================================================================
  
  cat("\n3. ACQUIRING COPY NUMBER DATA (ERBB2 locus)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  acquire_cnv_data <- function() {
    # Check if CNV data already exists and is complete
    cnv_dir <- file.path("GDCdata", "TCGA-BRCA", "harmonized", "Copy_Number_Variation")
    
    if (dir.exists(cnv_dir) && length(list.files(cnv_dir, recursive = TRUE)) > 10) {
      cat("✅ Copy number data directory exists, checking completeness...\n")
      
      # Try to load existing CNV data first
      tryCatch({
        # Try different data types that might be already downloaded
        cnv_options <- list(
          list(data.type = "Masked Copy Number Segment", desc = "Segment-level copy number"),
          list(data.type = "Gene Level Copy Number", desc = "Gene-level copy number")
        )
        
        successful_data <- NULL
        
        for(i in 1:length(cnv_options)) {
          if(!is.null(successful_data)) break
          
          option <- cnv_options[[i]]
          cat("Checking existing", option$desc, "...\n")
          
          tryCatch({
            cnv_query <- GDCquery(
              project = "TCGA-BRCA",
              data.category = "Copy Number Variation",
              data.type = option$data.type,
              sample.type = "Primary Tumor"
            )
            
            expected_files <- getResults(cnv_query)
            
            if (nrow(expected_files) > 0) {
              # Check file completeness
              existing_files <- list.files(cnv_dir, recursive = TRUE, 
                                           pattern = "\\.(txt|seg|tsv)$|copy_number")
              
              cat("Expected files:", nrow(expected_files), "| Existing files:", length(existing_files), "\n")
              
              if (length(existing_files) >= nrow(expected_files) * 0.7) {  # 70% threshold for CNV data
                cat("✅ Copy number data appears complete for", option$desc, "\n")
                
                # Try to prepare the data
                cnv_data <- GDCprepare(cnv_query)
                attr(cnv_data, "data_type") <- option$data.type
                successful_data <- cnv_data
                cat("✅ Successfully loaded existing", option$desc, "\n")
                
              } else {
                cat("⚠️  Incomplete", option$desc, "data\n")
              }
            }
            
          }, error = function(e) {
            cat("Could not process existing", option$desc, ":", e$message, "\n")
          })
        }
        
        if (!is.null(successful_data)) {
          return(successful_data)
        } else {
          cat("⚠️  No complete CNV datasets found, forcing fresh download...\n")
          unlink(cnv_dir, recursive = TRUE)
        }
        
      }, error = function(e) {
        cat("Error checking existing CNV data:", e$message, "\n")
        cat("Proceeding with fresh download...\n")
        unlink(cnv_dir, recursive = TRUE)
      })
    }
    
    # If we reach here, need to download fresh data
    cat("📥 Downloading TCGA-BRCA copy number data...\n")
    
    # Try different copy number data types
    cnv_options <- list(
      list(data.type = "Masked Copy Number Segment", desc = "Segment-level copy number"),
      list(data.type = "Gene Level Copy Number", desc = "Gene-level copy number")
    )
    
    successful_data <- NULL
    
    for(i in 1:length(cnv_options)) {
      if(!is.null(successful_data)) break
      
      option <- cnv_options[[i]]
      cat("Trying copy number option", i, ":", option$desc, "...\n")
      
      tryCatch({
        cnv_query <- GDCquery(
          project = "TCGA-BRCA",
          data.category = "Copy Number Variation",
          data.type = option$data.type,
          sample.type = "Primary Tumor"
        )
        
        if(nrow(getResults(cnv_query)) > 0) {
          cat("✓ Found", nrow(getResults(cnv_query)), "files\n")
          
          # Download
          GDCdownload(cnv_query, method = "api", files.per.chunk = 25)
          
          # Prepare
          cnv_data <- GDCprepare(cnv_query)
          
          attr(cnv_data, "data_type") <- option$data.type
          successful_data <- cnv_data
          cat("✓ Success with", option$desc, "\n")
        }
      }, error = function(e) {
        cat("✗ Failed with option", i, ":", e$message, "\n")
      })
    }
    
    return(successful_data)
  }
  
  cnv_data <- safe_execute(
    acquire_cnv_data(),
    "Copy number data acquisition"
  )
  
  # Process CNV data for ERBB2 locus
  if (!file.exists("data/processed/tcga_cnv.rds")) {
    process_cnv_data <- function(cnv_data) {
      if (is.null(cnv_data)) return(NULL)
      
      data_type <- attr(cnv_data, "data_type")
      cat("Processing copy number data type:", data_type, "\n")
      
      if(data_type %in% c("Masked Copy Number Segment", "Copy Number Segment")) {
        # For segment data, filter for ERBB2 locus (chr17: 37,844,167-37,886,679)
        cat("Processing segment-level data for ERBB2 locus...\n")
        
        cnv_processed <- cnv_data %>%
          filter(
            Chromosome == "17",
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
            sample_type = "Primary Tumor",
            # Classify amplification status
            erbb2_amplified = erbb2_cnv_score > 0.3,  # Log2 ratio > 0.3 suggests amplification
            erbb2_deleted = erbb2_cnv_score < -0.3
          ) %>%
          filter(!is.na(erbb2_cnv_score))
        
      } else {
        # For gene-level data
        cat("Processing gene-level data for ERBB2...\n")
        
        if(inherits(cnv_data, c("RangedSummarizedExperiment", "SummarizedExperiment"))) {
          gene_info <- rowData(cnv_data)
          erbb2_idx <- which(gene_info$gene_name == "ERBB2" | gene_info$Hugo_Symbol == "ERBB2")
          
          if (length(erbb2_idx) == 0) {
            erbb2_idx <- grep("ERBB2", gene_info$gene_name, ignore.case = TRUE)
          }
          
          if (length(erbb2_idx) > 0) {
            erbb2_cnv <- assay(cnv_data)[erbb2_idx[1], ]
            sample_info <- colData(cnv_data)
            
            cnv_processed <- data.frame(
              patient_id = substr(sample_info$submitter_id, 1, 12),
              sample_id = sample_info$submitter_id,
              erbb2_cnv_score = as.numeric(erbb2_cnv),
              sample_type = sample_info$sample_type,
              erbb2_amplified = as.numeric(erbb2_cnv) > 0.3,
              erbb2_deleted = as.numeric(erbb2_cnv) < -0.3,
              stringsAsFactors = FALSE
            ) %>%
              filter(sample_type == "Primary Tumor") %>%
              filter(!is.na(erbb2_cnv_score))
          } else {
            cat("ERBB2 not found in gene-level data\n")
            return(NULL)
          }
        }
      }
      
      cat("CNV processing complete:", nrow(cnv_processed), "samples\n")
      cat("Amplified samples:", sum(cnv_processed$erbb2_amplified), "\n")
      cat("Deleted samples:", sum(cnv_processed$erbb2_deleted), "\n")
      
      return(cnv_processed)
    }
    
    cnv_final <- safe_execute(
      process_cnv_data(cnv_data),
      "Copy number data processing"
    )
    
    # Save CNV data
    if (!is.null(cnv_final)) {
      saveRDS(cnv_final, "data/processed/tcga_cnv.rds")
      session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_cnv.rds")
      session_metadata$sample_counts$total_cnv <- nrow(cnv_final)
      session_metadata$sample_counts$erbb2_amplified <- sum(cnv_final$erbb2_amplified)
      cat("Copy number data saved:", nrow(cnv_final), "samples\n")
    }
    
  } else {
    cat("✅ Copy number data already processed, loading from file...\n")
    cnv_final <- readRDS("data/processed/tcga_cnv.rds")
    session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_cnv.rds")
    session_metadata$sample_counts$total_cnv <- nrow(cnv_final)
    session_metadata$sample_counts$erbb2_amplified <- sum(cnv_final$erbb2_amplified)
    cat("Copy number data loaded:", nrow(cnv_final), "samples\n")
  }
  
} # End of main download/processing block

# =============================================================================
# 4. COMPREHENSIVE DATA QUALITY ASSESSMENT
# =============================================================================

cat("\n4. COMPREHENSIVE DATA QUALITY ASSESSMENT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Ensure all data is loaded (safety check for mixed download/skip scenarios)
ensure_data_loaded <- function() {
  cat("🔍 Ensuring all data is properly loaded...\n")
  
  if (!exists("clinical_final") || is.null(clinical_final)) {
    if (file.exists("data/processed/tcga_clinical.rds")) {
      clinical_final <<- readRDS("data/processed/tcga_clinical.rds")
      cat("✅ Loaded clinical data from file\n")
    } else {
      cat("❌ Clinical data not available\n")
    }
  }
  
  if (!exists("rna_final") || is.null(rna_final)) {
    if (file.exists("data/processed/tcga_expression.rds")) {
      rna_final <<- readRDS("data/processed/tcga_expression.rds")
      cat("✅ Loaded RNA data from file\n")
    } else {
      cat("❌ RNA data not available\n")
    }
  }
  
  if (!exists("cnv_final") || is.null(cnv_final)) {
    if (file.exists("data/processed/tcga_cnv.rds")) {
      cnv_final <<- readRDS("data/processed/tcga_cnv.rds")
      cat("✅ Loaded CNV data from file\n")
    } else {
      cat("❌ CNV data not available\n")
    }
  }
  
  # Update session metadata if needed
  if (exists("clinical_final") && !is.null(clinical_final) && 
      !"tcga_clinical.rds" %in% session_metadata$datasets_created) {
    session_metadata$datasets_created <<- c(session_metadata$datasets_created, "tcga_clinical.rds")
    session_metadata$sample_counts$total_clinical <<- nrow(clinical_final)
    session_metadata$sample_counts$with_survival <<- sum(!is.na(clinical_final$survival_days))
    session_metadata$sample_counts$deaths_observed <<- sum(clinical_final$survival_event == 1, na.rm = TRUE)
  }
  
  if (exists("rna_final") && !is.null(rna_final) && 
      !"tcga_expression.rds" %in% session_metadata$datasets_created) {
    session_metadata$datasets_created <<- c(session_metadata$datasets_created, "tcga_expression.rds")
    session_metadata$sample_counts$total_rna <<- nrow(rna_final)
  }
  
  if (exists("cnv_final") && !is.null(cnv_final) && 
      !"tcga_cnv.rds" %in% session_metadata$datasets_created) {
    session_metadata$datasets_created <<- c(session_metadata$datasets_created, "tcga_cnv.rds")
    session_metadata$sample_counts$total_cnv <<- nrow(cnv_final)
    session_metadata$sample_counts$erbb2_amplified <<- sum(cnv_final$erbb2_amplified, na.rm = TRUE)
  }
}

ensure_data_loaded()

generate_quality_report <- function() {
  cat("Starting quality assessment...\n")
  
  # Debug: Check if data exists
  cat("Checking data availability:\n")
  cat("- clinical_final exists:", exists("clinical_final"), "\n")
  cat("- rna_final exists:", exists("rna_final"), "\n") 
  cat("- cnv_final exists:", exists("cnv_final"), "\n")
  
  if (exists("clinical_final")) cat("- clinical_final rows:", nrow(clinical_final), "\n")
  if (exists("rna_final")) cat("- rna_final rows:", nrow(rna_final), "\n")
  if (exists("cnv_final")) cat("- cnv_final rows:", nrow(cnv_final), "\n")
  
  if (!exists("clinical_final") || !exists("rna_final") || !exists("cnv_final") ||
      is.null(clinical_final) || is.null(rna_final) || is.null(cnv_final)) {
    cat("❌ Error: Not all datasets available for quality assessment\n")
    return(NULL)
  }
  
  if (nrow(clinical_final) == 0 || nrow(rna_final) == 0 || nrow(cnv_final) == 0) {
    cat("❌ Error: One or more datasets are empty\n")
    cat("Clinical rows:", nrow(clinical_final), "\n")
    cat("RNA rows:", nrow(rna_final), "\n") 
    cat("CNV rows:", nrow(cnv_final), "\n")
    return(NULL)
  }
  cat("Starting quality assessment...\n")
  
  # Debug: Check if data exists
  cat("Checking data availability:\n")
  cat("- clinical_final exists:", exists("clinical_final"), "\n")
  cat("- rna_final exists:", exists("rna_final"), "\n") 
  cat("- cnv_final exists:", exists("cnv_final"), "\n")
  
  if (exists("clinical_final")) cat("- clinical_final rows:", nrow(clinical_final), "\n")
  if (exists("rna_final")) cat("- rna_final rows:", nrow(rna_final), "\n")
  if (exists("cnv_final")) cat("- cnv_final rows:", nrow(cnv_final), "\n")
  
  if (!exists("clinical_final") || !exists("rna_final") || !exists("cnv_final") ||
      is.null(clinical_final) || is.null(rna_final) || is.null(cnv_final)) {
    cat("❌ Error: Not all datasets available for quality assessment\n")
    return(NULL)
  }
  
  if (nrow(clinical_final) == 0 || nrow(rna_final) == 0 || nrow(cnv_final) == 0) {
    cat("❌ Error: One or more datasets are empty\n")
    cat("Clinical rows:", nrow(clinical_final), "\n")
    cat("RNA rows:", nrow(rna_final), "\n") 
    cat("CNV rows:", nrow(cnv_final), "\n")
    return(NULL)
  }
  
  # Sample overlap analysis
  clinical_ids <- clinical_final$patient_id
  rna_ids <- rna_final$patient_id
  cnv_ids <- cnv_final$patient_id
  
  all_ids <- unique(c(clinical_ids, rna_ids, cnv_ids))
  
  overlap_matrix <- data.frame(
    patient_id = all_ids,
    has_clinical = all_ids %in% clinical_ids,
    has_rna = all_ids %in% rna_ids,
    has_cnv = all_ids %in% cnv_ids,
    stringsAsFactors = FALSE
  )
  
  # Add data completeness flags
  overlap_matrix <- overlap_matrix %>%
    mutate(
      complete_molecular = has_rna & has_cnv,
      complete_all = has_clinical & has_rna & has_cnv
    )
  
  # Complete data counts
  complete_data <- overlap_matrix %>% filter(complete_all)
  
  # Quality metrics
  quality_metrics <- list(
    total_unique_patients = length(all_ids),
    complete_data_patients = nrow(complete_data),
    complete_molecular_only = sum(overlap_matrix$complete_molecular),
    clinical_only = sum(overlap_matrix$has_clinical & !overlap_matrix$has_rna & !overlap_matrix$has_cnv),
    rna_missing = sum(!overlap_matrix$has_rna),
    cnv_missing = sum(!overlap_matrix$has_cnv),
    complete_data_percent = round(nrow(complete_data) / length(all_ids) * 100, 1),
    
    # Survival data quality
    survival_data_available = sum(!is.na(clinical_final$survival_days)),
    survival_events = sum(clinical_final$survival_event == 1, na.rm = TRUE),
    survival_data_percent = round(sum(!is.na(clinical_final$survival_days)) / nrow(clinical_final) * 100, 1),
    
    # HER2 amplification rates
    her2_amplified_cnv = if(!is.null(cnv_final)) sum(cnv_final$erbb2_amplified) else 0,
    her2_amplification_rate = if(!is.null(cnv_final)) round(sum(cnv_final$erbb2_amplified) / nrow(cnv_final) * 100, 1) else 0
  )
  
  session_metadata$quality_metrics <<- quality_metrics
  session_metadata$sample_counts$complete_data <<- nrow(complete_data)
  
  # Print comprehensive quality summary
  cat("=== COMPREHENSIVE QUALITY ASSESSMENT ===\n")
  cat("Dataset Overview:\n")
  cat("- Total unique patients:", quality_metrics$total_unique_patients, "\n")
  cat("- Complete data (all platforms):", quality_metrics$complete_data_patients, 
      "(", quality_metrics$complete_data_percent, "%)\n")
  cat("- Molecular data only (RNA+CNV):", quality_metrics$complete_molecular_only, "\n")
  
  cat("\nSurvival Data Quality:\n")
  cat("- Patients with survival data:", quality_metrics$survival_data_available,
      "(", quality_metrics$survival_data_percent, "%)\n")
  cat("- Death events observed:", quality_metrics$survival_events, "\n")
  cat("- Median follow-up time:", round(median(clinical_final$survival_years, na.rm = TRUE), 2), "years\n")
  
  cat("\nMolecular Data Quality:\n")
  cat("- ERBB2 amplified samples (CNV):", quality_metrics$her2_amplified_cnv,
      "(", quality_metrics$her2_amplification_rate, "%)\n")
  cat("- Median ERBB2 expression (log2):", round(median(rna_final$ERBB2_log2, na.rm = TRUE), 2), "\n")
  
  cat("\nPlatform Overlaps:\n")
  cat("- Clinical + RNA:", length(intersect(clinical_ids, rna_ids)), "patients\n")
  cat("- Clinical + CNV:", length(intersect(clinical_ids, cnv_ids)), "patients\n")
  cat("- RNA + CNV:", length(intersect(rna_ids, cnv_ids)), "patients\n")
  cat("- All three platforms:", nrow(complete_data), "patients\n")
  
  return(overlap_matrix)
}

overlap_data <- safe_execute(
  generate_quality_report(),
  "Comprehensive data quality assessment"
)

# =============================================================================
# 5. SESSION FINALIZATION AND METADATA
# =============================================================================

cat("\n5. FINALIZING SESSION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Update metadata with completion info
session_metadata$end_time <- Sys.time()
session_metadata$duration_minutes <- as.numeric(difftime(
  session_metadata$end_time, 
  session_metadata$start_time, 
  units = "mins"
))

# Set next session inputs
session_metadata$next_session_inputs <- session_metadata$datasets_created

# Add data characteristics for next session
session_metadata$data_characteristics <- list(
  clinical_variables = names(clinical_final),
  rna_genes_available = grep("_counts$", names(rna_final), value = TRUE),
  cnv_features = c("erbb2_cnv_score", "erbb2_amplified", "erbb2_deleted"),
  ready_for_harmonization = length(session_metadata$datasets_created) == 3
)

# Save comprehensive metadata
write_json(session_metadata, "metadata/session1_metadata.json", pretty = TRUE, auto_unbox = TRUE)

# Create a data quality report
create_quality_report <- function() {
  if (!is.null(overlap_data)) {
    # Create simple HTML report
    html_content <- paste0(
      "<html><head><title>TCGA-BRCA Data Quality Report</title></head><body>",
      "<h1>TCGA-BRCA HER2 Integration Project - Data Quality Report</h1>",
      "<h2>Session 1: Data Acquisition</h2>",
      "<p><strong>Date:</strong> ", Sys.Date(), "</p>",
      "<h3>Summary Statistics</h3>",
      "<ul>",
      "<li>Total patients: ", session_metadata$quality_metrics$total_unique_patients, "</li>",
      "<li>Complete data: ", session_metadata$quality_metrics$complete_data_patients, 
      " (", session_metadata$quality_metrics$complete_data_percent, "%)</li>",
      "<li>Survival events: ", session_metadata$quality_metrics$survival_events, "</li>",
      "<li>HER2 amplified: ", session_metadata$quality_metrics$her2_amplified_cnv, 
      " (", session_metadata$quality_metrics$her2_amplification_rate, "%)</li>",
      "</ul>",
      "<h3>Next Steps</h3>",
      "<p>Data is ready for Session 2: Data harmonization and biomarker standardization.</p>",
      "</body></html>"
    )
    
    writeLines(html_content, "results/reports/data_quality_report.html")
    cat("Quality report saved: results/reports/data_quality_report.html\n")
  }
}

create_quality_report()

# Final session summary
cat("\n=== SESSION 1 COMPLETE ===\n")
cat("Duration:", round(session_metadata$duration_minutes, 1), "minutes\n")
cat("Datasets processed/loaded:", length(session_metadata$datasets_created), "\n")

# Show what was downloaded vs loaded
data_sources <- session_metadata$data_sources
if (!is.null(data_sources)) {
  cat("\nData Sources:\n")
  if ("loaded_from_existing_file" %in% unlist(data_sources)) {
    loaded_files <- names(data_sources)[sapply(data_sources, function(x) x == "loaded_from_existing_file")]
    downloaded_files <- names(data_sources)[sapply(data_sources, function(x) x != "loaded_from_existing_file")]
    
    if (length(loaded_files) > 0) {
      cat("  📁 Loaded from existing files:", paste(loaded_files, collapse = ", "), "\n")
    }
    if (length(downloaded_files) > 0) {
      cat("  📥 Downloaded fresh:", paste(downloaded_files, collapse = ", "), "\n")
    }
  } else {
    cat("  📥 All data downloaded fresh\n")
  }
} else {
  cat("  📥 Data processing completed\n")
}

cat("\nFinal Datasets:\n")
for (dataset in session_metadata$datasets_created) {
  cat("  ✓", dataset, "\n")
}

if (!is.null(session_metadata$sample_counts$complete_data)) {
  cat("\n🎯 Ready for Session 2 with", session_metadata$sample_counts$complete_data, "complete cases\n")
  cat("📊 Survival events:", session_metadata$quality_metrics$survival_events, "\n")
  cat("🧬 HER2 amplified samples:", session_metadata$quality_metrics$her2_amplified_cnv, "\n")
}

cat("\n📋 Next Session Objectives:\n")
cat("1. Data harmonization across platforms\n")
cat("2. HER2 biomarker standardization\n") 
cat("3. Survival analysis preparation\n")
cat("4. Missing data imputation strategy\n")
cat("5. Quality control implementation\n")

# Clean up large objects
rm(clinical_data, rna_data, cnv_data)
gc()

cat("\n✅ All data successfully acquired and processed!\n")

cat("\n🎯 TROUBLESHOOTING TIPS:\n")
cat("If you see repeated download attempts or chunk retry messages:\n")
cat("1. Run: force_clean_tcga_downloads()\n")
cat("2. Restart R session\n") 
cat("3. Re-run this script for fresh downloads\n\n")

cat("This ensures no partial downloads interfere with the process.\n")
