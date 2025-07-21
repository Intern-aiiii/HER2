# HER2 Integration Project - Session 1: Data Acquisition
# Author: Research Team
# Date: January 2025
# Objective: Download and initial processing of TCGA-BRCA multi-platform data

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
  objective = "Data acquisition and initial processing",
  start_time = Sys.time(),
  datasets_created = character(),
  sample_counts = list(),
  quality_metrics = list(),
  next_session_inputs = character()
)

cat("=== HER2 Integration Project - Session 1 ===\n")
cat("Objective: TCGA-BRCA data acquisition\n")
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

# 1. CLINICAL DATA ACQUISITION
cat("\n1. ACQUIRING CLINICAL DATA\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

acquire_clinical_data <- function() {
  cat("Querying TCGA-BRCA clinical data...\n")
  
  # Query for clinical data
  clinical_query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
  )
  
  # Download clinical data
  GDCdownload(clinical_query, method = "api")
  
  # Prepare clinical data
  clinical_data <- GDCprepare_clinic(clinical_query, clinical.info = "patient")
  
  return(clinical_data)
}

clinical_data <- safe_execute(
  acquire_clinical_data(),
  "Clinical data acquisition"
)

# Process clinical data for HER2 status
process_clinical_data <- function(clinical_data) {
  if (is.null(clinical_data)) return(NULL)
  
  cat("Processing clinical data...\n")
  
  # Extract key variables and create processed dataset
  clinical_processed <- clinical_data %>%
    select(
      patient_id = submitter_id,
      age_at_diagnosis,
      gender,
      race,
      ethnicity,
      vital_status,
      days_to_death,
      days_to_last_follow_up,
      tumor_stage,
      pathologic_stage,
      histological_type,
      primary_diagnosis
    ) %>%
    # Convert age from days to years
    mutate(
      age_years = as.numeric(age_at_diagnosis) / 365.25,
      survival_time = ifelse(
        vital_status == "Dead",
        as.numeric(days_to_death),
        as.numeric(days_to_last_follow_up)
      ),
      event = ifelse(vital_status == "Dead", 1, 0)
    ) %>%
    filter(!is.na(patient_id))
  
  # Add placeholder for HER2 and other biomarker status (will extract in Session 2)
  clinical_processed$her2_ihc_status <- NA
  clinical_processed$her2_fish_status <- NA
  clinical_processed$er_status <- NA
  clinical_processed$pr_status <- NA
  clinical_processed$tumor_grade <- NA
  
  return(clinical_processed)
}

clinical_final <- safe_execute(
  process_clinical_data(clinical_data),
  "Clinical data processing"
)

# Save clinical data
if (!is.null(clinical_final)) {
  saveRDS(clinical_final, "data/processed/tcga_clinical.rds")
  session_metadata$datasets_created <- c(session_metadata$datasets_created, "tcga_clinical.rds")
  session_metadata$sample_counts$total_clinical <- nrow(clinical_final)
  cat("Clinical data saved:", nrow(clinical_final), "patients\n")
}

# 2. RNA EXPRESSION DATA ACQUISITION  
cat("\n2. ACQUIRING RNA EXPRESSION DATA\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

acquire_rna_data <- function() {
  # Query for RNA-Seq data (ERBB2 gene)
  rna_query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor")
  )
  
  # Download RNA data
  GDCdownload(rna_query, method = "api")
  
  # Prepare RNA data
  rna_data <- GDCprepare(rna_query)
  
  return(rna_data)
}

rna_data <- safe_execute(
  acquire_rna_data(),
  "RNA expression data acquisition"
)

# Process RNA data for ERBB2
process_rna_data <- function(rna_data) {
  if (is.null(rna_data)) return(NULL)
  
  # Extract ERBB2 expression
  gene_info <- rowData(rna_data)
  erbb2_idx <- which(gene_info$gene_name == "ERBB2")
  
  if (length(erbb2_idx) == 0) {
    cat("Warning: ERBB2 gene not found, searching for alternatives...\n")
    erbb2_idx <- grep("ERBB2|HER2", gene_info$gene_name, ignore.case = TRUE)
  }
  
  if (length(erbb2_idx) > 0) {
    # Extract expression data
    erbb2_expression <- assay(rna_data)[erbb2_idx[1], ]
    
    # Get sample information
    sample_info <- colData(rna_data)
    
    # Create RNA dataset
    rna_processed <- data.frame(
      patient_id = substr(sample_info$submitter_id, 1, 12),
      sample_id = sample_info$submitter_id,
      erbb2_counts = as.numeric(erbb2_expression),
      erbb2_log2 = log2(as.numeric(erbb2_expression) + 1),
      sample_type = sample_info$sample_type,
      stringsAsFactors = FALSE
    ) %>%
      filter(sample_type == "Primary Tumor") %>%
      filter(!is.na(erbb2_counts) & erbb2_counts > 0)
    
    return(rna_processed)
  } else {
    cat("Error: ERBB2 gene not found in expression data\n")
    return(NULL)
  }
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

# 3. COPY NUMBER DATA ACQUISITION
cat("\n3. ACQUIRING COPY NUMBER DATA\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

acquire_cnv_data <- function() {
  # Try different copy number data types and download methods
  cnv_options <- list(
    list(data.type = "Masked Copy Number Segment", desc = "Segment-level copy number"),
    list(data.type = "Copy Number Segment", desc = "Raw copy number segments"),
    list(data.type = "Gene Level Copy Number", desc = "Gene-level without scores suffix")
  )
  
  # Download methods to try (in order of preference)
  download_methods <- list(
    list(method = "client", files.per.chunk = 50, desc = "Client with 50 files/chunk"),
    list(method = "api", files.per.chunk = 25, desc = "API with 25 files/chunk"),
    list(method = "api", files.per.chunk = 10, desc = "API with 10 files/chunk")
  )
  
  successful_data <- NULL
  
  for(i in 1:length(cnv_options)) {
    if(!is.null(successful_data)) break  # Stop if we got data
    
    option <- cnv_options[[i]]
    cat("Testing copy number option", i, ":", option$desc, "...\n")
    
    tryCatch({
      test_query <- GDCquery(
        project = "TCGA-BRCA",
        data.category = "Copy Number Variation",
        data.type = option$data.type,
        sample.type = "Primary Tumor"
      )
      
      if(nrow(getResults(test_query)) > 0) {
        cat("✓ Found", nrow(getResults(test_query)), "files with", option$data.type, "\n")
        
        # Try different download methods for this data type
        for(j in 1:length(download_methods)) {
          if(!is.null(successful_data)) break  # Stop if we got data
          
          method <- download_methods[[j]]
          cat("  Trying download method:", method$desc, "...\n")
          
          tryCatch({
            # Download with current method
            GDCdownload(
              query = test_query,
              method = method$method,
              files.per.chunk = method$files.per.chunk
            )
            
            # Prepare the data
            cnv_data <- GDCprepare(test_query)
            
            # Store metadata about successful method
            attr(cnv_data, "data_type") <- option$data.type
            attr(cnv_data, "download_method") <- method$method
            attr(cnv_data, "files_per_chunk") <- method$files.per.chunk
            
            successful_data <- cnv_data
            cat("  ✓ Success with", method$desc, "\n")
            
          }, error = function(e) {
            cat("  ✗ Failed with", method$desc, ":", e$message, "\n")
          })
        }
      } else {
        cat("✗ No files found for", option$data.type, "\n")
      }
    }, error = function(e) {
      cat("✗ Query failed for option", i, ":", e$message, "\n")
    })
  }
  
  if(!is.null(successful_data)) {
    cat("✓ Copy number data acquired successfully!\n")
    return(successful_data)
  } else {
    stop("All copy number download methods failed")
  }
}

cnv_data <- safe_execute(
  acquire_cnv_data(),
  "Copy number data acquisition"
)

# Process CNV data for ERBB2
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
        n_segments = n(),
        .groups = "drop"
      ) %>%
      mutate(
        patient_id = substr(Sample, 1, 12),
        sample_id = Sample,
        sample_type = "Primary Tumor"
      ) %>%
      filter(!is.na(erbb2_cnv_score))
    
  } else {
    # For gene-level data, try to extract ERBB2 directly
    cat("Processing gene-level data for ERBB2...\n")
    
    # Check if this is a SummarizedExperiment with rowData
    if(class(cnv_data)[1] == "RangedSummarizedExperiment" || class(cnv_data)[1] == "SummarizedExperiment") {
      gene_info <- rowData(cnv_data)
      erbb2_idx <- which(gene_info$gene_name == "ERBB2" | gene_info$Hugo_Symbol == "ERBB2")
      
      if (length(erbb2_idx) == 0) {
        erbb2_idx <- grep("ERBB2|HER2", gene_info$gene_name, ignore.case = TRUE)
      }
      
      if (length(erbb2_idx) > 0) {
        erbb2_cnv <- assay(cnv_data)[erbb2_idx[1], ]
        sample_info <- colData(cnv_data)
        
        cnv_processed <- data.frame(
          patient_id = substr(sample_info$submitter_id, 1, 12),
          sample_id = sample_info$submitter_id,
          erbb2_cnv_score = as.numeric(erbb2_cnv),
          sample_type = sample_info$sample_type,
          stringsAsFactors = FALSE
        ) %>%
          filter(sample_type == "Primary Tumor") %>%
          filter(!is.na(erbb2_cnv_score))
      } else {
        cat("ERBB2 gene not found in gene-level data\n")
        return(NULL)
      }
    } else {
      cat("Unexpected data structure for gene-level data\n")
      return(NULL)
    }
  }
  
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
  cat("Copy number data saved:", nrow(cnv_final), "samples\n")
}

# 4. DATA QUALITY ASSESSMENT
cat("\n4. DATA QUALITY ASSESSMENT\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

generate_quality_report <- function() {
  # Sample overlap analysis
  if (!is.null(clinical_final) && !is.null(rna_final) && !is.null(cnv_final)) {
    
    clinical_ids <- clinical_final$patient_id
    rna_ids <- rna_final$patient_id
    cnv_ids <- cnv_final$patient_id
    
    # Find overlaps
    all_ids <- unique(c(clinical_ids, rna_ids, cnv_ids))
    
    overlap_matrix <- data.frame(
      patient_id = all_ids,
      has_clinical = all_ids %in% clinical_ids,
      has_rna = all_ids %in% rna_ids,
      has_cnv = all_ids %in% cnv_ids
    )
    
    # Complete data counts
    complete_data <- overlap_matrix %>%
      filter(has_clinical & has_rna & has_cnv)
    
    # Quality metrics
    quality_metrics <- list(
      total_unique_patients = length(all_ids),
      complete_data_patients = nrow(complete_data),
      clinical_only = sum(overlap_matrix$has_clinical & !overlap_matrix$has_rna & !overlap_matrix$has_cnv),
      rna_missing = sum(!overlap_matrix$has_rna),
      cnv_missing = sum(!overlap_matrix$has_cnv),
      complete_data_percent = round(nrow(complete_data) / length(all_ids) * 100, 1)
    )
    
    session_metadata$quality_metrics <<- quality_metrics
    session_metadata$sample_counts$complete_data <<- nrow(complete_data)
    
    # Print quality summary
    cat("Quality Assessment Summary:\n")
    cat("- Total unique patients:", quality_metrics$total_unique_patients, "\n")
    cat("- Complete data (all platforms):", quality_metrics$complete_data_patients, "\n")
    cat("- Complete data percentage:", quality_metrics$complete_data_percent, "%\n")
    cat("- RNA data missing:", quality_metrics$rna_missing, "patients\n")
    cat("- CNV data missing:", quality_metrics$cnv_missing, "patients\n")
    
    # Platform overlap details
    cat("\nPlatform Overlaps:\n")
    cat("- Clinical + RNA:", length(intersect(clinical_ids, rna_ids)), "patients\n")
    cat("- Clinical + CNV:", length(intersect(clinical_ids, cnv_ids)), "patients\n")
    cat("- RNA + CNV:", length(intersect(rna_ids, cnv_ids)), "patients\n")
    cat("- All three platforms:", nrow(complete_data), "patients\n")
    
    return(overlap_matrix)
  } else {
    cat("Warning: Not all datasets available for quality assessment\n")
    
    # Try to load missing datasets
    if(is.null(clinical_final) && file.exists("data/processed/tcga_clinical.rds")) {
      clinical_final <<- readRDS("data/processed/tcga_clinical.rds")
      cat("Loaded clinical data from file\n")
    }
    if(is.null(rna_final) && file.exists("data/processed/tcga_expression.rds")) {
      rna_final <<- readRDS("data/processed/tcga_expression.rds") 
      cat("Loaded RNA data from file\n")
    }
    if(is.null(cnv_final) && file.exists("data/processed/tcga_cnv.rds")) {
      cnv_final <<- readRDS("data/processed/tcga_cnv.rds")
      cat("Loaded CNV data from file\n")
    }
    
    # Try quality assessment again if we loaded data
    if(!is.null(clinical_final) && !is.null(rna_final) && !is.null(cnv_final)) {
      cat("Retrying quality assessment with loaded data...\n")
      return(generate_quality_report())
    }
    
    return(NULL)
  }
}

overlap_data <- safe_execute(
  generate_quality_report(),
  "Data quality assessment"
)

# 5. SAVE SESSION METADATA
cat("\n5. FINALIZING SESSION\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Update metadata with completion info
session_metadata$end_time <- Sys.time()
session_metadata$duration_minutes <- as.numeric(difftime(
  session_metadata$end_time, 
  session_metadata$start_time, 
  units = "mins"
))

# Set next session inputs
session_metadata$next_session_inputs <- session_metadata$datasets_created

# Save metadata
write_json(session_metadata, "metadata/session1_metadata.json", pretty = TRUE)

# Print session summary
cat("Session 1 Complete!\n")
cat("Duration:", round(session_metadata$duration_minutes, 1), "minutes\n")
cat("Datasets created:", length(session_metadata$datasets_created), "\n")
for (dataset in session_metadata$datasets_created) {
  cat("  -", dataset, "\n")
}

if (!is.null(session_metadata$sample_counts$complete_data)) {
  cat("Ready for Session 2 with", session_metadata$sample_counts$complete_data, "complete cases\n")
}

cat("\nNext steps for Session 2:\n")
cat("1. Data harmonization and standardization\n")
cat("2. HER2 IHC status extraction from clinical files\n") 
cat("3. Quality control procedures\n")
cat("4. Missing data imputation strategy\n")

# Clean up large objects to save memory
rm(clinical_data, rna_data, cnv_data)
gc()

cat("\n=== Session 1 Complete ===\n")