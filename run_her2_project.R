# HER2 Integration Project - Master Workflow Script
# Run this script to execute the complete analysis pipeline

cat("=== HER2 Integration Project - Master Workflow ===\n")
cat("Starting complete analysis pipeline...\n\n")

# Session control - set which sessions to run
run_session_1 <- TRUE   # Data acquisition
run_session_2 <- TRUE  # Data harmonization  
run_session_3 <- TRUE  # Bayesian modeling
run_session_4 <- FALSE  # Model validation
run_session_5 <- FALSE  # Clinical analysis
run_session_6 <- FALSE  # Publication figures

generate_reports <- TRUE  # Generate quality reports after each session

# SESSION 1: Data Acquisition
if(run_session_1) {
  cat("=== SESSION 1: Data Acquisition ===\n")
  start_time <- Sys.time()
  
  # Check if data already exists
  if(file.exists("data/processed/tcga_clinical.rds") && 
     file.exists("data/processed/tcga_expression.rds")) {
    cat("Data files already exist. Skipping download.\n")
    cat("To re-download, delete files in data/processed/ first.\n")
  } else {
    cat("Running data acquisition...\n")
    source("01_data_acquisition.R")
  }
  
  # Generate quality report
  if(generate_reports) {
    cat("\nGenerating Session 1 quality report...\n")
    if(file.exists("data_quality_report.Rmd")) {
      rmarkdown::render("data_quality_report.Rmd", 
                        output_file = "results/reports/session1_quality_report.html")
      cat("✓ Quality report saved to: results/reports/session1_quality_report.html\n")
    } else {
      cat("⚠ data_quality_report.Rmd not found\n")
    }
  }
  
  duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
  cat("Session 1 completed in", duration, "minutes\n\n")
}

# SESSION 2: Data Harmonization
if(run_session_2) {
  cat("=== SESSION 2: Data Harmonization ===\n")
  
  # Check prerequisites
  if(!file.exists("data/processed/tcga_clinical.rds")) {
    stop("Session 1 data not found. Run Session 1 first.")
  }
  
  start_time <- Sys.time()
  source("02_data_harmonization.R")
  duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
  cat("Session 2 completed in", duration, "minutes\n\n")
}

# SESSION 3: Bayesian Model Implementation
if(run_session_3) {
  cat("=== SESSION 3: Bayesian Model Implementation ===\n")
  
  if(!file.exists("data/processed/harmonized_dataset.rds")) {
    stop("Session 2 harmonized data not found. Run Session 2 first.")
  }
  
  start_time <- Sys.time()
  source("03_bayesian_model.R")
  duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
  cat("Session 3 completed in", duration, "minutes\n\n")
}

# SESSION 4: Model Validation
if(run_session_4) {
  cat("=== SESSION 4: Model Validation ===\n")
  
  if(!file.exists("data/results/fitted_model.rds")) {
    stop("Session 3 model not found. Run Session 3 first.")
  }
  
  start_time <- Sys.time()
  source("04_model_validation.R")
  duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
  cat("Session 4 completed in", duration, "minutes\n\n")
}

# SESSION 5: Clinical Analysis
if(run_session_5) {
  cat("=== SESSION 5: Clinical Analysis ===\n")
  
  if(!file.exists("data/results/validation_results.rds")) {
    stop("Session 4 validation not found. Run Session 4 first.")
  }
  
  start_time <- Sys.time()
  source("05_clinical_analysis.R")
  duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
  cat("Session 5 completed in", duration, "minutes\n\n")
}

# SESSION 6: Publication Figures
if(run_session_6) {
  cat("=== SESSION 6: Publication Figures ===\n")
  
  if(!file.exists("data/results/survival_results.rds")) {
    stop("Session 5 results not found. Run Session 5 first.")
  }
  
  start_time <- Sys.time()
  source("06_publication_figures.R")
  duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
  cat("Session 6 completed in", duration, "minutes\n\n")
}

# Final summary
cat("=== PROJECT COMPLETION SUMMARY ===\n")
completed_sessions <- sum(c(run_session_1, run_session_2, run_session_3, 
                            run_session_4, run_session_5, run_session_6))
cat("Completed sessions:", completed_sessions, "of 6\n")

if(completed_sessions == 6) {
  cat("🎉 Complete HER2 Integration project finished!\n")
  cat("📊 Check results/figures/ for publication figures\n")
  cat("📋 Check results/tables/ for manuscript tables\n")
} else {
  cat("📝 To continue, set run_session_", completed_sessions + 1, " <- TRUE and re-run\n")
}

cat("=== Master Workflow Complete ===\n")