# HER2 Integration Project - Session 3: Bayesian Model Implementation
# Author: Research Team
# Date: January 2025
# Objective: Implement Bayesian latent variable model for multi-platform HER2 integration

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(ggplot2)
  library(MASS)
  library(boot)
  library(nloptr)
  library(Matrix)
  library(mvtnorm)
  library(optimx)
  library(numDeriv)
  library(gridExtra)
})

# Session metadata initialization
session_metadata <- list(
  session_id = 3,
  date = Sys.Date(),
  objective = "Bayesian model implementation",
  start_time = Sys.time(),
  model_components = c("latent_variable", "IHC_model", "RNA_model", "CNV_model"),
  convergence_achieved = FALSE,
  parameter_estimates = list(),
  bootstrap_iterations = 0,
  model_fit_metrics = list(),
  next_session_inputs = character()
)

cat("=== HER2 Integration Project - Session 3 ===\n")
cat("Objective: Bayesian latent variable model implementation\n")
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

# 1. LOAD HARMONIZED DATA FROM SESSION 2
cat("\n1. LOADING HARMONIZED DATA FROM SESSION 2\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

load_harmonized_data <- function() {
  # Check prerequisites
  required_files <- c("data/processed/harmonized_dataset.rds")
  
  missing_files <- required_files[!file.exists(required_files)]
  if(length(missing_files) > 0) {
    stop("Missing required files from Session 2: ", paste(missing_files, collapse = ", "))
  }
  
  # Load harmonized dataset
  harmonized_data <- readRDS("data/processed/harmonized_dataset.rds")
  
  # Filter for complete cases needed for modeling
  modeling_data <- harmonized_data %>%
    filter(
      !is.na(her2_ihc_status) & her2_ihc_status != "Unknown",
      !is.na(erbb2_log2),
      !is.na(erbb2_cnv_score),
      !is.na(age_years),
      qc_pass == TRUE
    ) %>%
    mutate(
      # Convert HER2 IHC to ordered factor
      her2_ihc_ordered = factor(her2_ihc_status, levels = c("0", "1+", "2+", "3+"), ordered = TRUE),
      her2_ihc_numeric = as.numeric(her2_ihc_ordered) - 1,  # 0, 1, 2, 3
      
      # Standardize continuous variables for modeling
      age_std = scale(age_years)[,1],
      rna_std = scale(erbb2_log2)[,1],
      cnv_std = scale(erbb2_cnv_score)[,1]
    ) %>%
    filter(!is.na(her2_ihc_ordered))
  
  cat("Modeling data prepared:\n")
  cat("- Total samples:", nrow(modeling_data), "\n")
  cat("- HER2 IHC distribution:\n")
  print(table(modeling_data$her2_ihc_status))
  
  # Check data completeness
  cat("- RNA expression range:", round(range(modeling_data$erbb2_log2), 2), "\n")
  cat("- Copy number range:", round(range(modeling_data$erbb2_cnv_score), 2), "\n")
  cat("- Age range:", round(range(modeling_data$age_years), 1), "\n")
  
  return(modeling_data)
}

modeling_data <- safe_execute(load_harmonized_data(), "Loading harmonized data")

if(is.null(modeling_data)) {
  stop("Cannot proceed without harmonized data. Please run Session 2 first.")
}

# 2. BAYESIAN LATENT VARIABLE MODEL SPECIFICATION
cat("\n2. BAYESIAN LATENT VARIABLE MODEL SPECIFICATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Model specification:
# ζᵢ ~ Beta(αᵢ, βᵢ) where logit(E[ζᵢ]) = β₀ + β₁(age) + β₂(stage) + ...
# IHC_score ~ OrderedLogistic(f(ζᵢ), cutpoints)
# RNA_expr ~ Normal(γ₀ + γ₁ζᵢ, σ²_RNA)  
# CNV_score ~ Normal(δ₀ + δ₁ζᵢ, σ²_CNV)

specify_bayesian_model <- function() {
  cat("Specifying Bayesian latent variable model components...\n")
  
  # Model parameters structure
  model_spec <- list(
    # Latent variable parameters
    latent_params = list(
      beta_intercept = 0,    # β₀: intercept for latent variable prior
      beta_age = 0,          # β₁: age effect on latent variable
      beta_stage = 0,        # β₂: stage effect on latent variable
      beta_histology = 0,    # β₃: histology effect on latent variable
      concentration = 1      # concentration parameter for Beta distribution
    ),
    
    # IHC measurement model parameters
    ihc_params = list(
      alpha = 0,             # α: effect of latent variable on IHC
      cutpoints = c(-2, 0, 2) # τ₁, τ₂, τ₃: cutpoints for ordered logistic
    ),
    
    # RNA measurement model parameters
    rna_params = list(
      gamma_0 = 0,           # γ₀: RNA intercept
      gamma_1 = 1,           # γ₁: effect of latent variable on RNA
      sigma_rna = 1          # σ_RNA: RNA measurement error
    ),
    
    # CNV measurement model parameters
    cnv_params = list(
      delta_0 = 0,           # δ₀: CNV intercept  
      delta_1 = 1,           # δ₁: effect of latent variable on CNV
      sigma_cnv = 1          # σ_CNV: CNV measurement error
    )
  )
  
  cat("Model specification complete:\n")
  cat("- Latent variable: Beta distribution with covariate-dependent parameters\n")
  cat("- IHC model: Ordered logistic regression\n")
  cat("- RNA model: Linear Gaussian regression\n")
  cat("- CNV model: Linear Gaussian regression\n")
  
  return(model_spec)
}

model_spec <- safe_execute(specify_bayesian_model(), "Model specification")

# 3. LIKELIHOOD FUNCTIONS
cat("\n3. IMPLEMENTING LIKELIHOOD FUNCTIONS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

implement_likelihood_functions <- function() {
  cat("Implementing likelihood functions for each component...\n")
  
  # Log-likelihood for IHC data (ordered logistic)
  ihc_log_likelihood <- function(zeta, ihc_scores, alpha, cutpoints) {
    # Linear predictor
    eta <- alpha * zeta
    
    # Calculate probabilities for each category
    n_cats <- length(cutpoints) + 1
    probs <- matrix(0, length(zeta), n_cats)
    
    # P(Y = 0) = logit⁻¹(τ₁ - η)
    probs[, 1] <- plogis(cutpoints[1] - eta)
    
    # P(Y = k) = logit⁻¹(τₖ₊₁ - η) - logit⁻¹(τₖ - η) for k = 1, ..., K-1
    for(k in 2:(n_cats-1)) {
      probs[, k] <- plogis(cutpoints[k] - eta) - plogis(cutpoints[k-1] - eta)
    }
    
    # P(Y = K) = 1 - logit⁻¹(τₖ - η)
    probs[, n_cats] <- 1 - plogis(cutpoints[n_cats-1] - eta)
    
    # Ensure probabilities are positive
    probs <- pmax(probs, 1e-10)
    
    # Calculate log-likelihood
    ll <- 0
    for(i in 1:length(ihc_scores)) {
      if(!is.na(ihc_scores[i]) && ihc_scores[i] >= 0 && ihc_scores[i] < n_cats) {
        ll <- ll + log(probs[i, ihc_scores[i] + 1])
      }
    }
    
    return(ll)
  }
  
  # Log-likelihood for RNA data (Gaussian)
  rna_log_likelihood <- function(zeta, rna_expr, gamma_0, gamma_1, sigma_rna) {
    mu <- gamma_0 + gamma_1 * zeta
    ll <- sum(dnorm(rna_expr, mean = mu, sd = sigma_rna, log = TRUE), na.rm = TRUE)
    return(ll)
  }
  
  # Log-likelihood for CNV data (Gaussian)
  cnv_log_likelihood <- function(zeta, cnv_scores, delta_0, delta_1, sigma_cnv) {
    mu <- delta_0 + delta_1 * zeta
    ll <- sum(dnorm(cnv_scores, mean = mu, sd = sigma_cnv, log = TRUE), na.rm = TRUE)
    return(ll)
  }
  
  # Prior log-density for latent variable
  latent_log_prior <- function(zeta, X, beta, concentration) {
    # Linear predictor for mean of Beta distribution
    eta <- X %*% beta
    mu <- plogis(eta)  # logit⁻¹(X β)
    
    # Beta distribution parameters
    alpha <- mu * concentration
    beta_param <- (1 - mu) * concentration
    
    # Ensure parameters are positive
    alpha <- pmax(alpha, 0.1)
    beta_param <- pmax(beta_param, 0.1)
    
    # Calculate log-density
    ll <- sum(dbeta(zeta, alpha, beta_param, log = TRUE), na.rm = TRUE)
    return(ll)
  }
  
  # Complete log-likelihood
  complete_log_likelihood <- function(params, data) {
    n <- nrow(data)
    
    # Extract parameters
    beta <- params[1:4]  # [intercept, age, stage, histology]
    concentration <- exp(params[5])  # ensure positive
    alpha <- params[6]
    cutpoints <- params[7:9]
    gamma_0 <- params[10]
    gamma_1 <- params[11] 
    sigma_rna <- exp(params[12])  # ensure positive
    delta_0 <- params[13]
    delta_1 <- params[14]
    sigma_cnv <- exp(params[15])  # ensure positive
    zeta <- params[16:(15+n)]  # latent variables
    
    # Ensure zeta is in [0,1]
    zeta <- pmax(pmin(zeta, 0.999), 0.001)
    
    # Design matrix for latent variable
    X <- cbind(1, data$age_std, as.numeric(data$stage_simplified == "III" | data$stage_simplified == "IV"),
               as.numeric(data$histology_simplified == "Lobular"))
    
    # Calculate log-likelihood components
    ll_latent <- latent_log_prior(zeta, X, beta, concentration)
    ll_ihc <- ihc_log_likelihood(zeta, data$her2_ihc_numeric, alpha, cutpoints)
    ll_rna <- rna_log_likelihood(zeta, data$rna_std, gamma_0, gamma_1, sigma_rna)
    ll_cnv <- cnv_log_likelihood(zeta, data$cnv_std, delta_0, delta_1, sigma_cnv)
    
    total_ll <- ll_latent + ll_ihc + ll_rna + ll_cnv
    
    return(total_ll)
  }
  
  cat("Likelihood functions implemented:\n")
  cat("- IHC: Ordered logistic regression likelihood\n")
  cat("- RNA: Gaussian likelihood\n") 
  cat("- CNV: Gaussian likelihood\n")
  cat("- Latent variable: Beta prior with covariate dependence\n")
  cat("- Complete: Joint likelihood for all components\n")
  
  return(list(
    ihc_ll = ihc_log_likelihood,
    rna_ll = rna_log_likelihood,
    cnv_ll = cnv_log_likelihood,
    latent_prior = latent_log_prior,
    complete_ll = complete_log_likelihood
  ))
}

likelihood_functions <- safe_execute(implement_likelihood_functions(), "Likelihood functions")

# 4. EM ALGORITHM IMPLEMENTATION
cat("\n4. IMPLEMENTING EM ALGORITHM\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

implement_em_algorithm <- function(data, likelihood_funcs, max_iter = 100, tol = 1e-6) {
  cat("Implementing EM algorithm for parameter estimation...\n")
  
  n <- nrow(data)
  
  # Initialize parameters
  params <- list(
    beta = c(0, 0, 0, 0),          # latent variable coefficients
    concentration = 2,              # Beta concentration parameter
    alpha = 1,                     # IHC slope
    cutpoints = c(-1, 0, 1),       # IHC cutpoints
    gamma_0 = 0,                   # RNA intercept
    gamma_1 = 1,                   # RNA slope
    sigma_rna = 1,                 # RNA error SD
    delta_0 = 0,                   # CNV intercept  
    delta_1 = 1,                   # CNV slope
    sigma_cnv = 1,                 # CNV error SD
    zeta = rep(0.5, n)             # latent variables
  )
  
  # Storage for convergence tracking
  ll_history <- numeric(max_iter)
  param_history <- list()
  
  cat("Starting EM iterations...\n")
  
  for(iter in 1:max_iter) {
    if(iter %% 10 == 0) cat("EM iteration", iter, "...\n")
    
    # E-step: Update latent variables given current parameters
    X <- cbind(1, data$age_std, 
               as.numeric(data$stage_simplified == "III" | data$stage_simplified == "IV"),
               as.numeric(data$histology_simplified == "Lobular"))
    
    # Update each latent variable individually
    for(i in 1:n) {
      # Define objective function for individual i
      objective <- function(z) {
        if(z <= 0 || z >= 1) return(-Inf)
        
        # Prior contribution
        eta <- X[i, ] %*% params$beta
        mu <- plogis(eta)
        alpha_beta <- mu * params$concentration
        beta_beta <- (1 - mu) * params$concentration
        ll_prior <- dbeta(z, alpha_beta, beta_beta, log = TRUE)
        
        # IHC likelihood
        eta_ihc <- params$alpha * z
        probs_ihc <- c(
          plogis(params$cutpoints[1] - eta_ihc),
          plogis(params$cutpoints[2] - eta_ihc) - plogis(params$cutpoints[1] - eta_ihc),
          plogis(params$cutpoints[3] - eta_ihc) - plogis(params$cutpoints[2] - eta_ihc),
          1 - plogis(params$cutpoints[3] - eta_ihc)
        )
        probs_ihc <- pmax(probs_ihc, 1e-10)
        
        ll_ihc <- 0
        if(!is.na(data$her2_ihc_numeric[i])) {
          ll_ihc <- log(probs_ihc[data$her2_ihc_numeric[i] + 1])
        }
        
        # RNA likelihood
        mu_rna <- params$gamma_0 + params$gamma_1 * z
        ll_rna <- dnorm(data$rna_std[i], mu_rna, params$sigma_rna, log = TRUE)
        
        # CNV likelihood
        mu_cnv <- params$delta_0 + params$delta_1 * z
        ll_cnv <- dnorm(data$cnv_std[i], mu_cnv, params$sigma_cnv, log = TRUE)
        
        return(ll_prior + ll_ihc + ll_rna + ll_cnv)
      }
      
      # Optimize to find best zeta[i]
      opt_result <- tryCatch({
        optimize(objective, interval = c(0.001, 0.999), maximum = TRUE)
      }, error = function(e) {
        list(maximum = params$zeta[i])  # Keep current value if optimization fails
      })
      
      params$zeta[i] <- opt_result$maximum
    }
    
    # M-step: Update parameters given current latent variables
    # Update latent variable coefficients (beta)
    beta_objective <- function(beta_vec) {
      eta <- X %*% beta_vec
      mu <- plogis(eta)
      alpha_beta <- mu * params$concentration
      beta_beta <- (1 - mu) * params$concentration
      alpha_beta <- pmax(alpha_beta, 0.1)
      beta_beta <- pmax(beta_beta, 0.1)
      ll <- sum(dbeta(params$zeta, alpha_beta, beta_beta, log = TRUE))
      return(-ll)  # negative for minimization
    }
    
    beta_opt <- tryCatch({
      optim(params$beta, beta_objective, method = "BFGS")
    }, error = function(e) {
      list(par = params$beta)  # Keep current values if optimization fails
    })
    params$beta <- beta_opt$par
    
    # Update RNA parameters
    rna_lm <- lm(data$rna_std ~ params$zeta)
    params$gamma_0 <- coef(rna_lm)[1]
    params$gamma_1 <- coef(rna_lm)[2]
    params$sigma_rna <- summary(rna_lm)$sigma
    
    # Update CNV parameters
    cnv_lm <- lm(data$cnv_std ~ params$zeta)
    params$delta_0 <- coef(cnv_lm)[1]
    params$delta_1 <- coef(cnv_lm)[2]
    params$sigma_cnv <- summary(cnv_lm)$sigma
    
    # Calculate current log-likelihood
    param_vec <- c(params$beta, log(params$concentration), params$alpha, params$cutpoints,
                   params$gamma_0, params$gamma_1, log(params$sigma_rna),
                   params$delta_0, params$delta_1, log(params$sigma_cnv), params$zeta)
    
    ll_current <- tryCatch({
      likelihood_funcs$complete_ll(param_vec, data)
    }, error = function(e) {
      -Inf
    })
    
    ll_history[iter] <- ll_current
    param_history[[iter]] <- params
    
    # Check convergence
    if(iter > 1) {
      if(abs(ll_history[iter] - ll_history[iter-1]) < tol) {
        cat("EM algorithm converged at iteration", iter, "\n")
        session_metadata$convergence_achieved <<- TRUE
        break
      }
    }
  }
  
  cat("EM algorithm completed:\n")
  cat("- Final log-likelihood:", round(ll_current, 2), "\n")
  cat("- Convergence achieved:", session_metadata$convergence_achieved, "\n")
  
  # Store parameter estimates
  session_metadata$parameter_estimates <<- list(
    beta = params$beta,
    concentration = params$concentration,
    alpha = params$alpha,
    cutpoints = params$cutpoints,
    gamma_0 = params$gamma_0,
    gamma_1 = params$gamma_1,
    sigma_rna = params$sigma_rna,
    delta_0 = params$delta_0,
    delta_1 = params$delta_1,
    sigma_cnv = params$sigma_cnv,
    final_ll = ll_current
  )
  
  return(list(
    final_params = params,
    ll_history = ll_history[1:iter],
    param_history = param_history[1:iter],
    converged = session_metadata$convergence_achieved,
    iterations = iter
  ))
}

em_results <- safe_execute(
  implement_em_algorithm(modeling_data, likelihood_functions),
  "EM algorithm estimation"
)

# 5. BOOTSTRAP UNCERTAINTY QUANTIFICATION
cat("\n5. BOOTSTRAP UNCERTAINTY QUANTIFICATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

perform_bootstrap_uncertainty <- function(data, em_results, n_bootstrap = 100) {
  cat("Performing bootstrap for uncertainty quantification...\n")
  
  n <- nrow(data)
  bootstrap_params <- list()
  bootstrap_success <- 0
  
  cat("Running", n_bootstrap, "bootstrap iterations...\n")
  
  for(b in 1:n_bootstrap) {
    if(b %% 20 == 0) cat("Bootstrap iteration", b, "...\n")
    
    tryCatch({
      # Bootstrap sample
      boot_indices <- sample(1:n, n, replace = TRUE)
      boot_data <- data[boot_indices, ]
      
      # Run EM on bootstrap sample (fewer iterations for speed)
      boot_em <- implement_em_algorithm(boot_data, likelihood_functions, 
                                        max_iter = 30, tol = 1e-4)
      
      if(boot_em$converged) {
        bootstrap_params[[bootstrap_success + 1]] <- boot_em$final_params
        bootstrap_success <- bootstrap_success + 1
      }
      
    }, error = function(e) {
      # Skip failed bootstrap samples
    })
  }
  
  cat("Bootstrap completed:\n")
  cat("- Successful iterations:", bootstrap_success, "out of", n_bootstrap, "\n")
  
  session_metadata$bootstrap_iterations <<- bootstrap_success
  
  # Calculate confidence intervals
  if(bootstrap_success > 10) {
    # Extract parameter estimates from successful bootstrap samples
    beta_boots <- sapply(bootstrap_params, function(x) x$beta)
    gamma_boots <- sapply(bootstrap_params, function(x) c(x$gamma_0, x$gamma_1))
    delta_boots <- sapply(bootstrap_params, function(x) c(x$delta_0, x$delta_1))
    
    # Calculate 95% confidence intervals
    confidence_intervals <- list(
      beta = apply(beta_boots, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE),
      gamma = apply(gamma_boots, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE),
      delta = apply(delta_boots, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    )
    
    cat("- 95% confidence intervals calculated\n")
    
    return(list(
      bootstrap_params = bootstrap_params,
      confidence_intervals = confidence_intervals,
      n_successful = bootstrap_success
    ))
  } else {
    cat("- Insufficient successful bootstrap samples for confidence intervals\n")
    return(NULL)
  }
}

bootstrap_results <- safe_execute(
  perform_bootstrap_uncertainty(modeling_data, em_results, n_bootstrap = 50),
  "Bootstrap uncertainty quantification"
)

# 6. MODEL FIT DIAGNOSTICS
cat("\n6. MODEL FIT DIAGNOSTICS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

perform_model_diagnostics <- function(data, em_results, bootstrap_results) {
  cat("Performing model fit diagnostics...\n")
  
  final_params <- em_results$final_params
  
  # Calculate AIC and BIC
  n_params <- length(final_params$beta) + 1 + 1 + length(final_params$cutpoints) + 
    2 + 1 + 2 + 1  # Total number of parameters
  n_obs <- nrow(data)
  
  final_ll <- tail(em_results$ll_history, 1)
  aic <- -2 * final_ll + 2 * n_params
  bic <- -2 * final_ll + log(n_obs) * n_params
  
  # Store model fit metrics
  session_metadata$model_fit_metrics <<- list(
    log_likelihood = final_ll,
    AIC = aic,
    BIC = bic,
    n_parameters = n_params,
    n_observations = n_obs
  )
  
  cat("Model fit diagnostics:\n")
  cat("- Log-likelihood:", round(final_ll, 2), "\n")
  cat("- AIC:", round(aic, 2), "\n")
  cat("- BIC:", round(bic, 2), "\n")
  cat("- Number of parameters:", n_params, "\n")
  
  # Generate convergence plots
  pdf("results/figures/convergence_diagnostics.pdf", width = 12, height = 8)
  
  # Plot 1: Log-likelihood convergence
  par(mfrow = c(2, 2))
  
  plot(em_results$ll_history, type = "l", main = "EM Algorithm Convergence",
       xlab = "Iteration", ylab = "Log-Likelihood", lwd = 2, col = "blue")
  grid()
  
  # Plot 2: Latent variable distribution
  hist(final_params$zeta, breaks = 20, main = "Latent Variable Distribution",
       xlab = "Latent HER2 Biology (ζ)", ylab = "Frequency", 
       col = "lightblue", border = "black")
  
  # Plot 3: Latent variable vs RNA expression
  plot(final_params$zeta, data$rna_std, main = "Latent Variable vs RNA Expression",
       xlab = "Latent HER2 Biology (ζ)", ylab = "Standardized RNA Expression",
       pch = 16, col = alpha("red", 0.6))
  abline(lm(data$rna_std ~ final_params$zeta), col = "blue", lwd = 2)
  
  # Plot 4: Latent variable vs Copy Number
  plot(final_params$zeta, data$cnv_std, main = "Latent Variable vs Copy Number",
       xlab = "Latent HER2 Biology (ζ)", ylab = "Standardized Copy Number",
       pch = 16, col = alpha("green", 0.6))
  abline(lm(data$cnv_std ~ final_params$zeta), col = "blue", lwd = 2)
  
  dev.off()
  
  cat("- Convergence diagnostics saved to convergence_diagnostics.pdf\n")
  
  return(list(
    aic = aic,
    bic = bic,
    convergence_plot = "convergence_diagnostics.pdf"
  ))
}

diagnostics_results <- safe_execute(
  perform_model_diagnostics(modeling_data, em_results, bootstrap_results),
  "Model fit diagnostics"
)

# 7. SAVE FITTED MODEL
cat("\n7. SAVING FITTED MODEL\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

save_fitted_model <- function(data, em_results, bootstrap_results, diagnostics) {
  cat("Saving fitted Bayesian model...\n")
  
  # Compile complete model object
  fitted_model <- list(
    model_type = "Bayesian Latent Variable Model",
    data_summary = list(
      n_observations = nrow(data),
      n_complete_cases = nrow(data),
      her2_ihc_distribution = table(data$her2_ihc_status),
      platforms = c("IHC", "RNA", "CNV")
    ),
    parameters = em_results$final_params,
    parameter_estimates = session_metadata$parameter_estimates,
    convergence = list(
      converged = em_results$converged,
      iterations = em_results$iterations,
      ll_history = em_results$ll_history
    ),
    uncertainty = bootstrap_results,
    model_fit = session_metadata$model_fit_metrics,
    diagnostics = diagnostics,
    session_info = session_metadata
  )
  
  # Save main fitted model
  saveRDS(fitted_model, "data/results/fitted_model.rds")
  
  # Save parameter estimates separately for easy access
  saveRDS(em_results$final_params, "data/results/parameter_estimates.rds")
  
  # Save latent variable scores for downstream analysis
  latent_scores <- data.frame(
    patient_id = data$patient_id,
    latent_her2_score = em_results$final_params$zeta,
    her2_ihc_status = data$her2_ihc_status,
    rna_expression = data$erbb2_log2,
    cnv_score = data$erbb2_cnv_score
  )
  saveRDS(latent_scores, "data/results/latent_her2_scores.rds")
  
  session_metadata$next_session_inputs <<- c(
    "fitted_model.rds",
    "parameter_estimates.rds", 
    "latent_her2_scores.rds",
    "harmonized_dataset.rds"
  )
  
  cat("Fitted model saved:\n")
  cat("- Main model: fitted_model.rds\n")
  cat("- Parameters: parameter_estimates.rds\n") 
  cat("- Latent scores: latent_her2_scores.rds\n")
  cat("- Ready for Session 4 validation\n")
  
  return(fitted_model)
}

fitted_model <- safe_execute(
  save_fitted_model(modeling_data, em_results, bootstrap_results, diagnostics_results),
  "Saving fitted model"
)

# 8. FINALIZE SESSION 3
cat("\n8. FINALIZING SESSION 3\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Update metadata with completion info
session_metadata$end_time <- Sys.time()
session_metadata$duration_minutes <- as.numeric(difftime(
  session_metadata$end_time, 
  session_metadata$start_time, 
  units = "mins"
))

# Save session metadata
write_json(session_metadata, "metadata/session3_metadata.json", pretty = TRUE)

# Print session summary
cat("Session 3 Complete!\n")
cat("Duration:", round(session_metadata$duration_minutes, 1), "minutes\n")
cat("Model components implemented:", length(session_metadata$model_components), "\n")
for(component in session_metadata$model_components) {
  cat("  -", component, "\n")
}

cat("Convergence achieved:", session_metadata$convergence_achieved, "\n")
cat("Bootstrap iterations:", session_metadata$bootstrap_iterations, "\n")

if(!is.null(session_metadata$model_fit_metrics$AIC)) {
  cat("Model fit - AIC:", round(session_metadata$model_fit_metrics$AIC, 2), "\n")
  cat("Model fit - BIC:", round(session_metadata$model_fit_metrics$BIC, 2), "\n")
}

cat("\nNext steps for Session 4:\n")
cat("1. Cross-validation performance assessment\n")
cat("2. Model comparison vs individual platforms\n") 
cat("3. Uncertainty calibration assessment\n")
cat("4. Predictive accuracy evaluation\n")

# Clean up large objects to save memory
rm(modeling_data, likelihood_functions)
gc()

cat("\n=== Session 3 Complete ===\n")