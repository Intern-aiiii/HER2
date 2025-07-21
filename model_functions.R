
# Fixed ordered logistic likelihood
ordered_logistic_loglik <- function(y, eta, cutpoints) {
  n <- length(y)
  loglik <- 0
  
  # Ensure cutpoints are ordered
  cutpoints <- sort(cutpoints)
  
  for (i in 1:n) {
    if (is.na(y[i])) {
      next  # Skip missing values
    }
    
    if (y[i] == 0) {
      prob <- plogis(cutpoints[1] - eta[i])
    } else if (y[i] == 1) {
      prob <- plogis(cutpoints[2] - eta[i]) - plogis(cutpoints[1] - eta[i])
    } else if (y[i] == 2) {
      prob <- plogis(cutpoints[3] - eta[i]) - plogis(cutpoints[2] - eta[i])
    } else if (y[i] == 3) {
      prob <- 1 - plogis(cutpoints[3] - eta[i])
    } else {
      # Handle unexpected values
      prob <- 1e-10
    }
    
    # Add small constant to avoid log(0)
    prob <- pmax(prob, 1e-10)
    loglik <- loglik + log(prob)
  }
  
  return(loglik)
}

# Beta distribution parameters from mean and precision
beta_params_from_mean_precision <- function(mu, phi) {
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  return(list(alpha = alpha, beta = beta))
}

# Additional model functions would go here...

