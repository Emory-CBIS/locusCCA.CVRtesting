#' Compute Overall BIC for Canonical Directions
#'
#' This function calculates the Bayesian Information Criterion (BIC) for the entire set of
#' canonical components, based on the Gaussian likelihood between `XU[,j]` and `YV[,j]`.
#' The number of nonzero entries is determined using a threshold of 10% of the absolute quantiles.
#'
#' @param X A numeric matrix of predictor variables (n x p).
#' @param Y A numeric matrix of response variables (n x m).
#' @param U A numeric matrix of canonical weights for X (p x m).
#' @param V A numeric matrix of canonical weights for Y (q x m).
#' @return A single numeric value representing the overall BIC across all components.
#' @export
BIC_calc <- function(X, Y, U, V) {
  n <- nrow(X)  # Number of observations
  q <- ncol(U)  # Number of canonical components

  total_BIC <- 0  # Initialize overall BIC value

  for (j in 1:q) {
    # Project X and Y onto U and V respectively
    XU_j <- X %*% U[, j]
    YV_j <- Y %*% V[, j]

    # Compute residual variance using Gaussian likelihood
    model <- lm(YV_j ~ XU_j)
    sigma2 <- sum(model$residuals^2) / n  # Estimate variance

    # Determine nonzero entries based on 10% quantile threshold
    threshold <- quantile(abs(U[, j]), 0.10)  # Compute 10% quantile of absolute values
    nonzero_count <- sum(abs(U[, j]) > threshold)  # Count nonzero entries

    # BIC calculation: log-likelihood with penalty for nonzero parameters
    log_likelihood <- -0.5 * n * log(2 * pi * sigma2) - sum(model$residuals^2) / (2 * sigma2)
    penalty <- 0.5 * log(n) * nonzero_count  # BIC penalty term

    # Sum up BIC values for all components
    total_BIC <- total_BIC + (-2 * log_likelihood + penalty)
  }

  return(total_BIC)  # Return single summed BIC value
}
