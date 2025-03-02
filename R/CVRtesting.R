#' CVR Score Testing for Canonical Direction Weights
#'
#' Conducts a statistical test to assess the significance of extracted canonical components in characterizing some responses related to mental health.
#'
#' @param U A numeric matrix of canonical weights (p x m).
#' @param X A numeric matrix of predictor variables (n x p).
#' @param z A numeric response vector (n x 1).
#' @param lambda1 A numeric value for Lasso penalty in regression. If `NULL`, it is determined using cross-validation.
#' @param lambda2 A numeric value for thresholding. If `NULL`, it is  determined by our procedure.
#' @return A numeric vector of test T-statistics for each canonical component.
#' @importFrom glmnet glmnet cv.glmnet
#' @export
CVR_testing <- function(U, X, z, lambda1 = NULL, lambda2 = NULL) {
  q = dim(U)[2]
  n = dim(X)[1]

  # Compute factor matrix F_hat
  F_hat = X %*% U %*% diag(1/apply(X %*% U, 2, sd))
  R_hat = X - F_hat %*% solve(t(F_hat) %*% F_hat) %*% t(F_hat) %*% X

  # Estimate error using Lasso
  beta_est1 = cv.glmnet(x = X, y = z, intercept = TRUE, nfolds = 3, family = 'gaussian', alpha = 1, standardize = TRUE)
  beta_est = glmnet(x = X, y = z, intercept = TRUE, family = 'gaussian', alpha = 1, lambda = beta_est1$lambda.min, standardize = TRUE)
  fitted = predict(beta_est, newx = X, type = "response")

  # Initialize test statistics vector
  T_stats = c()

  for (j in 1:q) {

    sigma_sqr = mean((z - fitted)^2)
    epsilon = quantile(abs(U[, j]), 0.05)
    f_hat_j = F_hat[, j]
    index = abs(U[, j]) > epsilon

    # Determine lambda1 using cross-validation if not provided
    if (is.null(lambda1)) {
      glmnet_j_cv = cv.glmnet(x = cbind(f_hat_j, X[, !index]), y = z, intercept = TRUE, family = 'gaussian',
                              nfolds = 3, alpha = 1, standardize = TRUE, penalty.factor = c(0, rep(1, sum(!index))))
      lambda1 = glmnet_j_cv$lambda.1se
    }

    # Fit Lasso regression using lambda1
    glmnet_j = glmnet(x = cbind(f_hat_j, X[, !index]), y = z, intercept = TRUE, family = 'gaussian',
                      alpha = 1, lambda = lambda1, standardize = TRUE, penalty.factor = c(0, rep(1, sum(!index))), relax = TRUE)

    # Determine lambda2 using a thresholding rule if not provided
    if (is.null(lambda2)) {
      lambda2 = sqrt(log(length(glmnet_j$beta)) / (2 * n))
    }

    temp = 0.01
    while (TRUE) {
      w_null = glmnet(x = X[, !index], intercept = TRUE, y = f_hat_j, family = 'gaussian', alpha = 1,
                      standardize = TRUE, lambda = temp, relax = TRUE)
      condition = max(1 / n * abs(t(X[, !index]) %*% (f_hat_j - predict(w_null, newx = X[, !index], type = 'response'))))
      if (condition > lambda2) {
        temp = temp / 1.5
      } else {
        break
      }
    }

    # Compute test statistic
    S = 1 / (n * sigma_sqr) * t(z - predict(glmnet_j, newx = cbind(0, X[, !index]), type = 'response')) %*%
      (f_hat_j - predict(w_null, newx = X[, !index], type = 'response'))

    I = 1 / (n * sigma_sqr) * (sum(f_hat_j^2) - sum(f_hat_j * predict(w_null, newx = X[, !index], type = 'response')))
    #if (I <= 0)
    T_stat = sqrt(n) * 1 / sqrt(abs(I)) * S
    T_stats = c(T_stats, T_stat)
  }

  return(T_stats)
}
