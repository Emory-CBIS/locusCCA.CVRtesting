#' Helper Functions for Low-Rank CCA and Score Testing
#'
#' These functions include matrix transformations and the SCAD penalty function.

#' Transform Vector into a Symmetric Matrix
#'
#' @param x A numeric vector representing the upper triangular part of a matrix.
#' @param V An integer specifying the dimension of the square matrix.
#' @param d Logical, if `TRUE`, includes diagonal elements.
#' @param d_vec A numeric value for diagonal scaling (default = 1).
#' @return A symmetric matrix of dimension (V x V).
#' @export
Ltrinv <- function(x, V, d = TRUE, d_vec = 1) {
  Y = matrix(0, ncol = V, nrow = V)
  if (d) {
    Y[upper.tri(Y, FALSE)] = x
    return(Y + t(Y) + diag(d_vec * rep(1, V)))
  } else {
    Y[upper.tri(Y, FALSE)] = x
    return(Y + t(Y))
  }
}

#' Extract Upper Triangular Part of a Matrix
#'
#' @param X A square matrix.
#' @param d Logical, if `TRUE`, includes diagonal elements.
#' @return A numeric vector of the upper triangular part of `X`.
#' @export
Ltrans <- function(X, d = FALSE) {
  X[upper.tri(X, d)]
}

#' Smoothly Clipped Absolute Deviation (SCAD) Penalty Function
#'
#' @param yvpen A numeric vector of values to be penalized.
#' @param lambda_ch A numeric value for the SCAD penalty tuning parameter (default = 0.01).
#' @param gamma A numeric value for SCAD adjustment (default = 3). Must be greater than 2.
#' @return A numeric vector with SCAD-penalized values.
#' @export
SCAD_func <- function(yvpen, lambda_ch = 0.01, gamma = 3) {
  if (gamma <= 2) {
    gamma = 2.01
    warning("Gamma needs to be > 2!")
  }

  ynew = sign(yvpen) * (abs(yvpen) - lambda_ch) * (abs(yvpen) >= lambda_ch) * (abs(yvpen) <= 2 * lambda_ch) +
    yvpen * (abs(yvpen) > gamma * lambda_ch) +
    ((gamma - 1) * yvpen - sign(yvpen) * gamma * lambda_ch) / (gamma - 2) *
    (abs(yvpen) <= gamma * lambda_ch) * (abs(yvpen) > 2 * lambda_ch)

  if (sd(ynew) < 1e-7) {
    warning("Parameters are not correctly specified!")
    return(ynew)
  }

  return(ynew)
}
