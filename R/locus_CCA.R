#' Locus Canonical Correlation Analysis (locus_CCA)
#'
#' This function performs locus_CCA with sparse regularization.
#'
#' @param X A matrix of predictor variables (n x p).
#' @param Y A matrix of response variables (n x m).
#' @param node An integer indicating the node size.
#' @param m An integer specifying the number of canonical components.
#' @param rho A numeric value regulating the penalty term applied to canonical weights on brain connectivity.
#' @param gamma A numeric value (default = 2.1) used in SCAD penalty. Must be > 2.
#' @param penalt A string specifying the penalty type. Options:
#'        - "Hardthreshold": Applies hard thresholding.
#'        - "SCAD": Uses SCAD regularization.
#'        - "L1": Applies L1 regularization.
#'        - Other: No penalty.
#' @param proportion A numeric value (default = 0.9) indicating the threshold proportion for the connectivity traits, determining the ranks of our low rank structure.
#' @param silent Logical. If TRUE, suppresses output messages (default = FALSE).
#' @param tol A numeric tolerance level for convergence (default = 1e-3).
#' @return A list containing:
#' \item{U}{Canonical variable matrix for `X` (p x m).}
#' \item{V}{Canonical variable matrix for `Y` (q x m).}
#' \item{CC}{Canonical correlations between the each corresponding projection, m by m matrices.}
#' \item{R}{Vector of selected rank for each canonical component.}
#' @importFrom PMA CCA
#' @export
Locus_CCA <- function(X, Y, node, m, rho, gamma = 2.1,
                         penalt = 'L1', proportion = 0.9
                         , silent = FALSE, tol = 1e-3) {
  initial_sparse = 0.5
  imput_method = "No"  # Default imputation method
  Y_c = scale(Y)  # Standardize Y
  At = eigen(t(Y_c) %*% Y_c)
  transT = At$vectors[,1:m] %*% diag(1/sqrt(At$values[1:m]))
  Y = Y_c %*% transT

  # Initial CCA
  cca = CCA(X, Y, typex = 'standard', typez = 'standard', K = m,
            standardize = FALSE, penaltyx = initial_sparse, trace = FALSE)

  V = far::orthonormalization(cca$v)
  V_0 = far::orthonormalization(cca$v)
  U = cca$u
  U_0 = cca$u
  V_pre = 0
  U_pre = U
  R = c()
  Iter = 0
  sd1 = apply(U, 2, sd)

  # Low-rank approximation
  for (j in 1:m) {
    eigenSl = eigen(Ltrinv(U[, j], V = node, F))
    orderEigen = order(abs(eigenSl$values), decreasing = TRUE)
    for (r in 2:length(eigenSl$values)) {
      eigenset = orderEigen[1:r]
      imgeRL = eigenSl$vectors[, eigenset] %*% diag(eigenSl$values[eigenset]) %*% t(eigenSl$vectors[, eigenset])
      if ((cor(Ltrans(imgeRL, FALSE), U[, j]) >= proportion) || (r > 10)) {
        R = c(R, r)
        break
      }
    }
    U[, j] = Ltrans(imgeRL)
  }

  U_lr = U
  J_l = list()
  eigen_l = list()

  # Iterative CCA refinement
  while (sum((U - U_pre)^2) / sum(U_pre^2) + sum((V - V_pre)^2) / sum(V_pre^2) > tol) {
    if (Iter > 100) break
    Iter = Iter + 1
    V_pre = V
    U_pre = U

    cca3 = cc(X %*% cbind(U), Y)
    V1 = far::orthonormalization(cca3$ycoef)

    # Align V with previous iterations
    V_temp = V1
    for (j in 1:m) {
      cor_j = which.max(abs(cor(V1, V[, j])))
      V_tar = if (j == m) V1 else V1[, cor_j]
      if (j <= m - 1) V1 = V1[, -cor_j]
      V_temp[, j] = if (cor(V_tar, V[, j]) > 0) V_tar else -V_tar
    }
    V = 0.7 * V + 0.3 * V_temp
    V = far::orthonormalization(V)

    cca2 = CCA(X, Y %*% V, typex = 'standard', typez = 'standard', K = m,
               standardize = TRUE, penaltyx = initial_sparse, trace = FALSE)
    U1 = cca2$u

    # Align U with previous iterations
    U_temp = U1
    for (j in 1:m) {
      cor_j = which.max(abs(cor(U1, U[, j])))
      U_tar = if (j == m) U1 else U1[, cor_j]
      if (j <= m - 1) U1 = U1[, -cor_j]
      U_temp[, j] = if (cor(U_tar, U[, j]) > 0) U_tar else -U_tar
    }
    U = 0.7 * U + 0.3 * U_temp

    # Apply penalty if specified
    if (is.null(penalt)) {
      U_thresh = U
    } else if (penalt == "SCAD") {
      U_thresh = U
      if (gamma <= 2) {
        warning("Gamma needs to be > 2!")
        gamma = 2.01
      }
      for (j in 1:m) {
        U_thresh[, j] = U_thresh[, j] / sd(U_thresh[, j]) * sd1[j]
        U_thresh[, j] = SCAD_func(U_thresh[, j], lambda_ch = rho, gamma = gamma)
      }
    } else if (penalt == "Hardthreshold") {
      U_thresh = U * (abs(U) >= rho)
      for (j in 1:m) {
        if (sd(U_thresh[, j]) == 0) next
        U_thresh[, j] = U_thresh[, j] / sd(U_thresh[, j]) * sd1[j]
      }
    } else if (penalt == "L1") {
      U_thresh = sign(U) * (abs(U) - rho) * (abs(U) >= rho)
      for (j in 1:m) {
        U_thresh[, j] = U_thresh[, j] / sd(U_thresh[, j]) * sd1[j]
      }
    } else {
      print("No valid penalty method provided! No penalty was added")
      U_thresh = U
      for (j in 1:m) {
        U_thresh[, j] = U_thresh[, j] / sd(U_thresh[, j]) * sd1[j]
      }
    }

    # Update U using imputation method
    for (j in 1:m) {
      if (imput_method == "Previous") {
        remat = Ltrinv(U_thresh[, j], voxel, F) + diag(apply(abs(Ltrinv(U_thresh[, j], voxel, F)), 2, max))
      } else if (imput_method == "Average") {
        remat = Ltrinv(U_thresh[, j], voxel, F) + diag(apply(abs(Ltrinv(U_thresh[, j], voxel, F)), 2, mean))
      } else if (imput_method == "Max") {
        remat = Ltrinv(U_thresh[, j], voxel, F) + diag(apply(abs(Ltrinv(U_thresh[, j], voxel, F)), 2, max))
      } else if (imput_method == "No") {
        remat = Ltrinv(U_thresh[, j], voxel, F)
      } else {
        stop("No valid imputation method provided!")
      }

      eigenSl1 = eigen(remat)
      orderEigen = order(abs(eigenSl1$values), decreasing = TRUE)
      eigenset = orderEigen[1:R[j]]
      eigen_l[[j]] = eigenSl1$values[eigenset]
      J_l[[j]] = eigenSl1$vectors[, eigenset]
      U[, j] = Ltrans(eigenSl1$vectors[, eigenset] %*% diag(eigenSl1$values[eigenset]) %*% t(eigenSl1$vectors[, eigenset]))
    }

    if (!silent) {
      message(sprintf("Iter %d; Percentage change on U: %f; Percentage change on V: %f.",
                      Iter, sum((U - U_pre)^2) / sum(U^2), sum((V - V_pre)^2) / sum(V^2)))
    }
  }

  return(list(U = U, V = transT %*% V, R = R, CC = cor(Y %*% V, X %*% U)))
}
