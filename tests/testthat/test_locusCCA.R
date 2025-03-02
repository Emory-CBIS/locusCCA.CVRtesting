library(testthat)
library(locusCCA.CVRtesting)
library(MASS)  # For ginv() function


test_that("Low_rank_CCA runs correctly with real dataset", {

  # Load the dataset
  S_real_agg <- readRDS(system.file("data", "S_real_agg.rds", package = "locusCCA.CVRtesting"))
  original_Y <- readRDS(system.file("data", "original_Y.rds", package = "locusCCA.CVRtesting"))

  # Define parameters
  n <- 300  # Sample size
  m <- 6    # Number of canonical components
  p <- ncol(S_real_agg)  # Number of predictors in X
  q <- 10   # Response variables in Y

  # Generate synthetic signals based on the real dataset
  U <- t(S_real_agg[ 1:m,]) / 1000
  sample1 <- sample(2:13, q)
  eigen_Y <- eigen(t(original_Y[, sample1]) %*% original_Y[, sample1])
  V <- eigen_Y$vectors[ sample(1:10, q),sample(1:10, m)]

  # Simulate X and Y using known structures
  set.seed(111)
  fx <- matrix(rnorm(n * m, 0, 1), nrow = n)
  fy <- fx + matrix(rnorm(n * m, 0, 0.6), nrow = n)
  X <- 500 * fx[, 1:m] %*% (ginv(U)) + matrix(rnorm(n * p, 0, 0.01), nrow = n)
  Y <- fy[, 1:m] %*% ginv(V) + matrix(rnorm(n * q, 0, 0.01), nrow = n)

  # Run Low_rank_CCA on real data
  result <- Locus_CCA(X, Y, voxel = 264, m = m, lambda = 0.008,
                         penalt = "Hardthreshold", proportion = 0.95,
                          silent = FALSE, tol = 1e-3)

  # Assertions for correctness
  expect_true(is.list(result))
  expect_true("U" %in% names(result))
  expect_true("V" %in% names(result))
  expect_true("CC" %in% names(result))
  expect_equal(dim(result$U), c(p, m))
  expect_equal(dim(result$V), c(q, m))
})

test_that("score_testing runs correctly with real dataset", {

  # Load the dataset
  S_real_agg <- readRDS(system.file("data", "S_real_agg.rds", package = "locusCCA.CVRtesting"))
  original_Y <- readRDS(system.file("data", "original_Y.rds", package = "locusCCA.CVRtesting"))

  # Define parameters
  n <- 300
  q <- 10
  p <- ncol(S_real_agg)
  m <- 6
  # Simulate X and Y using known structures
  # Generate synthetic signals based on the real dataset
  U <- t(S_real_agg[ 1:m,]) / 1000
  sample1 <- sample(2:13, q)
  eigen_Y <- eigen(t(original_Y[, sample1]) %*% original_Y[, sample1])
  V <- eigen_Y$vectors[ sample(1:10, q),sample(1:10, m)]

  # Simulate X and Y using known structures
  set.seed(111)
  fx <- matrix(rnorm(n * m, 0, 1), nrow = n)
  fy <- fx + matrix(rnorm(n * m, 0, 0.6), nrow = n)
  X <- 500 * fx[, 1:m] %*% (ginv(U)) + matrix(rnorm(n * p, 0, 0.01), nrow = n)
  Y <- fy[, 1:m] %*% ginv(V) + matrix(rnorm(n * q, 0, 0.01), nrow = n)
  result <- Locus_CCA(X, Y, voxel = 264, m = m, lambda = 0.008,
                      penalt = "Hardthreshold", proportion = 0.95,
                      silent = FALSE, tol = 1e-3)

  component = sample(1:6,2)
  weights = rnorm(2,1,0.1)
  beta =2000*apply((U[,component] %*% diag(weights)), 1, sum)
  z = X %*% beta + rnorm(n,sd = 0.1)

  scores <- CVR_testing(result$U, X, z)

  # Assertions for correctness
  expect_length(scores, m)
  expect_true(all(is.finite(scores)))  # Ensure no NaN/Inf values
})
