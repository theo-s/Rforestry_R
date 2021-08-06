library(testthat)
test_that("tests functions needed for the Ridge splitting", {
  context("Ridge Splitting")
  
  set.seed(309814)
  feat <- matrix(rnorm(80), ncol = 4)
  lambda <- .2
  feat_current <- feat[-20, ]
  feat_current_c <- cbind(feat_current, 1)
  A_inv_current <- solve(t(feat_current_c) %*% feat_current_c +
                           lambda * diag(rep(1, ncol(feat_current_c))))
  new_x <- feat[20, ]
  A_inv_algo <- update_A_inv_R(A_inv_current = A_inv_current, new_x = new_x,
                             leftnode = TRUE)
  feat_c <- cbind(feat, 1)
  A_inv_truth <- solve(t(feat_c) %*% feat_c +
                         lambda * diag(rep(1, ncol(feat_c)))) # truth

  expect_equal(A_inv_algo, A_inv_truth)


  set.seed(309814)
  feat <- matrix(rnorm(80), ncol = 4)
  y <- rnorm(20)
  linear.idx <- 1:4
  current.splitting.idx <- 2
  lambda <- .2
  best_split <-
    Find_ridge_best_split_linear_MASS(
      feat = feat,
      y = y,
      current.splitting.idx = current.splitting.idx,
      linear.idx = linear.idx,
      lambda = lambda
    )
  expect_equal(best_split,
               c(
                 'split_val_best' = 1.1028234,
                 'MSE_best' = 0.1785655
               ), tolerance = 1e-7)



  set.seed(309814)
  n <- 100
  feat <- matrix(rnorm(4 * n), ncol = 4)
  y <- rnorm(n)
  linear.idx <- 1:4
  current.splitting.idx <- 2
  lambda <- .2

  split_pts <-
    Find_ridge_best_split_slow(
      feat = feat,
      y = y,
      linear.idx = linear.idx,
      current.splitting.idx = current.splitting.idx,
      lambda = lambda
    )
  expect_equal(split_pts,  c(
    'split_val_best' = -0.6049398,
    'MSE_best' = 0.9093351
  ), tolerance = 1e-7)


  set.seed(309814)
  n <- 1000
  feat <- matrix(rnorm(4 * n), ncol = 4)
  feat[2,] <- feat[which.min(feat[,2]),]
  feat[3,] <- feat[which.min(feat[,2]),]
  feat[4,] <- feat[which.min(feat[,2]),]
  feat[5,] <- feat[which.min(feat[,2]),] - 1
  y <- rnorm(n)
  linear.idx <- 1:4
  current.splitting.idx <- 2
  lambda <- .2
  fast_solver_sol <- Find_ridge_best_split(
    feat = feat,
    y = y,
    linear.idx = linear.idx,
    current.splitting.idx = current.splitting.idx,
    lambda = lambda
  )
  expect_equal(fast_solver_sol,  c(
    'split_val_best' = 1.32333,
    'shifted_RSS_best' = -17.33089
  ), tolerance = 1e-6)

  
  # ----------------------------------------------------------------------------
  context('SelectBestFeatureLinear')
  
  set.seed(309814)
  feat <- matrix(rnorm(80), ncol = 4)
  lambda <- .2
  
  x = feat
  y = x[, 1] + ifelse(x[, 2] > .1,-10, 0)
  featureList_spit = c(1, 3, 4)
  featureList_lin = c(1, 3)
  lambda = 1
  sampleIndex = list(
    "averagingSampleIndex" = 1:length(y),
    "splittingSampleIndex" = 1:length(y)
  )
  nodesize = list("splittingNodeSize" = 5,
                  "averagingNodeSize" = 5)
  categoricalFeatureCols = list()
  
  selectbf <- selectBestFeatureLinear(
    x,
    y,
    featureList_spit,
    featureList_lin,
    lambda,
    sampleIndex = list(
      "averagingSampleIndex" = 1:length(y),
      "splittingSampleIndex" = 1:length(y)
    ),
    nodesize = list(
      "splittingNodeSize" = 5,
      "averagingNodeSize" = 5
    ),
    categoricalFeatureCols = list()
  )
  
  expect_equal(selectbf, list("bestSplitFeature" = 3,
                              "bestSplitValue" = 0.1691263), tolerance = 1e-6)
    
})

