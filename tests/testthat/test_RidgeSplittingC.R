library(testthat)
test_that("tests helper functions for cpp RidgeRSS calculator", {
  
  context("Ridge Splitting C")


  set.seed(3789)
  feat <- matrix(rnorm(80), ncol = 4)
  lambda <- .2
  feat_current <- feat[-20, ]
  feat_current_c <- cbind(feat_current, 1)
  A_inv_current <- solve(t(feat_current_c) %*% feat_current_c +
                           lambda * diag(rep(1, ncol(feat_current_c))))
  new_x <- feat[20, ]
  feat_c <- cbind(feat, 1)

  A_inv_truth <- solve(t(feat_c) %*% feat_c +
                         lambda * diag(rep(1, ncol(feat_c))))

  A_inv_C <- update_A_inv_C(
    a = A_inv_current, 
    new_x = matrix(new_x), 
    leftNode = TRUE)

  expect_equal(A_inv_C, A_inv_truth)




  set.seed(3789)
  prev <- matrix(rnorm(25), ncol = 5)
  nex <- matrix(rnorm(4), ncol = 1)
  nex <- rbind(nex,1)
  g_K_test <- update_G_k(prev, nex, TRUE)
  g_K_actual <- as.vector(nex) %*% t(as.vector(nex)) + prev
  expect_equal(g_K_actual, g_K_test, tolerance = 1e-6)

  g_K_test <- update_G_k(prev, nex, FALSE)
  g_K_actual <- prev - as.vector(nex) %*% t(as.vector(nex))
  expect_equal(g_K_actual, g_K_test, tolerance = 1e-6)



  set.seed(3789)
  prev <- matrix(rnorm(5), ncol = 1)
  nex <- matrix(rnorm(4), ncol = 1)
  s_K_test <- update_S_k(prev, nex, 59.54, TRUE)
  s_K_actual <- prev + (59.54 * rbind(nex,1))
  expect_equal(s_K_test, s_K_actual, tolerance = 1e-6)

  s_K_test <- update_S_k(prev, nex, 59.54, FALSE)
  s_K_actual <- prev - (59.54 * rbind(nex,1))
  expect_equal(s_K_test, s_K_actual, tolerance = 1e-6)

})

