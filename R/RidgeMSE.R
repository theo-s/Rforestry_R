#' @useDynLib forestryR

#' @title update_A_inv
#' @description This function is needed in the Ridge
#' @param A_inv_current current A matrix
#' @param new_x new feature to be added
#' @examples
#' set.seed(309814)
#' feat <- matrix(rnorm(80), ncol = 4)
#' lambda <- .2
#' feat_current <- feat[-20, ]
#' feat_current_c <- cbind(feat_current, 1)
#' A_inv_current <- solve(t(feat_current_c) %*% feat_current_c +
#'                          lambda * diag(rep(1, ncol(feat_current_c))))
#' new_x <- feat[20, ]
#' update_A_inv(A_inv_current = A_inv_current, new_x = new_x)
#' feat_c <- cbind(feat, 1)
#' solve(t(feat_c) %*% feat_c + lambda * diag(rep(1, ncol(feat_c)))) # truth
#' @return the updated A matrix
update_A_inv_R <- function(A_inv_current, new_x, leftnode) {
  if (!is.matrix(A_inv_current) ||
      nrow(A_inv_current) != ncol(A_inv_current)) {
    stop('A_inv_current should be a diagonal matrix')
  }
  if (nrow(A_inv_current) != length(new_x) + 1) {
    stop("new_x and A_inv_current should be dimensionwise compatible")
  }
  z <- A_inv_current %*% c(new_x, 1)
  g_L <- z %*% t(z) / as.numeric(1 + t(c(new_x, 1)) %*% z)
  g_R <- z %*% t(z) / as.numeric(1 - t(c(new_x, 1)) %*% z)
  if (leftnode)  return(A_inv_current - g_L)
  if (!leftnode) return(A_inv_current + g_R)
}



##
## Find_ridge_best_split_linear_MASS
##
#' @title Find_ridge_best_split_linear_MASS
#' @description Loops through all possible splits of feature
#' current.splitting.idx. This is just a test function which uses the ridge
#' regression canned ridge function of MASS
#' @param feat all features
#' @param y outcome vector
#' @param linear.idx feature indexes where a Ridge Regression should be fitted
#' @param current.splitting.idx the index of the feature we are currently
#' testing for splits.
#' @return It returns a vector of two doubles which determines what the best
#' splitting point is, and what it's MSE is.
#'   \item{\code{split_val_best}}{the best splitting value}
#'   \item{\code{MSE_best}}{the best MSE}
#' @examples
#' set.seed(309814)
#' feat <- matrix(rnorm(80), ncol = 4)
#' y <- rnorm(20)
#' linear.idx <- 1:4
#' current.splitting.idx <- 2
#' lambda <- .2
#' Find_ridge_best_split_linear_MASS(
#'   feat = feat,
#'   y = y,
#'   current.splitting.idx = current.splitting.idx,
#'   linear.idx = linear.idx,
#'   lambda = lambda
#' )
Find_ridge_best_split_linear_MASS <-
  function(feat,
           y,
           linear.idx,
           current.splitting.idx,
           lambda) {
    # sort by current.splitting.idx:
    split_feat_order <- order(feat[, current.splitting.idx])
    feat <- feat[split_feat_order, ]
    y <- y[split_feat_order]

    # sort by current.splitting.idx:
    split_val_best <- NA
    MSE_best <- Inf
    split_val_grid <- unique(feat[, current.splitting.idx])
    split_val_old <- split_val_grid[2]
    for (split_val_current in split_val_grid[-c(1:2, length(split_val_grid))]) {
      # split_val_current = split_val_grid[8]
      idx_L <-
        feat[, current.splitting.idx] < split_val_current
      idx_R <- !idx_L

      y_L <- y[idx_L]
      y_R <- y[idx_R]

      feat_L <- feat[idx_L, linear.idx]
      feat_R <- feat[idx_R, linear.idx]

      ridge_L <- MASS::lm.ridge(
        formula = y ~ .,
        data = data.frame(y = y_L, feat_L),
        lambda = lambda
      )
      ridge_R <- MASS::lm.ridge(
        formula = y ~ .,
        data = data.frame(y = y_R, feat_R),
        lambda = lambda
      )

      MSE_current <-
        mean(c((
          y_L - as.matrix(cbind(const = 1, feat_L)) %*% coef(ridge_L)
        ) ^ 2,
        (
          y_R - as.matrix(cbind(const = 1, feat_R)) %*% coef(ridge_R)
        ) ^ 2))

      if (MSE_best > MSE_current) {
        split_val_best <- stats::runif(1, split_val_old, split_val_current)
        MSE_best <- MSE_current
        bestSplitCount <- 1
      } else if (MSE_best == MSE_current) {
        bestSplitCount <- bestSplitCount + 1
        # Only update with probability 1/nseen
        if (stats::runif(1, 0, bestSplitCount) <= 1) {
          split_val_best <- stats::runif(1, split_val_old, split_val_current)
        }
      }
      split_val_old <- split_val_current
    }
    return(c('split_val_best' = split_val_best, 'MSE_best' = MSE_best))
  }

#' @title Find_ridge_best_split_slow
#' @description This is the fast impelementation of the best split function when
#' the leaf prediciton is coming from ridge regression.
#' @param feat all features
#' @param y outcome vector
#' @param linear.idx feature indexes where a Ridge Regression should be fitted
#' @param current.splitting.idx the index of the feature we are currently
#' testing for splits.
#' @return It returns a vector of two doubles which determines what the best
#' splitting point is, and what it's MSE is.
#'   \item{\code{split_val_best}}{the best splitting value}
#'   \item{\code{MSE_best}}{the best MSE}
#' @examples
#' set.seed(309814)
#' feat <- matrix(rnorm(80), ncol = 4)
#' y <- rnorm(20)
#' linear.idx <- 1:4
#' current.splitting.idx <- 2
#' lambda <- .2
#' Find_ridge_best_split_slow_linear_MASS(
#'   feat = feat,
#'   y = y,
#'   current.splitting.idx = current.splitting.idx,
#'   linear.idx = linear.idx,
#'   lambda = lambda
#' )
Find_ridge_best_split_slow  <-
  function(feat,
           y,
           linear.idx,
           current.splitting.idx,
           lambda) {
    # sort by current.splitting.idx:
    split_feat_order <- order(feat[, current.splitting.idx])
    feat <- as.matrix(feat[split_feat_order, ])
    y <- y[split_feat_order]
    # sort by current.splitting.idx:
    split_val_best <- NA
    MSE_best <- Inf
    split_val_grid <- unique(feat[, current.splitting.idx])
    split_val_old <- split_val_grid[2]
    for (split_val_current in split_val_grid[-c(1:2, length(split_val_grid))]) {
      # split_val_current = split_val_grid[17]
      idx_L <- feat[, current.splitting.idx] < split_val_current
      idx_R <- !idx_L

      y_L <- y[idx_L]
      y_R <- y[idx_R]

      feat_L <- feat[idx_L, linear.idx]
      feat_R <- feat[idx_R, linear.idx]
      #
      compute_coef <- function(feat, y) {
        A11 <- t(feat) %*% feat + lambda * diag(rep(1, ncol(feat)))
        A12 <- apply(feat, 2, sum)
        A22 <- nrow(feat)
        k <- ncol(feat)
        A <- matrix(NA, nrow = k + 1, ncol = k + 1)
        A[1:k, 1:k] <- A11
        A[k + 1, 1:k] <- A12
        A[1:k, k + 1] <- A12
        A[k + 1, k + 1] <- A22

        feat_c_L <- cbind(feat, 1)

        coef <- as.numeric(solve(A) %*% t(feat_c_L) %*% y)

        names(coef) <- c(paste0('X', 1:ncol(feat)), 'intercept')
        return(coef)
      }
      coef_L <- compute_coef(feat = feat_L, y = y_L)
      coef_R <- compute_coef(feat = feat_R, y = y_R)

      MSE_current <-
        mean(c((
          y_L - as.matrix(cbind(feat_L, const = 1)) %*% coef_L
        ) ^ 2,
        (
          y_R - as.matrix(cbind(feat_R, const = 1)) %*% coef_R
        ) ^ 2))

      if (MSE_best > MSE_current) {
        split_val_best <- stats::runif(1, split_val_old, split_val_current)
        MSE_best <- MSE_current
        bestSplitCount <- 1
      } else if (MSE_best == MSE_current) {
        bestSplitCount <- bestSplitCount + 1
        # Only update with probability 1/nseen
        if (stats::runif(1, 0, bestSplitCount) <= 1) {
          split_val_best <- stats::runif(1, split_val_old, split_val_current)
        }
      }
      split_val_old <- split_val_current
    }
    return(c('split_val_best' = split_val_best, 'MSE_best' = MSE_best))
  }


# ------------------------------------------------------------------------------
# helper function
#'@title Compute the A inverse matrix directly
compute_A_inv <- function(feat, lambda){
  A11 <- t(feat) %*% feat + lambda * diag(rep(1, ncol(feat)))
  A12 <- apply(feat, 2, sum)
  A22 <- nrow(feat)
  k <- ncol(feat)
  A <- matrix(NA, nrow = k + 1, ncol = k + 1)
  A[1:k, 1:k] <- A11
  A[k + 1, 1:k] <- A12
  A[1:k, k + 1] <- A12
  A[k + 1, k + 1] <- A22
  return(solve(A))
}


#'@title Compute the RSS based on the S, G and A terms
computeRSS <- function(S_L, S_R, A_inv_L, A_inv_R, G_L, G_R) {
  return(
    t(S_L) %*% A_inv_L %*% G_L %*% A_inv_L %*% S_L - 2 * t(S_L) %*% A_inv_L %*% S_L +
    t(S_R) %*% A_inv_R %*% G_R %*% A_inv_R %*% S_R - 2 * t(S_R) %*% A_inv_R %*% S_R
  )
}


#' @title Find_ridge_best_split
#' @description This is the fast impelementation of the best split function when
#' the leaf prediciton is coming from ridge regression.
#' @param feat all features
#' @param y outcome vector
#' @param linear.idx feature indexes where a Ridge Regression should be fitted
#' @param current.splitting.idx the index of the feature we are currently
#' testing for splits.
#' @return It returns a vector of two doubles which determines what the best
#' splitting point is, and what it's MSE is.
#'   \item{\code{split_val_best}}{the best splitting value}
#'   \item{\code{MSE_best}}{the best MSE}
#' @examples
#'set.seed(309814)
#'n <- 1000
#'feat <- matrix(rnorm(4 * n), ncol = 4)
#'feat[2,] <- feat[which.min(feat[,2]),]
#'feat[3,] <- feat[which.min(feat[,2]),]
#'feat[4,] <- feat[which.min(feat[,2]),]
#'feat[5,] <- feat[which.min(feat[,2]),] - 1
#'y <- rnorm(n)
#'linear.idx <- 1:4
#'current.splitting.idx <- 2
#'lambda <- .2
#'Find_ridge_best_split_linear_MASS(feat = feat, y = y, linear.idx = linear.idx,
#'                                  current.splitting.idx = current.splitting.idx,
#'                                  lambda = lambda)
#'Find_ridge_best_split_slow(feat = feat, y = y, linear.idx = linear.idx,
#'                           current.splitting.idx = current.splitting.idx,
#'                           lambda = lambda)
#'Find_ridge_best_split(feat = feat, y = y, linear.idx = linear.idx,
#'                      current.splitting.idx = current.splitting.idx,
#'                      lambda = lambda)
#' # Note that fast split is not estimating the RSS or the MSE, but the RSS -
#' # sum(y^2).
Find_ridge_best_split  <-
  function(feat,
           y,
           linear.idx,
           current.splitting.idx,
           lambda,
           update_A_inv = update_A_inv_R, 
           nodesize = list("splittingNodeSize" = 1,
                           "averagingNodeSize" = 1)) {
    # ---Sort by current.splitting.idx------------------------------------------
    split_feat_order <- order(feat[, current.splitting.idx])
    feat <- as.matrix(feat[split_feat_order, ])
    y <- y[split_feat_order]

    # ---Setup for the first possible split-------------------------------------
    # Get all units which are in the left node for the first possible split
    first_left_idx <-
      feat[, current.splitting.idx] <= feat[1, current.splitting.idx]

    A_inv_L <- compute_A_inv(feat = feat[first_left_idx, linear.idx,
                                         drop = FALSE],
                             lambda = lambda)
    A_inv_R <- compute_A_inv(feat = feat[!first_left_idx, linear.idx,
                                         drop = FALSE],
                             lambda = lambda)
    S_L <- t(cbind(feat[first_left_idx, linear.idx, drop = FALSE], 1)) %*%
      y[first_left_idx]
    S_R <- t(cbind(feat[!first_left_idx, linear.idx, drop = FALSE],1)) %*%
      y[!first_left_idx]
    G_L <- t(cbind(feat[first_left_idx, linear.idx, drop = FALSE],1)) %*%
             cbind(feat[first_left_idx, linear.idx, drop = FALSE],1)
    G_R <- t(cbind(feat[!first_left_idx, linear.idx, drop = FALSE],1)) %*%
             cbind(feat[!first_left_idx, linear.idx, drop = FALSE],1)

    RSS_best <- computeRSS(S_L, S_R, A_inv_L, A_inv_R, G_L, G_R)
    split_val_grid <- unique(feat[, current.splitting.idx])
    split_val_best <- stats::runif(1, split_val_grid[1], split_val_grid[2])
    bestSplitCount <- 1
    split_val_old <- split_val_grid[1]

    # ---Loop through all possible spllits and compute the RSS online-----------
    row_ptr <- 1 + sum(first_left_idx)
    for (split_val_current in split_val_grid[-1]) {
      # split_val_current = split_val_grid[2]
      while (feat[row_ptr, current.splitting.idx] < split_val_current) {
        x_toadd <- feat[row_ptr, linear.idx, drop = FALSE]
        y_toadd <- y[row_ptr]
        # update A
        A_inv_L <- update_A_inv(A_inv_current = A_inv_L, new_x = x_toadd,
                                leftnode = TRUE)
        A_inv_R <- update_A_inv(A_inv_current = A_inv_R, new_x = x_toadd,
                                leftnode = FALSE)
        S_L <- S_L + y_toadd * c(x_toadd, 1)
        S_R <- S_R - y_toadd * c(x_toadd, 1)
        G_L <- G_L + c(x_toadd, 1) %*% t(c(x_toadd, 1))
        G_R <- G_R - c(x_toadd, 1) %*% t(c(x_toadd, 1))
        row_ptr = row_ptr + 1
      }
      RSS_current <- computeRSS(S_L, S_R, A_inv_L, A_inv_R, G_L, G_R)
      if (RSS_best > RSS_current) {
        split_val_best <- stats::runif(1, split_val_old, split_val_current)
        RSS_best <- RSS_current
        bestSplitCount <- 1
      } else if (RSS_best == RSS_current) {
        bestSplitCount <- bestSplitCount + 1
        # Only update with probability 1/nseen
        if (stats::runif(1, 0, bestSplitCount) <= 1) {
          split_val_best <- stats::runif(1, split_val_old, split_val_current)
        }
      }
      # print(RSS_best)
      split_val_old <- split_val_current
    }
    return(c('split_val_best' = split_val_best, 'shifted_RSS_best' = RSS_best))
  }
