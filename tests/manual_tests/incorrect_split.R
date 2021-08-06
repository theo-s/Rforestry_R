f <- function(x) {
  if(x < .3) {
    return(x)
  } else {
    return(-2*x)
  }
}
set.seed(34545)
n <- 1000
feat <- matrix(rnorm(4 * n), ncol = 4)
y <- sapply(feat[,1], f)
linear.idx <- 1:4
current.splitting.idx <- 1
lambda <- .2

Find_ridge_best_split(feat = feat, y = y, linear.idx = linear.idx,
                     current.splitting.idx = current.splitting.idx,
                     lambda = lambda)
