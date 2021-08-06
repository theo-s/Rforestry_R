set.seed(2343221)
n <- 100
feat <- data.frame(x1 = rnorm(n),
                   x2 = rnorm(n),
                   x3 = rnorm(n),
                   x4 = rnorm(n),
                   x5 = rnorm(n))
f <- function(x){
  x[x > .3] <- -2 * x[x > .3] 
  return(x)
}
y <- f(feat$x1)
lambda = .3
current.splitting.idx = 1
linear.idx = 1:4

Find_ridge_best_split(feat, y, linear.idx = linear.idx, 
                      current.splitting.idx = current.splitting.idx, 
                      lambda = lambda)
Find_ridge_best_split_linear_MASS(feat, y, linear.idx = linear.idx, 
                                  current.splitting.idx = current.splitting.idx, 
                                  lambda = lambda)
Find_ridge_best_split_slow(feat, y, linear.idx = linear.idx, 
                                  current.splitting.idx = current.splitting.idx, 
                                  lambda = lambda)



