n <- 1000
p <- 3

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta0 <- c(1, rep(0.5, p - 1))
Y <- as.vector(sin(X %*% beta0)) + rnorm(n, mean = 0, sd = 0.2)

test1 <- cumuSIR_rcpp(X = X, Y = Y, eps_inv = 1e-7)
test2 <- cumuSIR_R(X = X, Y = Y, eps.inv = 1e-7)
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "Rcpp" = cumuSIR_rcpp(X = X, Y = Y, eps_inv = 1e-7),
    "R" = cumuSIR_R(X = X, Y = Y, eps.inv = 1e-7)
  )
)


