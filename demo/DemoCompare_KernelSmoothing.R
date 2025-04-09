ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "Rcpp" = K2_Ep_rcpp(0.5),
    "Rpp_v1" = K2_Ep_rcpp_v1(0.5),
    "R" = K2_Ep(0.5)
  )
)

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- as.matrix(seq(-3, 3, 0.1))
w <- rexp(n = n, rate = 1)

test1 <- KDE_K2Ep_rcpp(X = X, x = x, h = 1.5)
test2 <- KDE_R(X = X, x = x, K = K2_Ep, h = 1.5)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Ep_rcpp(X = X, x = x, h = 1.5),
    R = KDE_R(X = X, x = x, K = K2_Ep, h = 1.5)
  )
)

test1 <- KDE_K2Ep_w_rcpp(X = X, x = x, h = 1.5, w = w)
test2 <- KDE_w_R(X = X, x = x, K = K2_Ep, h = 1.5, w = w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Ep_w_rcpp(X = X, x = x, h = 1.5, w = w),
    R = KDE_w_R(X = X, x = x, K = K2_Ep, h = 1.5, w = w)
  )
)

test1 <- KDE_K2Bw_rcpp(X = X, x = x, h = 1.5)
test2 <- KDE_R(X = X, x = x, K = K2_Bw, h = 1.5)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Bw_rcpp(X = X, x = x, h = 1.5),
    R = KDE_R(X = X, x = x, K = K2_Bw, h = 1.5)
  )
)

test1 <- KDE_K2Bw_w_rcpp(X = X, x = x, h = 1.5, w = w)
test2 <- KDE_w_R(X = X, x = x, K = K2_Bw, h = 1.5, w = w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Bw_w_rcpp(X = X, x = x, h = 1.5, w = w),
    R = KDE_w_R(X = X, x = x, K = K2_Bw, h = 1.5, w = w)
  )
)

test1 <- KDE_K4Bw_rcpp(X = X, x = x, h = 1.5)
test2 <- KDE_R(X = X, x = x, K = K4_Bw, h = 1.5)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K4Bw_rcpp(X = X, x = x, h = 1.5),
    R = KDE_R(X = X, x = x, K = K4_Bw, h = 1.5)
  )
)

test1 <- KDE_K4Bw_w_rcpp(X = X, x = x, h = 1.5, w = w)
test2 <- KDE_w_R(X = X, x = x, K = K4_Bw, h = 1.5, w = w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K4Bw_w_rcpp(X = X, x = x, h = 1.5, w = w),
    R = KDE_w_R(X = X, x = x, K = K4_Bw, h = 1.5, w = w)
  )
)
