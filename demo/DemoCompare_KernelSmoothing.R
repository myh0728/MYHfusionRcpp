ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "Rcpp" = K2_Ep_rcpp(0.5),
    "Rpp_v1" = K2_Ep_rcpp_v1(0.5),
    "R" = K2_Ep(0.5)
  )
)

n <- 200
p <- 2

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- matrix(rnorm(10 * p), nrow = n, ncol = p)
w <- rexp(n = n, rate = 1)
Y <- as.vector(sin(X[, 1] * 2)) + rnorm(n, mean = 0, sd = 0.2)
y <- seq(min(Y), max(Y), length.out = 50)

test1 <- KDE_K2Ep_rcpp(X = X, x = x, h = rep(1.5, p))
test2 <- KDE_R(X = X, x = x, K = K2_Ep, h = rep(1.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Ep_rcpp(X = X, x = x, h = rep(1.5, p)),
    R = KDE_R(X = X, x = x, K = K2_Ep, h = rep(1.5, p))
  )
)

test1 <- KDE_K2Ep_w_rcpp(X = X, x = x, h = rep(1.5, p), w = w)
test2 <- KDE_w_R(X = X, x = x, K = K2_Ep, h = rep(1.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Ep_w_rcpp(X = X, x = x, h = rep(1.5, p), w = w),
    R = KDE_w_R(X = X, x = x, K = K2_Ep, h = rep(1.5, p), w = w)
  )
)

test1 <- KDE_K2Bw_rcpp(X = X, x = x, h = rep(1.5, p))
test2 <- KDE_R(X = X, x = x, K = K2_Bw, h = rep(1.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Bw_rcpp(X = X, x = x, h = rep(1.5, p)),
    R = KDE_R(X = X, x = x, K = K2_Bw, h = rep(1.5, p))
  )
)

test1 <- KDE_K2Bw_w_rcpp(X = X, x = x, h = rep(1.5, p), w = w)
test2 <- KDE_w_R(X = X, x = x, K = K2_Bw, h = rep(1.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2Bw_w_rcpp(X = X, x = x, h = rep(1.5, p), w = w),
    R = KDE_w_R(X = X, x = x, K = K2_Bw, h = rep(1.5, p), w = w)
  )
)

test1 <- KDE_K4Bw_rcpp(X = X, x = x, h = rep(1.5, p))
test2 <- KDE_R(X = X, x = x, K = K4_Bw, h = rep(1.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K4Bw_rcpp(X = X, x = x, h = rep(1.5, p)),
    R = KDE_R(X = X, x = x, K = K4_Bw, h = rep(1.5, p))
  )
)

test1 <- KDE_K4Bw_w_rcpp(X = X, x = x, h = rep(1.5, p), w = w)
test2 <- KDE_w_R(X = X, x = x, K = K4_Bw, h = rep(1.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K4Bw_w_rcpp(X = X, x = x, h = rep(1.5, p), w = w),
    R = KDE_w_R(X = X, x = x, K = K4_Bw, h = rep(1.5, p), w = w)
  )
)

test1 <- KNW_K2Ep_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p))
test2 <- KNW_R(Y = Y, X = X, x = x, K = K2_Ep, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Ep_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p)),
    R = KNW_R(Y = Y, X = X, x = x, K = K2_Ep, h = rep(0.5, p))
  )
)

test1 <- KNW_K2Ep_w_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p), w = w)
test2 <- KNW_w_R(Y = Y, X = X, x = x, K = K2_Ep, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Ep_w_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p), w = w),
    R = KNW_w_R(Y = Y, X = X, x = x, K = K2_Ep, h = rep(0.5, p), w = w)
  )
)

test1 <- KNW_K2Bw_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p))
test2 <- KNW_R(Y = Y, X = X, x = x, K = K2_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Bw_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p)),
    R = KNW_R(Y = Y, X = X, x = x, K = K2_Bw, h = rep(0.5, p))
  )
)

test1 <- KNW_K2Bw_w_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p), w = w)
test2 <- KNW_w_R(Y = Y, X = X, x = x, K = K2_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Bw_w_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p), w = w),
    R = KNW_w_R(Y = Y, X = X, x = x, K = K2_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- KNW_K4Bw_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p))
test2 <- KNW_R(Y = Y, X = X, x = x, K = K4_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K4Bw_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p)),
    R = KNW_R(Y = Y, X = X, x = x, K = K4_Bw, h = rep(0.5, p))
  )
)

test1 <- KNW_K4Bw_w_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p), w = w)
test2 <- KNW_w_R(Y = Y, X = X, x = x, K = K4_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K4Bw_w_rcpp(Y = Y, X = X, x = x, h = rep(0.5, p), w = w),
    R = KNW_w_R(Y = Y, X = X, x = x, K = K4_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- KNWcdf_K2Ep_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p))
test2 <- KNWcdf_R(Y = Y, y = y , X = X, x = x, K = K2_Ep, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNWcdf_K2Ep_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p)),
    R = KNWcdf_R(Y = Y, y = y , X = X, x = x, K = K2_Ep, h = rep(0.5, p))
  )
)

test1 <- KNWcdf_K2Bw_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p))
test2 <- KNWcdf_R(Y = Y, y = y , X = X, x = x, K = K2_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNWcdf_K2Bw_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p)),
    R = KNWcdf_R(Y = Y, y = y , X = X, x = x, K = K2_Bw, h = rep(0.5, p))
  )
)

test1 <- KNWcdf_K4Bw_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p))
test2 <- KNWcdf_R(Y = Y, y = y , X = X, x = x, K = K4_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNWcdf_K4Bw_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p)),
    R = KNWcdf_R(Y = Y, y = y , X = X, x = x, K = K4_Bw, h = rep(0.5, p))
  )
)

test1 <- KNWcdf_K2Ep_w_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p), w = w)
test2 <- KNWcdf_w_R(Y = Y, y = y , X = X, x = x, K = K2_Ep, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNWcdf_K2Ep_w_rcpp(Y = Y, y = y , X = X, x = x, h = rep(0.5, p), w = w),
    R = KNWcdf_w_R(Y = Y, y = y , X = X, x = x, K = K2_Ep, h = rep(0.5, p), w = w)
  )
)

test1 <- CVKNW_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- CVKNW_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = CVKNW_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
  )
)

test1 <- CVKNW_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- CVKNW_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = CVKNW_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
  )
)

test1 <- CVKNW_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- CVKNW_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = CVKNW_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
  )
)

test1 <- CVKNW_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- CVKNW_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = CVKNW_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- CVKNW_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- CVKNW_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = CVKNW_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
  )
)

test1 <- CVKNW_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- CVKNW_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = CVKNW_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- CVKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- CVKNWcdf_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = CVKNWcdf_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
  )
)

test1 <- CVKNWcdf_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- CVKNWcdf_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNWcdf_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = CVKNWcdf_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
  )
)

test1 <- CVKNWcdf_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- CVKNWcdf_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNWcdf_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = CVKNWcdf_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
  )
)

test1 <- CVKNWcdf_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- CVKNWcdf_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNWcdf_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = CVKNWcdf_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- CVKNWcdf_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- CVKNWcdf_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNWcdf_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = CVKNWcdf_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
  )
)

test1 <- CVKNWcdf_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- CVKNWcdf_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNWcdf_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = CVKNWcdf_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- LSKNW_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- LSKNW_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNW_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = LSKNW_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
  )
)

test1 <- LSKNW_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- LSKNW_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNW_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = LSKNW_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
  )
)

test1 <- LSKNW_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- LSKNW_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNW_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = LSKNW_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
  )
)

test1 <- LSKNW_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- LSKNW_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNW_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = LSKNW_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- LSKNW_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- LSKNW_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNW_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = LSKNW_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
  )
)

test1 <- LSKNW_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- LSKNW_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNW_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = LSKNW_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- LSKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- LSKNWcdf_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = LSKNWcdf_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p))
  )
)

test1 <- LSKNWcdf_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- LSKNWcdf_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNWcdf_K2Ep_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = LSKNWcdf_w_R(Y = Y, X = X, K = K2_Ep, h = rep(0.5, p), w = w)
  )
)

test1 <- LSKNWcdf_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- LSKNWcdf_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNWcdf_K2Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = LSKNWcdf_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p))
  )
)

test1 <- LSKNWcdf_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- LSKNWcdf_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNWcdf_K2Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = LSKNWcdf_w_R(Y = Y, X = X, K = K2_Bw, h = rep(0.5, p), w = w)
  )
)

test1 <- LSKNWcdf_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p))
test2 <- LSKNWcdf_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNWcdf_K4Bw_rcpp(Y = Y, X = X, h = rep(0.5, p)),
    R = LSKNWcdf_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p))
  )
)

test1 <- LSKNWcdf_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w)
test2 <- LSKNWcdf_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = LSKNWcdf_K4Bw_w_rcpp(Y = Y, X = X, h = rep(0.5, p), w = w),
    R = LSKNWcdf_w_R(Y = Y, X = X, K = K4_Bw, h = rep(0.5, p), w = w)
  )
)

