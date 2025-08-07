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
Y <- as.vector(sin(X[, 1] * 2)) + rnorm(n, mean = 0, sd = 0.2)

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

test1 <- KNW_K2Ep_rcpp(Y = Y, X = X, x = x, h = 0.5)
test2 <- KNW_R(Y = Y, X = X, x = x, K = K2_Ep, h = 0.5)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Ep_rcpp(Y = Y, X = X, x = x, h = 0.5),
    R = KNW_R(Y = Y, X = X, x = x, K = K2_Ep, h = 0.5)
  )
)

test1 <- KNW_K2Ep_w_rcpp(Y = Y, X = X, x = x, h = 0.5, w = w)
test2 <- KNW_w_R(Y = Y, X = X, x = x, K = K2_Ep, h = 0.5, w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Ep_w_rcpp(Y = Y, X = X, x = x, h = 0.5, w = w),
    R = KNW_w_R(Y = Y, X = X, x = x, K = K2_Ep, h = 0.5, w = w)
  )
)

test1 <- KNW_K2Bw_rcpp(Y = Y, X = X, x = x, h = 0.5)
test2 <- KNW_R(Y = Y, X = X, x = x, K = K2_Bw, h = 0.5)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Bw_rcpp(Y = Y, X = X, x = x, h = 0.5),
    R = KNW_R(Y = Y, X = X, x = x, K = K2_Bw, h = 0.5)
  )
)

test1 <- KNW_K2Bw_w_rcpp(Y = Y, X = X, x = x, h = 0.5, w = w)
test2 <- KNW_w_R(Y = Y, X = X, x = x, K = K2_Bw, h = 0.5, w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K2Bw_w_rcpp(Y = Y, X = X, x = x, h = 0.5, w = w),
    R = KNW_w_R(Y = Y, X = X, x = x, K = K2_Bw, h = 0.5, w = w)
  )
)

test1 <- KNW_K4Bw_rcpp(Y = Y, X = X, x = x, h = 0.5)
test2 <- KNW_R(Y = Y, X = X, x = x, K = K4_Bw, h = 0.5)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K4Bw_rcpp(Y = Y, X = X, x = x, h = 0.5),
    R = KNW_R(Y = Y, X = X, x = x, K = K4_Bw, h = 0.5)
  )
)

test1 <- KNW_K4Bw_w_rcpp(Y = Y, X = X, x = x, h = 0.5, w = w)
test2 <- KNW_w_R(Y = Y, X = X, x = x, K = K4_Bw, h = 0.5, w = w)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KNW_K4Bw_w_rcpp(Y = Y, X = X, x = x, h = 0.5, w = w),
    R = KNW_w_R(Y = Y, X = X, x = x, K = K4_Bw, h = 0.5, w = w)
  )
)

test1 <- CVKNW_K2Ep_rcpp(Y = Y, X = X, h = 0.5)
test2 <- CVKNW_R(Y = Y, X = X, K = K2_Ep, h = 0.5)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K2Ep_rcpp(Y = Y, X = X, h = 0.5),
    R = CVKNW_R(Y = Y, X = X, K = K2_Ep, h = 0.5)
  )
)

test1 <- CVKNW_K2Bw_rcpp(Y = Y, X = X, h = 0.5)
test2 <- CVKNW_R(Y = Y, X = X, K = K2_Bw, h = 0.5)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K2Bw_rcpp(Y = Y, X = X, h = 0.5),
    R = CVKNW_R(Y = Y, X = X, K = K2_Bw, h = 0.5)
  )
)

test1 <- CVKNW_K4Bw_rcpp(Y = Y, X = X, h = 0.5)
test2 <- CVKNW_R(Y = Y, X = X, K = K4_Bw, h = 0.5)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNW_K4Bw_rcpp(Y = Y, X = X, h = 0.5),
    R = CVKNW_R(Y = Y, X = X, K = K4_Bw, h = 0.5)
  )
)

test1 <- CVKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = 0.5)
test2 <- CVKNWcdf_R(Y = Y, X = X, K = K2_Ep, h = 0.5)
sum(abs(test1 - test2))
ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = CVKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = 0.5),
    R = CVKNWcdf_R(Y = Y, X = X, K = K2_Ep, h = 0.5)
  )
)



#

cv_fun_K2Ep <- function(h.log) {

  cv_value <- CVKNW_K2Ep_rcpp(Y = Y, X = X, h = exp(h.log))

  return(cv_value)
}

nlminb(start = 0,
       objective = cv_fun_K2Ep)

h.grid <- seq(0.01, 1, 0.01)
cv.grid <- rep(0, length(h.grid))
for (l in 1:length(h.grid))
{
  cv.grid[l] <- CVKNW_K2Ep_rcpp(Y = Y, X = X, h = h.grid[l])
}
plot(h.grid, cv.grid, type = 'l')

#

cv_fun_K2Bw <- function(h.log) {

  cv_value <- CVKNW_K2Bw_rcpp(Y = Y, X = X, h = exp(h.log))

  return(cv_value)
}

nlminb(start = 0,
       objective = cv_fun_K2Bw)

h.grid <- seq(0.01, 1, 0.01)
cv.grid <- rep(0, length(h.grid))
for (l in 1:length(h.grid))
{
  cv.grid[l] <- CVKNW_K2Bw_rcpp(Y = Y, X = X, h = h.grid[l])
}
plot(h.grid, cv.grid, type = 'l')

#

cv_fun_K4Bw <- function(h.log) {

  cv_value <- CVKNW_K4Bw_rcpp(Y = Y, X = X, h = exp(h.log))

  return(cv_value)
}

nlminb(start = 0,
       objective = cv_fun_K4Bw)

h.grid <- seq(0.01, 2, 0.01)
cv.grid <- rep(0, length(h.grid))
for (l in 1:length(h.grid))
{
  cv.grid[l] <- CVKNW_K4Bw_rcpp(Y = Y, X = X, h = h.grid[l])
}
plot(h.grid, cv.grid, type = 'l', ylim = c(0, 0.2))

#






