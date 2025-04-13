##### Normal regression model #####

p <- 2

theta0 <- 1
theta1 <- 1
theta2 <- 1
sigma0 <- 0.5

n <- 100
N <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rnorm(n = n, mean = theta0 + X1 * theta1 + X2 * theta2, sd = sigma0)
test.data <- data.frame(response = Y,
                        covariate = X)

X1_shift <- rnorm(n = N, mean = 0, sd = 1)
X2_shift <- rbinom(n = N, size = 1, prob = 0.6)
X_shift <- cbind(X1_shift, X2_shift)
Y_shift <- rnorm(n = N, mean = theta0 + X1_shift * theta1 + X2_shift * theta2, sd = sigma0)

### auxiliary information: average of X

X <- as.matrix(X)
phi <- colMeans(X_shift)

test1 <- ADPPS_EX_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta = rep(0.1, length(phi)))

test2 <- ADPPS_EX_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EX_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta = rep(0.1, length(phi))),
    "Rcpp" = ADPPS_EX_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADPPS_EX_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EX_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, eta_initial = rep(0.1, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- ADPPS_EX_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EX_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EX_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, eta_initial = rep(0.1, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = ADPPS_EX_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1)

test2 <- ADPPS_EX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1),
    "Rcpp" = ADPPS_EX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1)
  )
)

test1 <- ADPPSvar_EX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta = rep(0.1, p))

test2 <- ADPPSvar_EX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, eta = rep(0.1, p))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPSvar_EX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta = rep(0.1, p)),
    "Rcpp" = ADPPSvar_EX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, eta = rep(0.1, p))
  )
)

### auxiliary information: average of Y

X <- as.matrix(X)
phi <- mean(Y_shift)

test1 <- AD_EY_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, eta = 0.1)

test2 <- AD_EY_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EY_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, eta = 0.1),
    "Rcpp" = AD_EY_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, eta = 0.1)
  )
)

test1 <- AD_EY_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EY_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, eta_initial = 0,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- AD_EY_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EY_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EY_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, eta_initial = 0,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = AD_EY_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- AD_EY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi)

test2 <- AD_EY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi),
    "Rcpp" = AD_EY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi)
  )
)

test1 <- ADvar_EY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, eta = 0.1)

test2 <- ADvar_EY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADvar_EY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, eta = 0.1),
    "Rcpp" = ADvar_EY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, eta = 0.1)
  )
)

test1 <- ADCS_EY_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), eta = 0.1)

test2 <- ADCS_EY_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EY_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), eta = 0.1),
    "Rcpp" = ADCS_EY_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), eta = 0.1)
  )
)

test1 <- ADCS_EY_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), eta.initial = 0.1,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EY_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), eta_initial = 0.1,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- ADCS_EY_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), eta.initial = 0.1,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EY_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), eta.initial = 0.1,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EY_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), eta_initial = 0.1,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = ADCS_EY_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), eta.initial = 0.1,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p))

test2 <- ADCS_EY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p)),
    "Rcpp" = ADCS_EY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p))
  )
)

test1 <- ADCSvar_EY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), eta = 0.1)

test2 <- ADCSvar_EY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCSvar_EY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), eta = 0.1),
    "Rcpp" = ADCSvar_EY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), eta = 0.1)
  )
)

test1 <- ADPPS_EY_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta = 0.1)

test2 <- ADPPS_EY_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EY_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta = 0.1),
    "Rcpp" = ADPPS_EY_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, eta = 0.1)
  )
)

test1 <- ADPPS_EY_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta.initial = 0.1,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EY_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, eta_initial = 0.1,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- ADPPS_EY_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta.initial = 0.1,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EY_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta.initial = 0.1,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EY_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, eta_initial = 0.1,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = ADPPS_EY_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta.initial = 0.1,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1)

test2 <- ADPPS_EY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1),
    "Rcpp" = ADPPS_EY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1)
  )
)

test1 <- ADPPSvar_EY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, eta = 0.1)

test2 <- ADPPSvar_EY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPSvar_EY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, eta = 0.1),
    "Rcpp" = ADPPSvar_EY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, eta = 0.1)
  )
)

### auxiliary information: average of X given Y

X <- as.matrix(X)

y.pts <- matrix(0, 2, 2)
y.pts[1, ] <- quantile(Y_shift, c(0.0, 0.5))
y.pts[2, ] <- quantile(Y_shift, c(0.5, 1.0))

phi <- matrix(0, 2, p)
phi[1, ] <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) & (Y_shift <= y.pts[1, 2]), ])
phi[2, ] <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) & (Y_shift <= y.pts[2, 2]), ])

test1 <- AD_EXsubY_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y.pts = y.pts, eta = rep(0.1, length(phi)))

test2 <- AD_EXsubY_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y_pts = y.pts, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EXsubY_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y.pts = y.pts, eta = rep(0.1, length(phi))),
    "Rcpp" = AD_EXsubY_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y_pts = y.pts, eta = rep(0.1, length(phi)))
  )
)

test1 <- AD_EXsubY_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y.pts = y.pts,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EXsubY_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y_pts = y.pts,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- AD_EXsubY_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y.pts = y.pts,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EXsubY_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y.pts = y.pts,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EXsubY_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y_pts = y.pts,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = AD_EXsubY_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y.pts = y.pts,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- AD_EXsubY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y.pts = y.pts)

test2 <- AD_EXsubY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y_pts = y.pts)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EXsubY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y.pts = y.pts),
    "Rcpp" = AD_EXsubY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y_pts = y.pts)
  )
)

test1 <- ADvar_EXsubY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y.pts = y.pts, eta = rep(0.1, length(phi)))

test2 <- ADvar_EXsubY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, y_pts = y.pts, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADvar_EXsubY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y.pts = y.pts, eta = rep(0.1, length(phi))),
    "Rcpp" = ADvar_EXsubY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, y_pts = y.pts, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADCS_EXsubY_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts, eta = rep(0.1, length(phi)))

test2 <- ADCS_EXsubY_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EXsubY_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts, eta = rep(0.1, length(phi))),
    "Rcpp" = ADCS_EXsubY_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADCS_EXsubY_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts,
  eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EXsubY_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts,
  eta_initial = rep(0.1, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- ADCS_EXsubY_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts,
  eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EXsubY_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts,
      eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EXsubY_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts,
      eta_initial = rep(0.1, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = ADCS_EXsubY_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts,
      eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EXsubY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts)

test2 <- ADCS_EXsubY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EXsubY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts),
    "Rcpp" = ADCS_EXsubY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts)
  )
)

test1 <- ADCSvar_EXsubY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts, eta = rep(0.1, length(phi)))

test2 <- ADCSvar_EXsubY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCSvar_EXsubY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts, eta = rep(0.1, length(phi))),
    "Rcpp" = ADCSvar_EXsubY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADPPS_EXsubY_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta = rep(0.1, length(phi)))

test2 <- ADPPS_EXsubY_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, y_pts = y.pts, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EXsubY_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta = rep(0.1, length(phi))),
    "Rcpp" = ADPPS_EXsubY_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, y_pts = y.pts, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADPPS_EXsubY_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EXsubY_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, y_pts = y.pts, eta_initial = rep(0.1, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- ADPPS_EXsubY_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EXsubY_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EXsubY_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, y_pts = y.pts, eta_initial = rep(0.1, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = ADPPS_EXsubY_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EXsubY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, y.pts = y.pts)

test2 <- ADPPS_EXsubY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, y_pts = y.pts)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EXsubY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, y.pts = y.pts),
    "Rcpp" = ADPPS_EXsubY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, y_pts = y.pts)
  )
)

test1 <- ADPPSvar_EXsubY_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta = rep(0.1, length(phi)))

test2 <- ADPPSvar_EXsubY_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, y_pts = y.pts, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPSvar_EXsubY_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, y.pts = y.pts, eta = rep(0.1, length(phi))),
    "Rcpp" = ADPPSvar_EXsubY_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, y_pts = y.pts, eta = rep(0.1, length(phi)))
  )
)

### auxiliary information: average of Y given X

X <- as.matrix(X)

phi <- c(mean(Y_shift[X_shift[, 1] > 0]),
         mean(Y_shift[X_shift[, 2] > 0]))
inclusion <- cbind(X[, 1] > 0,
                   X[, 2] <= 0)

test1 <- AD_EYsubX_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi)))

test2 <- AD_EYsubX_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EYsubX_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi))),
    "Rcpp" = AD_EYsubX_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi)))
  )
)

test1 <- AD_EYsubX_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EYsubX_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- AD_EYsubX_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EYsubX_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EYsubX_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = AD_EYsubX_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- AD_EYsubX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion)

test2 <- AD_EYsubX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EYsubX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion),
    "Rcpp" = AD_EYsubX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion)
  )
)

test1 <- ADvar_EYsubX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi)))

test2 <- ADvar_EYsubX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADvar_EYsubX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi))),
    "Rcpp" = ADvar_EYsubX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, inclusion = inclusion, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADCS_EYsubX_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion, eta = rep(0.1, length(phi)))

test2 <- ADCS_EYsubX_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EYsubX_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion, eta = rep(0.1, length(phi))),
    "Rcpp" = ADCS_EYsubX_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADCS_EYsubX_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
  eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EYsubX_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion,
  eta_initial = rep(0.1, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- ADCS_EYsubX_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
  eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EYsubX_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
      eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EYsubX_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion,
      eta_initial = rep(0.1, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = ADCS_EYsubX_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
      eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EYsubX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion)

test2 <- ADCS_EYsubX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EYsubX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion),
    "Rcpp" = ADCS_EYsubX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion)
  )
)

test1 <- ADCSvar_EYsubX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion, eta = rep(0, length(phi)))

test2 <- ADCSvar_EYsubX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion, eta = rep(0, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCSvar_EYsubX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion, eta = rep(0, length(phi))),
    "Rcpp" = ADCSvar_EYsubX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion, eta = rep(0, length(phi)))
  )
)

test1 <- ADPPS_EYsubX_normal_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi)))

test2 <- ADPPS_EYsubX_normal_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EYsubX_normal_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi))),
    "Rcpp" = ADPPS_EYsubX_normal_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi)))
  )
)

test1 <- ADPPS_EYsubX_normal_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, inclusion = inclusion,
  eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EYsubX_normal_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, inclusion = inclusion,
  eta_initial = rep(0.1, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- ADPPS_EYsubX_normal_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, inclusion = inclusion,
  eta.initial = rep(0.1, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EYsubX_normal_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, inclusion = inclusion,
      eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EYsubX_normal_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, inclusion = inclusion,
      eta_initial = rep(0.1, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = ADPPS_EYsubX_normal_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, inclusion = inclusion,
      eta.initial = rep(0.1, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EYsubX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, inclusion = inclusion)

test2 <- ADPPS_EYsubX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, inclusion = inclusion)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EYsubX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, inclusion = inclusion),
    "Rcpp" = ADPPS_EYsubX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, inclusion = inclusion)
  )
)

test1 <- ADPPSvar_EYsubX_normal(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS.beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi)))

test2 <- ADPPSvar_EYsubX_normal_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
  phi = phi, PPS_beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPSvar_EYsubX_normal(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS.beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi))),
    "Rcpp" = ADPPSvar_EYsubX_normal_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
      phi = phi, PPS_beta = 0.1, inclusion = inclusion, eta = rep(0.1, length(phi)))
  )
)

##### Logistic regression model #####

p <- 2

theta0 <- 1
theta1 <- 1
theta2 <- 1

n <- 100
N <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rbinom(n = n,
            size = 1,
            prob = 1 / (1 + exp(-(theta0 + X1 * theta1 + X2 * theta2))))
test.data <- data.frame(response = Y,
                        covariate = X)

X1_shift <- rnorm(n = n, mean = 0, sd = 1)
X2_shift <- rbinom(n = n, size = 1, prob = 0.6)
X_shift <- cbind(X1_shift, X2_shift)
Y_shift <- rbinom(n = n,
                  size = 1,
                  prob = 1 / (1 + exp(-(theta0 + X1_shift * theta1 + X2_shift * theta2))))

### auxiliary information: average of X

### auxiliary information: average of Y

X <- as.matrix(X)
phi <- mean(Y_shift)

test1 <- AD_EY_logistic_Lagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, eta = 0.1)

test2 <- AD_EY_logistic_Lagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EY_logistic_Lagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, eta = 0.1),
    "Rcpp" = AD_EY_logistic_Lagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, eta = 0.1)
  )
)

test1 <- AD_EY_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EY_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, eta_initial = 0,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

test3 <- AD_EY_logistic_SolveLagrange_v1(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EY_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EY_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, eta_initial = 0,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7),
    "R_Rcpp" = AD_EY_logistic_SolveLagrange_v1(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)
  )
)

test1 <- AD_EY_logistic(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi)

test2 <- AD_EY_logistic_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EY_logistic(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi),
    "Rcpp" = AD_EY_logistic_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi)
  )
)

test1 <- ADvar_EY_logistic(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, eta = 0.1)

test2 <- ADvar_EY_logistic_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, eta = 0.1)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADvar_EY_logistic(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, eta = 0.1),
    "Rcpp" = ADvar_EY_logistic_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, eta = 0.1)
  )
)

















### auxiliary information: average of X

X <- as.matrix(X)
phi <- colMeans(X_shift)

test1 <- ADPPS_EX_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, PPS.beta = 0.1, eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EX_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, PPS_beta = 0.1, eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EX_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, PPS.beta = 0.1, eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EX_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, PPS_beta = 0.1, eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

### auxiliary information: average of Y

X <- as.matrix(X)
phi <- mean(Y_shift)



test1 <- ADPPS_EY_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, PPS.beta = 0.1, eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EY_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, PPS_beta = 0.1, eta_initial = 0,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EY_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, PPS.beta = 0.1, eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EY_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, PPS_beta = 0.1, eta_initial = 0,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EY_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, CS.beta = rep(0.1, p), eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EY_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, CS_beta = rep(0.1, p), eta_initial = 0,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EY_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, CS.beta = rep(0.1, p), eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EY_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, CS_beta = rep(0.1, p), eta_initial = 0,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

### auxiliary information: average of X given Y

X <- as.matrix(X)

phi <- matrix(0, 2, p)
phi[1, ] <- colMeans(X_shift[Y_shift == 1, ])
phi[2, ] <- colMeans(X_shift[Y_shift == 0, ])

test1 <- AD_EXsubY_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EXsubY_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EXsubY_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EXsubY_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EXsubY_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
  CS.beta = rep(0.1, p), eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EXsubY_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
  CS_beta = rep(0.1, p), eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EXsubY_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
      CS.beta = rep(0.1, p), eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EXsubY_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), phi = phi,
      CS_beta = rep(0.1, p), eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

### auxiliary information: average of Y given X

X <- as.matrix(X)

phi <- c(mean(Y_shift[X_shift[, 1] > 0]),
         mean(Y_shift[X_shift[, 1] <= 0]))
inclusion <- cbind(X_shift[, 1] > 0,
                   X_shift[, 1] <= 0)

test1 <- AD_EYsubX_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EYsubX_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, inclusion = inclusion,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EYsubX_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EYsubX_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, inclusion = inclusion,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EYsubX_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, PPS.beta = 0.1, inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EYsubX_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, PPS_beta = 0.1, inclusion = inclusion,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EYsubX_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, PPS.beta = 0.1, inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EYsubX_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, PPS_beta = 0.1, inclusion = inclusion,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EYsubX_logistic_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EYsubX_logistic_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2),
  phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EYsubX_logistic_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EYsubX_logistic_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2),
      phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)









##### Gamma regression model #####

p <- 2

theta0 <- 1
theta1 <- 1
theta2 <- 1
nu0 <- 2
beta0 <- -0.2

n <- 100
N <- 1000
N_sim <- 10000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rgamma(n = n, shape = nu0, rate = nu0 / exp(theta0 + X1 * theta1 + X2 * theta2))
test.data <- data.frame(response = Y,
                        covariate = X)

X1_sim <- rnorm(n = N_sim, mean = 0, sd = 1)
X2_sim <- rbinom(n = N_sim, size = 1, prob = 0.6)
X_sim <- cbind(X1_sim, X2_sim)
Y_sim <- rgamma(n = N_sim,
                shape = nu0,
                rate = nu0 / exp(theta0 + X1_sim * theta1 + X2_sim * theta2))

Y_shift <- sample(x = Y_sim, size = N, replace = TRUE, prob = exp(beta0 * Y_sim))
X_shift <- matrix(0, nrow = N, ncol = 2)
for (i in 1:N)
{
  w_i <- dgamma(Y_shift[i], shape = nu0,
                rate = nu0 / exp(theta0 + X1_sim * theta1 + X2_sim * theta2))
  X_shift[i, ] <- X_sim[sample(1:N_sim, size = 1, prob = w_i), ]
}
X1_shift <- X_shift[, 1]
X2_shift <- X_shift[, 2]

### auxiliary information: average of X

X <- as.matrix(X)
phi <- colMeans(X_shift)

test1 <- ADPPS_EX_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS.beta = beta0, eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EX_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS_beta = beta0, eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EX_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS.beta = beta0, eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EX_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS_beta = beta0, eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

### auxiliary information: average of Y

X <- as.matrix(X)
phi <- mean(Y_shift)

test1 <- AD_EY_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EY_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, eta_initial = 0,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EY_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EY_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, eta_initial = 0,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EY_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS.beta = beta0, eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EY_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS_beta = beta0, eta_initial = 0,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EY_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS.beta = beta0, eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EY_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS_beta = beta0, eta_initial = 0,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EY_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, CS.beta = rep(0.1, p), eta.initial = 0,
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EY_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, CS_beta = rep(0.1, p), eta_initial = 0,
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EY_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, CS.beta = rep(0.1, p), eta.initial = 0,
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EY_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, CS_beta = rep(0.1, p), eta_initial = 0,
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

### auxiliary information: average of X given Y

X <- as.matrix(X)

y.pts <- matrix(0, 2, 2)
y.pts[1, ] <- quantile(Y_shift, c(0.0, 0.5))
y.pts[2, ] <- quantile(Y_shift, c(0.5, 1.0))

phi <- matrix(0, 2, p)
phi[1, ] <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) & (Y_shift <= y.pts[1, 2]), ])
phi[2, ] <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) & (Y_shift <= y.pts[2, 2]), ])

test1 <- AD_EXsubY_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, y.pts = y.pts,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EXsubY_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, y_pts = y.pts,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EXsubY_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, y.pts = y.pts,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EXsubY_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, y_pts = y.pts,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EXsubY_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS.beta = beta0, y.pts = y.pts,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EXsubY_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS_beta = beta0, y_pts = y.pts,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EXsubY_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS.beta = beta0, y.pts = y.pts,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EXsubY_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS_beta = beta0, y_pts = y.pts,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EXsubY_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EXsubY_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EXsubY_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, CS.beta = rep(0.1, p), y.pts = y.pts,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EXsubY_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, CS_beta = rep(0.1, p), y_pts = y.pts,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

### auxiliary information: average of Y given X

X <- as.matrix(X)

phi <- c(mean(Y_shift[X_shift[, 1] > 0]),
         mean(Y_shift[X_shift[, 1] <= 0]))
inclusion <- cbind(X_shift[, 1] > 0,
               X_shift[, 1] <= 0)

test1 <- AD_EYsubX_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- AD_EYsubX_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, inclusion = inclusion,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = AD_EYsubX_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = AD_EYsubX_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, inclusion = inclusion,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADPPS_EYsubX_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS.beta = beta0, inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADPPS_EYsubX_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, PPS_beta = beta0, inclusion = inclusion,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADPPS_EYsubX_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS.beta = beta0, inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADPPS_EYsubX_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, PPS_beta = beta0, inclusion = inclusion,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

test1 <- ADCS_EYsubX_Gamma_SolveLagrange(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
  eta.initial = rep(0, length(phi)),
  iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7)

test2 <- ADCS_EYsubX_Gamma_SolveLagrange_rcpp(
  X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
  phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion,
  eta_initial = rep(0, length(phi)),
  iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ADCS_EYsubX_Gamma_SolveLagrange(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, CS.beta = rep(0.1, p), inclusion = inclusion,
      eta.initial = rep(0, length(phi)),
      iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-7),
    "Rcpp" = ADCS_EYsubX_Gamma_SolveLagrange_rcpp(
      X = X, alpha = theta0, beta = c(theta1, theta2), nu = nu0,
      phi = phi, CS_beta = rep(0.1, p), inclusion = inclusion,
      eta_initial = rep(0, length(phi)),
      iter_max = 10, step_rate = 2, step_max = 5, tol = 1e-7)
  )
)

