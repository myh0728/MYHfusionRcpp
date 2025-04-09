##### Normal regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1
sigma0 <- 0.5

n <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rnorm(n = n, mean = theta0 + X1 * theta1 + X2 * theta2, sd = sigma0)
W <- rexp(n = n, rate = 1)

test1 <- lL_normal(X = X, Y = Y,
                   alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
test2 <- lL_normal_rcpp(X = X, Y = Y,
                        alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
test3 <- lL_normal_rcpp_fun(X = X, Y = Y,
                            alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)

abs(test1 - test2)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_normal(X = X, Y = Y,
                    alpha = theta0, beta = c(theta1, theta2), sigma = sigma0),
    "Rcpp" = lL_normal_rcpp(X = X, Y = Y,
                            alpha = theta0, beta = c(theta1, theta2), sigma = sigma0),
    "Rcpp_fun" = lL_normal_rcpp_fun(X = X, Y = Y,
                                    alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
  )
)

test1 <- lL_normal(X = X, Y = Y,
                   alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
                   wi.boot = W)
test2 <- lL_normal_w_rcpp(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
                          w = W)
test3 <- lL_normal_w_rcpp_fun(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
                              w = W)
abs(test1 - test2)
abs(test1 - test3)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_normal(X = X, Y = Y,
                    alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
                    wi.boot = W),
    "Rcpp" = lL_normal_w_rcpp(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
                              w = W),
    "Rcpp_fun" = lL_normal_w_rcpp_fun(X = X, Y = Y,
                                      alpha = theta0, beta = c(theta1, theta2), sigma = sigma0,
                                      w = W)
  )
)

test1 <- diff_lL_normal(X = X, Y = Y,
                        alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
test2 <- diff_lL_normal_rcpp(X = X, Y = Y,
                             alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
test3 <- diff_lL_normal_rcpp_v1(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
sum(abs(test1$gradient - test2$gradient))
sum(abs(test1$hessian - test2$hessian))
sum(abs(test1$gradient - test3$gradient))
sum(abs(test1$hessian - test3$hessian))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_normal(X = X, Y = Y,
                         alpha = theta0, beta = c(theta1, theta2), sigma = sigma0),
    "Rcpp" = diff_lL_normal_rcpp(X = X, Y = Y,
                                 alpha = theta0, beta = c(theta1, theta2), sigma = sigma0),
    "Rcpp_v1" = diff_lL_normal_rcpp_v1(X = X, Y = Y,
                                       alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
  )
)

##### Logistic regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1

n <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rbinom(n = n, size = 1, prob = 1 / (1 + exp(-(theta0 + X1 * theta1 + X2 * theta2))))
W <- rexp(n = n, rate = 1)

test1 <- lL_logistic(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2))
test2 <- lL_logistic_rcpp(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2))
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_logistic(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2)),
    "Rcpp" = lL_logistic_rcpp(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2))
  )
)

test1 <- lL_logistic(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2),
                     wi.boot = W)
test2 <- lL_logistic_w_rcpp(X = X, Y = Y,
                            alpha = theta0, beta = c(theta1, theta2),
                            w = W)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_logistic(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2),
                      wi.boot = W),
    "Rcpp" = lL_logistic_w_rcpp(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2),
                                w = W)
  )
)

test1 <- diff_lL_logistic(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2))
test2 <- diff_lL_logistic_rcpp(X = X, Y = Y,
                               alpha = theta0, beta = c(theta1, theta2))
sum(abs(test1$gradient - test2$gradient))
sum(abs(test1$hessian - test2$hessian))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_logistic(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2)),
    "Rcpp" = diff_lL_logistic_rcpp(X = X, Y = Y,
                                   alpha = theta0, beta = c(theta1, theta2))
  )
)

##### Gamma regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1
nu <- 2

n <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rgamma(n = n, shape = nu, rate = nu * exp(-(theta0 + X1 * theta1 + X2 * theta2)))
W <- rexp(n = n, rate = 1)

test1 <- lL_Gamma(X = X, Y = Y,
                  alpha = theta0, beta = c(theta1, theta2), nu = nu)
test2 <- lL_Gamma_rcpp(X = X, Y = Y,
                       alpha = theta0, beta = c(theta1, theta2), nu = nu)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Gamma(X = X, Y = Y,
                   alpha = theta0, beta = c(theta1, theta2), nu = nu),
    "Rcpp" = lL_Gamma_rcpp(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2), nu = nu)
  )
)

test1 <- lL_Gamma(X = X, Y = Y,
                  alpha = theta0, beta = c(theta1, theta2), nu = nu,
                  wi.boot = W)
test2 <- lL_Gamma_w_rcpp(X = X, Y = Y,
                         alpha = theta0, beta = c(theta1, theta2), nu = nu,
                         w = W)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Gamma(X = X, Y = Y,
                   alpha = theta0, beta = c(theta1, theta2), nu = nu,
                   wi.boot = W),
    "Rcpp" = lL_Gamma_w_rcpp(X = X, Y = Y,
                             alpha = theta0, beta = c(theta1, theta2), nu = nu,
                             w = W)
  )
)

test1 <- diff_lL_Gamma(X = X, Y = Y,
                       alpha = theta0, beta = c(theta1, theta2), nu = nu)
test2 <- diff_lL_Gamma_rcpp(X = X, Y = Y,
                            alpha = theta0, beta = c(theta1, theta2), nu = nu)
sum(abs(test1$gradient - test2$gradient))
sum(abs(test1$hessian - test2$hessian))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_Gamma(X = X, Y = Y,
                        alpha = theta0, beta = c(theta1, theta2), nu = nu),
    "Rcpp" = diff_lL_Gamma_rcpp(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), nu = nu)
  )
)

##### Binomial regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1
N <- 5

n <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rbinom(n = n, size = N, prob = 1 / (1 + exp(-(theta0 + X1 * theta1 + X2 * theta2))))
W <- rexp(n = n, rate = 1)

test1 <- lL_binomial(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2), N = N)
test2 <- lL_binomial_rcpp(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2), N = N)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_binomial(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2), N = N),
    "Rcpp" = lL_binomial_rcpp(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2), N = N)
  )
)

test1 <- lL_binomial(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2), N = N,
                     wi.boot = W)
test2 <- lL_binomial_w_rcpp(X = X, Y = Y,
                            alpha = theta0, beta = c(theta1, theta2), N = N,
                            w = W)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_binomial(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2), N = N,
                      wi.boot = W),
    "Rcpp" = lL_binomial_w_rcpp(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), N = N,
                                w = W)
  )
)

##### Poisson regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1

n <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rpois(n = n, lambda = exp(theta0 + X1 * theta1 + X2 * theta2))
W <- rexp(n = n, rate = 1)

test1 <- lL_Poisson(X = X, Y = Y,
                    alpha = theta0, beta = c(theta1, theta2))
test2 <- lL_Poisson_rcpp(X = X, Y = Y,
                         alpha = theta0, beta = c(theta1, theta2))
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Poisson(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2)),
    "Rcpp" = lL_Poisson_rcpp(X = X, Y = Y,
                             alpha = theta0, beta = c(theta1, theta2))
  )
)

test1 <- lL_Poisson(X = X, Y = Y,
                    alpha = theta0, beta = c(theta1, theta2),
                    wi.boot = W)
test2 <- lL_Poisson_w_rcpp(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2),
                           w = W)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Poisson(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2),
                     wi.boot = W),
    "Rcpp" = lL_Poisson_w_rcpp(X = X, Y = Y,
                               alpha = theta0, beta = c(theta1, theta2),
                               w = W)
  )
)
