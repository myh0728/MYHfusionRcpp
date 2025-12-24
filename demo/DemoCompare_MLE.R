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

all.equal(as.vector(lL_normal(X = X, Y = Y,
                              alpha = theta0,
                              beta = c(theta1, theta2),
                              sigma = sigma0)),
          as.vector(lL_normal_rcpp(X = X, Y = Y,
                                   alpha = theta0,
                                   beta = c(theta1, theta2),
                                   sigma = sigma0)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_normal(X = X, Y = Y,
                    alpha = theta0,
                    beta = c(theta1, theta2),
                    sigma = sigma0),
    "Rcpp" = lL_normal_rcpp(X = X, Y = Y,
                            alpha = theta0,
                            beta = c(theta1, theta2),
                            sigma = sigma0)
  )
)

all.equal(as.vector(lL_normal(X = X, Y = Y,
                              alpha = theta0,
                              beta = c(theta1, theta2),
                              sigma = sigma0,
                              wi.boot = W)),
          as.vector(lL_normal_rcpp(X = X, Y = Y,
                                   alpha = theta0,
                                   beta = c(theta1, theta2),
                                   sigma = sigma0,
                                   w = W)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_normal(X = X, Y = Y,
                    alpha = theta0,
                    beta = c(theta1, theta2),
                    sigma = sigma0,
                    wi.boot = W),
    "Rcpp" = lL_normal_rcpp(X = X, Y = Y,
                            alpha = theta0,
                            beta = c(theta1, theta2),
                            sigma = sigma0,
                            w = W)
  )
)

res_R <- diff_lL_normal(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
res_Rcpp <- diff_lL_normal_rcpp(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_normal(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0),
    "Rcpp" = diff_lL_normal_rcpp(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
  )
)

res_R <- diff_lL_normal(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0, wi.boot = W)
res_Rcpp <- diff_lL_normal_rcpp(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0, w = W)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_normal(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0, wi.boot = W),
    "Rcpp" = diff_lL_normal_rcpp(X = X, Y = Y, alpha = theta0, beta = c(theta1, theta2), sigma = sigma0, w = W)
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

all.equal(as.vector(lL_logistic(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2))),
          as.vector(lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = theta0, beta = c(theta1, theta2))))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_logistic(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2)),
    "Rcpp" = lL_logistic_rcpp(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2))
  )
)

all.equal(as.vector(lL_logistic(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), wi.boot = W)),
          as.vector(lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = theta0, beta = c(theta1, theta2), w = W)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_logistic(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2), wi.boot = W),
    "Rcpp" = lL_logistic_rcpp(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2), w = W)
  )
)

res_R <- diff_lL_logistic(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2))
res_Rcpp <- diff_lL_logistic_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2))
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_logistic(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2)),
    "Rcpp" = diff_lL_logistic_rcpp(X = X, Y = Y,
                                   alpha = theta0, beta = c(theta1, theta2))
  )
)

res_R <- diff_lL_logistic(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2), wi.boot = W)
res_Rcpp <- diff_lL_logistic_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2), w = W)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_logistic(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2), wi.boot = W),
    "Rcpp" = diff_lL_logistic_rcpp(X = X, Y = Y,
                                   alpha = theta0, beta = c(theta1, theta2), w = W)
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

all.equal(as.vector(lL_Gamma(X = X, Y = Y,
                             alpha = theta0, beta = c(theta1, theta2), nu = nu)),
          as.vector(lL_Gamma_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2), nu = nu)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Gamma(X = X, Y = Y,
                   alpha = theta0, beta = c(theta1, theta2), nu = nu),
    "Rcpp" = lL_Gamma_rcpp(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2), nu = nu)
  )
)

all.equal(as.vector(lL_Gamma(X = X, Y = Y,
                             alpha = theta0, beta = c(theta1, theta2), nu = nu, wi.boot = W)),
          as.vector(lL_Gamma_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2), nu = nu, w = W)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Gamma(X = X, Y = Y,
                   alpha = theta0, beta = c(theta1, theta2), nu = nu, wi.boot = W),
    "Rcpp" = lL_Gamma_rcpp(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2), nu = nu, w = W)
  )
)

res_R <- diff_lL_Gamma(X = X, Y = Y,
                       alpha = theta0, beta = c(theta1, theta2), nu = nu)
res_Rcpp <- diff_lL_Gamma_rcpp(X = X, Y = Y,
                               alpha = theta0, beta = c(theta1, theta2), nu = nu)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_Gamma(X = X, Y = Y,
                        alpha = theta0, beta = c(theta1, theta2), nu = nu),
    "Rcpp" = diff_lL_Gamma_rcpp(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), nu = nu)
  )
)

res_R <- diff_lL_Gamma(X = X, Y = Y,
                       alpha = theta0, beta = c(theta1, theta2), nu = nu, wi.boot = W)
res_Rcpp <- diff_lL_Gamma_rcpp(X = X, Y = Y,
                               alpha = theta0, beta = c(theta1, theta2), nu = nu, w = W)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_Gamma(X = X, Y = Y,
                        alpha = theta0, beta = c(theta1, theta2), nu = nu, wi.boot = W),
    "Rcpp" = diff_lL_Gamma_rcpp(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), nu = nu, w = W)
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

all.equal(as.vector(lL_binomial(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), N = N)),
          as.vector(lL_binomial_rcpp(X = X, Y = Y,
                                     alpha = theta0, beta = c(theta1, theta2), N = N)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_binomial(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2), N = N),
    "Rcpp" = lL_binomial_rcpp(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2), N = N)
  )
)

all.equal(as.vector(lL_binomial(X = X, Y = Y,
                                alpha = theta0, beta = c(theta1, theta2), N = N, wi.boot = W)),
          as.vector(lL_binomial_rcpp(X = X, Y = Y,
                                     alpha = theta0, beta = c(theta1, theta2), N = N, w = W)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_binomial(X = X, Y = Y,
                      alpha = theta0, beta = c(theta1, theta2), N = N, wi.boot = W),
    "Rcpp" = lL_binomial_rcpp(X = X, Y = Y,
                              alpha = theta0, beta = c(theta1, theta2), N = N, w = W)
  )
)

res_R <- diff_lL_binomial(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2), N = N)
res_Rcpp <- diff_lL_binomial_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2), N = N)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_binomial(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2), N = N),
    "Rcpp" = diff_lL_binomial_rcpp(X = X, Y = Y,
                                   alpha = theta0, beta = c(theta1, theta2), N = N)
  )
)

res_R <- diff_lL_binomial(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2), N = N, wi.boot = W)
res_Rcpp <- diff_lL_binomial_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2), N = N, w = W)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_binomial(X = X, Y = Y,
                           alpha = theta0, beta = c(theta1, theta2), N = N, wi.boot = W),
    "Rcpp" = diff_lL_binomial_rcpp(X = X, Y = Y,
                                   alpha = theta0, beta = c(theta1, theta2), N = N, w = W)
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

all.equal(as.vector(lL_Poisson(X = X, Y = Y,
                               alpha = theta0, beta = c(theta1, theta2))),
          as.vector(lL_Poisson_rcpp(X = X, Y = Y,
                                    alpha = theta0, beta = c(theta1, theta2))))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Poisson(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2)),
    "Rcpp" = lL_Poisson_rcpp(X = X, Y = Y,
                             alpha = theta0, beta = c(theta1, theta2))
  )
)

all.equal(as.vector(lL_Poisson(X = X, Y = Y,
                               alpha = theta0, beta = c(theta1, theta2), wi.boot = W)),
          as.vector(lL_Poisson_rcpp(X = X, Y = Y,
                                    alpha = theta0, beta = c(theta1, theta2), w = W)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = lL_Poisson(X = X, Y = Y,
                     alpha = theta0, beta = c(theta1, theta2), wi.boot = W),
    "Rcpp" = lL_Poisson_rcpp(X = X, Y = Y,
                             alpha = theta0, beta = c(theta1, theta2), w = W)
  )
)

res_R <- diff_lL_Poisson(X = X, Y = Y,
                         alpha = theta0, beta = c(theta1, theta2))
res_Rcpp <- diff_lL_Poisson_rcpp(X = X, Y = Y,
                                 alpha = theta0, beta = c(theta1, theta2))
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_Poisson(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2)),
    "Rcpp" = diff_lL_Poisson_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2))
  )
)

res_R <- diff_lL_Poisson(X = X, Y = Y,
                         alpha = theta0, beta = c(theta1, theta2), wi.boot = W)
res_Rcpp <- diff_lL_Poisson_rcpp(X = X, Y = Y,
                                 alpha = theta0, beta = c(theta1, theta2), w = W)
check_grad <- all.equal(as.vector(res_R$gradient), as.vector(res_Rcpp$gradient))
check_hess <- all.equal(as.matrix(res_R$hessian), as.matrix(res_Rcpp$hessian), check.attributes = FALSE)
cat("Gradient Consistent:", isTRUE(check_grad), "\n")
if(!isTRUE(check_grad)) print(check_grad)
cat("Hessian Consistent:", isTRUE(check_hess), "\n")
if(!isTRUE(check_hess)) print(check_hess)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = diff_lL_Poisson(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2), wi.boot = W),
    "Rcpp" = diff_lL_Poisson_rcpp(X = X, Y = Y,
                                  alpha = theta0, beta = c(theta1, theta2), w = W)
  )
)





