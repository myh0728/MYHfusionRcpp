rm(list = ls())

##### Normal regression model #####

n <- 100
p <- 3
alpha0 <- 1
beta0 <- rep(0.2, p)
sigma0 <- 0.1

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)
test.data <- simGLM.normal(Xi = X,
                           alpha0 = alpha0,
                           beta0 = beta0,
                           sigma0 = sigma0,
                           seed = 123)

test1 <- GLM.MLE(data = test.data,
                 X.name = paste("covariate", 1:p, sep="."),
                 Y.name = "response",
                 distribution = "normal",
                 initial = c(-1, rep(1, p), 1),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

test2 <- GLM.MLE(X = test.data[, paste("covariate", 1:p, sep=".")],
                 Y = test.data[, "response"],
                 distribution = "normal",
                 initial = c(-1, rep(1, p), 1),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

##### Logistic regression model

n <- 100
p <- 3
alpha0 <- 0.1
beta0 <- rep(0.2, p)

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)

test.data <- simGLM.logistic(Xi = X,
                             alpha0 = alpha0,
                             beta0 = beta0,
                             seed = 123)

test1 <- GLM.MLE(data = test.data,
                 X.name = paste("covariate", 1:p, sep="."),
                 Y.name = "response",
                 distribution = "Bernoulli",
                 initial = c(-1, rep(1, p)),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

test2 <- GLM.MLE(X = test.data[, paste("covariate", 1:p, sep=".")],
                 Y = test.data[, "response"],
                 distribution = "Bernoulli",
                 initial = c(-1, rep(1, p)),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

##### Gamma regression model

n <- 100
p <- 3
alpha0 <- 1
beta0 <- rep(0.2, p)
nu0 <- 2

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)
test.data <- simGLM.Gamma(Xi = X,
                          alpha0 = alpha0,
                          beta0 = beta0,
                          nu0 = nu0,
                          seed = 123)

test1 <- GLM.MLE(data = test.data,
                 X.name = paste("covariate", 1:p, sep="."),
                 Y.name = "response",
                 distribution = "Gamma",
                 initial = c(-1, rep(1, p), 1),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

test2 <- GLM.MLE(X = test.data[, paste("covariate", 1:p, sep=".")],
                 Y = test.data[, "response"],
                 distribution = "Gamma",
                 initial = c(-1, rep(1, p), 1),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

##### Binomial regression model

n <- 100
p <- 3
alpha0 <- 0.1
beta0 <- rep(0.2, p)

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)

test.data <- simGLM.binomial(Xi = X, N.binomial = 10,
                             alpha0 = alpha0,
                             beta0 = beta0,
                             seed = 123)

test1 <- GLM.MLE(data = test.data,
                 X.name = paste("covariate", 1:p, sep="."),
                 Y.name = "response",
                 distribution = "binomial",
                 N.binomial = 10,
                 initial = c(-1, rep(1, p)),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

test2 <- GLM.MLE(X = test.data[, paste("covariate", 1:p, sep=".")],
                 Y = test.data[, "response"],
                 distribution = "binomial",
                 N.binomial = 10,
                 initial = c(-1, rep(1, p)),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

##### Poisson regression model

n <- 100
p <- 3
alpha0 <- 0.1
beta0 <- rep(0.2, p)

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)

test.data <- simGLM.Poisson(Xi = X,
                            alpha0 = alpha0,
                            beta0 = beta0,
                            seed = 123)

test1 <- GLM.MLE(data = test.data,
                 X.name = paste("covariate", 1:p, sep="."),
                 Y.name = "response",
                 distribution = "Poisson",
                 initial = c(-1, rep(1, p)),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))

test2 <- GLM.MLE(X = test.data[, paste("covariate", 1:p, sep=".")],
                 Y = test.data[, "response"],
                 distribution = "Poisson",
                 initial = c(-1, rep(1, p)),
                 do.SE = TRUE,
                 X.future = matrix(0, nrow = 1, ncol = p))





