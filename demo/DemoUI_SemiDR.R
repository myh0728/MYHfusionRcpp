n <- 100
p <- 5
k <- 10

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- matrix(rnorm(k * p), nrow = k, ncol = p)
w <- rexp(n = n, rate = 1)
beta0 <- c(1, rep(0.5, p - 1))
Y <- as.vector(sin(X %*% beta0)) + rnorm(n, mean = 0, sd = 0.2)
y <- seq(min(Y), max(Y), length.out = 10)

test <- SemiSI(X = X, Y = Y, type = "mean", kernel = "K2.Epanechnikov")

test <- SemiSI(X = X, Y = Y, type = "mean", kernel = "K2.Epanechnikov", wi.boot = w)

test <- SemiSI(X = X, Y = Y, type = "distribution", kernel = "K2.Epanechnikov")

test <- SemiSI(X = X, Y = Y, type = "distribution", kernel = "K2.Epanechnikov", wi.boot = w)

test <- SemiSI(X = X, Y = Y, type = "mean", kernel = "K2.Biweight")

test <- SemiSI(X = X, Y = Y, type = "mean", kernel = "K2.Biweight", wi.boot = w)

test <- SemiSI(X = X, Y = Y, type = "distribution", kernel = "K2.Biweight")

test <- SemiSI(X = X, Y = Y, type = "distribution", kernel = "K2.Biweight", wi.boot = w)

test <- SemiSI(X = X, Y = Y, type = "mean", kernel = "K4.Biweight")

test <- SemiSI(X = X, Y = Y, type = "mean", kernel = "K4.Biweight", wi.boot = w)

test <- SemiSI(X = X, Y = Y, type = "distribution", kernel = "K4.Biweight")

test <- SemiSI(X = X, Y = Y, type = "distribution", kernel = "K4.Biweight", wi.boot = w)

test <- SemiMI(X = X, Y = Y, type = "mean", dimension = 2, kernel = "K2.Epanechnikov")

test <- SemiMI(X = X, Y = Y, type = "mean", dimension = 2, kernel = "K2.Epanechnikov", wi.boot = w)

test <- SemiMI(X = X, Y = Y, type = "distribution", dimension = 2, kernel = "K2.Epanechnikov")

test <- SemiMI(X = X, Y = Y, type = "distribution", dimension = 2, kernel = "K2.Epanechnikov", wi.boot = w)

test <- SemiMI(X = X, Y = Y, type = "mean", dimension = 2, kernel = "K2.Biweight")

test <- SemiMI(X = X, Y = Y, type = "mean", dimension = 2, kernel = "K2.Biweight", wi.boot = w)

test <- SemiMI(X = X, Y = Y, type = "distribution", dimension = 2, kernel = "K2.Biweight")

test <- SemiMI(X = X, Y = Y, type = "distribution", dimension = 2, kernel = "K2.Biweight", wi.boot = w)

test <- SemiMI(X = X, Y = Y, type = "mean", dimension = 2, kernel = "K4.Biweight")

test <- SemiMI(X = X, Y = Y, type = "mean", dimension = 2, kernel = "K4.Biweight", wi.boot = w)

test <- SemiMI(X = X, Y = Y, type = "distribution", dimension = 2, kernel = "K4.Biweight")

test <- SemiMI(X = X, Y = Y, type = "distribution", dimension = 2, kernel = "K4.Biweight", wi.boot = w)



