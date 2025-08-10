n <- 100
p <- 3
k <- 10

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- matrix(rnorm(k * p), nrow = k, ncol = p)
w <- rexp(n = n, rate = 1)
beta0 <- c(1, rep(0.5, p - 1))
Y <- as.vector(sin(X %*% beta0)) + rnorm(n, mean = 0, sd = 0.2)
y <- seq(min(Y), max(Y), length.out = 10)

test <- KernelSmooth(X = X, Y = Y, x = x,
                     type = "mean", kernel = "K2.Epanechnikov", bandwidth = 0.5)

test <- KernelSmooth(X = X, Y = Y, x = x,
                     type = "mean", kernel = "K2.Epanechnikov")

test <- KernelSmooth(X = X, Y = Y, x = x,
                     type = "mean", kernel = "K2.Epanechnikov", bandwidth = 0.5, wi.boot = w)

test <- KernelSmooth(X = X, Y = Y, x = x,
                     type = "mean", kernel = "K2.Epanechnikov", wi.boot = w)

test <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                     type = "distribution", kernel = "K2.Epanechnikov", bandwidth = 0.5)

test <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                     type = "distribution", kernel = "K2.Epanechnikov", wi.boot = w)

test <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                     type = "distribution", kernel = "K2.Biweight")

test <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                     type = "distribution", kernel = "K4.Biweight")

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




