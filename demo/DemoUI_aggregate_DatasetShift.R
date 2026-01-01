rm(list = ls())

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
Y <- rnorm(n = n,
           mean = theta0 + X1 * theta1 + X2 * theta2,
           sd = sigma0)
test.data <- data.frame(response = Y,
                        covariate = X)

X1_shift <- rnorm(n = N, mean = 0, sd = 1)
X2_shift <- rbinom(n = N, size = 1, prob = 0.6)
X_shift <- cbind(X1_shift, X2_shift)
Y_shift <- rnorm(n = N,
                 mean = theta0 + X1_shift * theta1 + X2_shift * theta2,
                 sd = sigma0)

### auxiliary information

phi.EX <- colMeans(X_shift)
phi.EY <- mean(Y_shift)

y.pts <- matrix(0, 2, 2)
y.pts[1, ] <- quantile(Y_shift, c(0.0, 0.5))
y.pts[2, ] <- quantile(Y_shift, c(0.5, 1.0))

phi.EXsubY <- matrix(0, 2, p)
phi.EXsubY[1, ] <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) &
                                      (Y_shift <= y.pts[1, 2]), ])
phi.EXsubY[2, ] <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) &
                                      (Y_shift <= y.pts[2, 2]), ])

phi.EYsubX.1 <- mean(Y_shift[X_shift[, 1] > 0])
phi.EYsubX.2 <- mean(Y_shift[X_shift[, 2] > 0])
phi.EYsubX <- c(phi.EYsubX.1, phi.EYsubX.2)
inclusion <- cbind(X[, 1] > 0, X[, 2] > 0)

#############################################################################

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EY = list(phi = phi.EY,
                 ext.size = N,
                 ext.var = 10),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N,
                     ext.var = 10 * diag(4)),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N,
                     ext.var = 10 * diag(2)))










