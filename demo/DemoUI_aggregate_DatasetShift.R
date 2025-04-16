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

### combining auxiliary information

data = NULL
X = test.data[, paste("covariate.X", 1:p, sep = "")]
Y = test.data$response
shift = "NS"
distribution = "normal"
method = "EL"
do.SE = TRUE
initial = NULL
initial.DRM = NULL
info.EX = list(phi = phi.EX,
               ext.size = N,
               ext.var = NULL)
info.EY = list(phi = phi.EY,
               ext.size = N,
               ext.var = NULL)
info.EXsubY = list(phi = phi.EXsubY,
                   y.pts = y.pts,
                   ext.size = N,
                   ext.var = NULL)
info.EYsubX = list(phi = phi.EYsubX,
                   inclusion = inclusion,
                   ext.size = N,
                   ext.var = NULL)
iter.max = 10
step.rate = 2
step.max = 5
tol = 1e-5
eps.inv = 1e-7

###

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EY = list(phi = phi.EY))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EY = list(phi = phi.EY,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX),
  info.EY = list(phi = phi.EY),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N),
  info.EY = list(phi = phi.EY,
                 ext.size = N),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EY = list(phi = phi.EY))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EY = list(phi = phi.EY,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX),
  info.EY = list(phi = phi.EY),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N),
  info.EY = list(phi = phi.EY,
                 ext.size = N),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EY = list(phi = phi.EY))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EY = list(phi = phi.EY,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX),
  info.EY = list(phi = phi.EY),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "EL", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N),
  info.EY = list(phi = phi.EY,
                 ext.size = N),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EY = list(phi = phi.EY))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EY = list(phi = phi.EY,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX),
  info.EY = list(phi = phi.EY),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "NS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N),
  info.EY = list(phi = phi.EY,
                 ext.size = N),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EY = list(phi = phi.EY))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EY = list(phi = phi.EY,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX),
  info.EY = list(phi = phi.EY),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "CS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N),
  info.EY = list(phi = phi.EY,
                 ext.size = N),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EY = list(phi = phi.EY))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EY = list(phi = phi.EY,
                 ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX),
  info.EY = list(phi = phi.EY),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion))

GLMcombineAD.DatasetShift(
  X = test.data[, paste("covariate.X", 1:p, sep = "")],
  Y = test.data$response,
  shift = "PPS", distribution = "normal",
  method = "fast", do.SE = TRUE,
  info.EX = list(phi = phi.EX,
                 ext.size = N),
  info.EY = list(phi = phi.EY,
                 ext.size = N),
  info.EXsubY = list(phi = phi.EXsubY,
                     y.pts = y.pts,
                     ext.size = N),
  info.EYsubX = list(phi = phi.EYsubX,
                     inclusion = inclusion,
                     ext.size = N))















CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "normal",
                       method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "normal",
                       ext.sample.size = N, method = "EL",
                       info.EY = list(phi = mean(Y_shift)),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "normal",
                       method = "fast",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "normal",
                       ext.sample.size = N, method = "fast",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "normal",
                       method = "EL",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "normal",
                       ext.sample.size = N, method = "EL",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "normal",
                       method = "fast",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "normal",
                       ext.sample.size = N, method = "fast",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "normal",
                       method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "normal",
                       ext.sample.size = N, method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "normal",
                       method = "fast",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "normal",
                       ext.sample.size = N, method = "fast",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

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

### auxiliary information

phi.EX <- colMeans(X_shift)
phi.EY <- mean(Y_shift)

y.pts <- matrix(0, 2, 2)
y.pts[1, ] <- quantile(Y_shift, c(0.0, 0.5))
y.pts[2, ] <- quantile(Y_shift, c(0.5, 1.0))

phi.EXsubY <- matrix(0, 2, p)
phi.EXsubY[1, ] <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) & (Y_shift <= y.pts[1, 2]), ])
phi.EXsubY[2, ] <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) & (Y_shift <= y.pts[2, 2]), ])

phi.EYsubX.1 <- mean(Y_shift[X_shift[, 1] > 0])
phi.EYsubX.2 <- mean(Y_shift[X_shift[, 1] <= 0])
phi.EYsubX <- c(phi.EYsubX.1, phi.EYsubX.2)
inclusion <- cbind(X[, 1] > 0, X[, 1] <= 0)

### combining auxiliary information

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "Gamma",
                       method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "Gamma",
                       ext.sample.size = N, method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "Gamma",
                       method = "EL",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "Gamma",
                       ext.sample.size = N, method = "EL",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "Gamma",
                       method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "Gamma",
                       ext.sample.size = N, method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY,
                                          y.pts = y.pts),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

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
Y <- rbinom(n = n, size = 1,
            prob = 1 / (1 + exp(-(theta0 + X1 * theta1 + X2 * theta2))))
test.data <- data.frame(response = Y,
                        covariate = X)

X1_shift <- rnorm(n = N, mean = 0, sd = 1)
X2_shift <- rbinom(n = N, size = 1, prob = 0.6)
X_shift <- cbind(X1_shift, X2_shift)
Y_shift <- rbinom(n = N, size = 1,
                  prob = 1 / (1 + exp(-(theta0 + X1_shift * theta1 + X2_shift * theta2))))

### auxiliary information

phi.EX <- colMeans(X_shift)
phi.EY <- mean(Y_shift)

phi.EXsubY <- matrix(0, 2, p)
phi.EXsubY[1, ] <- colMeans(X_shift[Y_shift == 1, ])
phi.EXsubY[2, ] <- colMeans(X_shift[Y_shift == 0, ])

phi.EYsubX.1 <- mean(Y_shift[X_shift[, 1] > 0])
phi.EYsubX.2 <- mean(Y_shift[X_shift[, 1] <= 0])
phi.EYsubX <- c(phi.EYsubX.1, phi.EYsubX.2)
inclusion <- cbind(X[, 1] > 0, X[, 1] <= 0)

### combining auxiliary information

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "Bernoulli",
                       method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "NS", distribution = "Bernoulli",
                       ext.sample.size = N, method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "Bernoulli",
                       method = "EL",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "PPS", distribution = "Bernoulli",
                       ext.sample.size = N, method = "EL",
                       info.EX = list(phi = phi.EX),
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "Bernoulli",
                       method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))

CombineAD.DatasetShift(X = test.data[, paste("covariate.X", 1:p, sep = "")],
                       Y = test.data$response,
                       shift = "CS", distribution = "Bernoulli",
                       ext.sample.size = N, method = "EL",
                       info.EY = list(phi = phi.EY),
                       info.EXsubY = list(phi = phi.EXsubY),
                       info.EYsubX = list(phi = phi.EYsubX,
                                          inclusion = inclusion))




