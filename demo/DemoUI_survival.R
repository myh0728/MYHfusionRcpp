rm(list = ls())

n <- 100
p <- 2
beta0 <- c(1, -1)

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
C <- 1 + rexp(n = n, rate = 1)
data.test <- simSurv.PH(Xi = X, Ci = C,
                        beta0 = beta0,
                        cumhazard.inv = function(t){t})

surv.data.test <- SurvData_process(data = data.test,
                                   X.name = paste0("covariate.X", 1:p),
                                   Y.name = "time.last",
                                   D.name = "is.event")

surv.data.test <- SurvData_process(data = data.test,
                                   Y.name = "time.last",
                                   D.name = "is.event")

KaplanMeier(data = data.test,
            Y.name = "time.last",
            D.name = "is.event")

KaplanMeier(data = data.test,
            Y.name = "time.last",
            D.name = "is.event",
            time.points = seq(0, 2, 0.01))
