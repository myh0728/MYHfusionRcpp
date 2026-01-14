rm(list = ls())

library(abind)

n <- 200
L <- 12

beta0 <- c(1, -1)

SN <- 10
results.simulation <- array(0, c(3, 5, SN))
dimnames(results.simulation) <- list(c("Frailty", "Exchangable", "AR1"),
                                     c("CV", "SS", "SSeff.CV", "SSeff.SS", "bandwidth"),
                                     paste0("sim", 1:SN))

link.evaluation <- seq(-2, 2, 0.01)
results.link <- array(0, c(length(link.evaluation), 3, 2, SN))
dimnames(results.link) <- list(link.evaluation,
                               c("Frailty", "Exchangable", "AR1"),
                               c("CV", "SS"),
                               paste0("sim", 1:SN))

for (sn in 1:SN) {

  set.seed(12345 + sn)

  Y <- X1 <- X2 <- matrix(0, nrow = n, ncol = L)
  Z <- rep(0, n)

  for (i in 1:n) {

    Z[i] <- runif(1, min = 0, max = 1)

    for (j in 1:L) {

      X1[i, j] <- rnorm(1, mean = 0, sd = 1)
      X2[i, j] <- rbinom(1, size = 1, prob = 0.5)
      SI <- sum(c(X1[i, j], X2[i, j]) * beta0)
      p_ij <- Z[i] / (1 + exp(-SI))
      Y[i, j] <- rbinom(1, size = 1, prob = p_ij)
    }
  }

  X <- abind(X1, X2, along = 3)

  esti <- Panel.SemiSI.Frailty(X = X, Y = Y, kernel = "K2.Biweight",
                               beta.initial = beta0,
                               link.eval = link.evaluation)

  results.simulation["Frailty", "CV", sn] <- esti$coef.CV[2]
  results.simulation["Frailty", "SS", sn] <- esti$coef.SS[2]
  results.simulation["Frailty", "SSeff.CV", sn] <- esti$coef.SSeff.CV[2]
  results.simulation["Frailty", "SSeff.SS", sn] <- esti$coef.SSeff.SS[2]
  results.simulation["Frailty", "bandwidth", sn] <- esti$bandwidth

  results.link[, "Frailty", "CV", sn] <- esti$link.CV
  results.link[, "Frailty", "SS", sn] <- esti$link.SS

  esti <- Panel.SemiSI.Exchangable(X = X, Y = Y, kernel = "K2.Biweight",
                                   beta.initial = beta0,
                                   link.eval = link.evaluation)

  results.simulation["Exchangable", "CV", sn] <- esti$coef.CV[2]
  results.simulation["Exchangable", "SS", sn] <- esti$coef.SS[2]
  results.simulation["Exchangable", "SSeff.CV", sn] <- esti$coef.SSeff.CV[2]
  results.simulation["Exchangable", "SSeff.SS", sn] <- esti$coef.SSeff.SS[2]
  results.simulation["Exchangable", "bandwidth", sn] <- esti$bandwidth

  results.link[, "Exchangable", "CV", sn] <- esti$link.CV
  results.link[, "Exchangable", "SS", sn] <- esti$link.SS

  esti <- Panel.SemiSI.AR1(X = X, Y = Y, kernel = "K2.Biweight",
                           beta.initial = beta0,
                           link.eval = link.evaluation)

  results.simulation["AR1", "CV", sn] <- esti$coef.CV[2]
  results.simulation["AR1", "SS", sn] <- esti$coef.SS[2]
  results.simulation["AR1", "SSeff.CV", sn] <- esti$coef.SSeff.CV[2]
  results.simulation["AR1", "SSeff.SS", sn] <- esti$coef.SSeff.SS[2]
  results.simulation["AR1", "bandwidth", sn] <- esti$bandwidth

  results.link[, "AR1", "CV", sn] <- esti$link.CV
  results.link[, "AR1", "SS", sn] <- esti$link.SS

  print(sn)
}

apply(results.simulation, c(1, 2), mean)
apply(results.simulation, c(1, 2), sd)



link.MC <- apply(results.link, c(1, 2, 3), mean)

par(mfrow = c(1,2))
plot(link.evaluation, link.MC[, "Frailty", "CV"], type = 'l')
lines(link.evaluation, 0.5 / (1 + exp(-link.evaluation)), col = 2)
plot(link.evaluation, link.MC[, "Frailty", "SS"], type = 'l')
lines(link.evaluation, 0.5 / (1 + exp(-link.evaluation)), col = 2)

