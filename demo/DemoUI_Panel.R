library(abind)

n <- 200
L <- 12

beta0 <- c(1, -1)

SN <- 100
results.simulation <- matrix(0, nrow = SN, ncol = 5)
dimnames(results.simulation) <- list(paste0("sim", 1:SN),
                                     c("CV", "SS", "SSeff.CV", "SSeff.SS", "bandwidth"))

for (sn in 1:SN) {

  set.seed(12345+sn)

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

  esti <- Panel.SemiSI(X = X, Y = Y, kernel = "K2.Biweight",
                       beta.initial = beta0)

  results.simulation[sn, "CV"] <- esti$coef.CV[2]
  results.simulation[sn, "SS"] <- esti$coef.SS[2]
  results.simulation[sn, "SSeff.CV"] <- esti$coef.SSeff.CV[2]
  results.simulation[sn, "SSeff.SS"] <- esti$coef.SSeff.SS[2]
  results.simulation[sn, "bandwidth"] <- esti$bandwidth

  print(sn)
}

summary(results.simulation)
apply(results.simulation, 2, sd)


