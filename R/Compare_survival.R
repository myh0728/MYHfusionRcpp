KME_R_outer <- function(time.last, is.event, time.event)
{
  time.last <- as.vector(time.last)
  is.event <- as.vector(is.event)
  time.event <- as.vector(time.event)

  dNit <- outer(time.last, time.event, FUN = "==") * is.event
  Yit <- outer(time.last, time.event, FUN = ">=")
  Nhat <- colSums(dNit)
  Dhat <- colSums(Yit)
  dLhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))
  Shat <- cumprod(1 - dLhat)
  Lhat <- cumsum(dLhat)

  results <- list(jumps = time.event,
                  survival = Shat,
                  hazard = dLhat,
                  cumhazard = Lhat)

  return(results)
}

SKME_R_outer <- function(time.last, is.event, time.event,
                       X, x, K, h)
{
  time.last <- as.vector(time.last)
  is.event <- as.vector(is.event)
  time.event <- as.vector(time.event)

  number_n <- length(time.last)
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  dNit <- outer(time.last, time.event, FUN = "==") * is.event
  Yit <- outer(time.last, time.event, FUN = ">=")

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h), c(2, 3), prod)

  Nhat <- t(Kik) %*% dNit
  Dhat <- t(Kik) %*% Yit
  dLhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))
  dLhat <- pmin(pmax(dLhat, 0), 1)
  Shat <- t(apply(1 - dLhat, 1, cumprod))
  Lhat <- t(apply(dLhat, 1, cumsum))

  results <- list(jumps = time.event,
                  survival = Shat,
                  hazard = dLhat,
                  cumhazard = Lhat,
                  Dhat = Dhat,
                  Nhat = Nhat)

  return(results)
}
