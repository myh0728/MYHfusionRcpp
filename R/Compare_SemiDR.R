cumuSIR_R <- function(X, Y, eps.inv)
{
  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  Y.CP <- outer_leq_rcpp(as.vector(Y), as.vector(Y))
  varX.eigen <- eigen_rcpp(var(X))
  varX.eigen.value <- as.vector(varX.eigen$value)
  varX.eigen.value[varX.eigen.value < eps.inv] <- eps.inv
  normalizing <- varX.eigen$vector %*% diag(varX.eigen.value ^ (-0.5)) %*% t(varX.eigen$vector)
  X.cs <- t(normalizing %*% (t(X) - colMeans(X)))

  m.y <- t(X.cs) %*% Y.CP / number_n
  Km <- m.y %*% t(m.y) / number_n

  RR <- eigen_rcpp(Km)
  Bhat <- normalizing %*% RR$vector
  dimnames(Bhat) <- list(paste("covariate", 1:number_p, sep=""),
                         paste("direction", 1:number_p, sep=""))
  results <- list(basis = Bhat,
                  values = as.vector(RR$value))

  return(results)
}

cumuSAVE_R <- function(X, Y, eps.inv)
{
  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  Y.CP <- outer_leq_rcpp(as.vector(Y), as.vector(Y))
  varX.eigen <- eigen_rcpp(var(X))
  varX.eigen.value <- as.vector(varX.eigen$value)
  varX.eigen.value[varX.eigen.value < eps.inv] <- eps.inv
  normalizing <- varX.eigen$vector %*% diag(varX.eigen.value ^ (-0.5)) %*% t(varX.eigen$vector)
  X.cs <- t(normalizing %*% (t(X) - colMeans(X)))
  Y.CP.cs <- t(t(Y.CP) - colMeans(Y.CP))

  m.y <- t(X.cs) %*% Y.CP / number_n
  M.y <- m.y[rep(1:number_p, times = number_p), ] *
    m.y[rep(1:number_p, each = number_p), ]
  dim(M.y) <- c(number_p, number_p, number_n)
  N.y <- t(X.cs[, rep(1:number_p, times = number_p)] *
             X.cs[, rep(1:number_p, each = number_p)]) %*%
    Y.CP.cs / number_n
  dim(N.y) <- c(number_p, number_p, number_n)
  D.y <- N.y - M.y
  Km <- apply(apply(
    aperm(D.y[rep(1:number_p, times = number_p), , ], c(2, 1, 3)) *
      D.y[, rep(1:number_p, each = number_p), ], c(2, 3), sum
  ), 1, sum) / number_n
  dim(Km) <- c(number_p, number_p)

  RR <- eigen_rcpp(Km)
  Bhat <- normalizing %*% RR$vector
  dimnames(Bhat) <- list(paste("covariate", 1:number_p, sep=""),
                         paste("direction", 1:number_p, sep=""))
  results <- list(basis = Bhat,
                  values = as.vector(RR$value))
}
