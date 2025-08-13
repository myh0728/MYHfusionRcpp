cumuSIR_R <- function(X, Y, eps.inv)
{
  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  Y.CP <- ctingP_uni_rcpp(as.vector(Y), as.vector(Y))
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

