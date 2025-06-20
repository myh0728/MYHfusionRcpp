ginv.sym.eigen <- function(M, eps.inv)
{
  number_m <- dim(M)[1]
  eigen.M <- eigen_rcpp(M)
  eigen.M.value <- as.vector(eigen.M$value)
  eigen.M.vector <- eigen.M$vector
  eigen.M.value[abs(eigen.M.value) < eps.inv] <- 0
  eigen.M.value.inv <- (eigen.M.value != 0) /
    (eigen.M.value + (eigen.M.value == 0))
  M.ginv <- eigen.M.vector %*%
    diag(eigen.M.value.inv, nrow = number_m, ncol = number_m) %*%
    t(eigen.M.vector)

  return(M.ginv)
}
