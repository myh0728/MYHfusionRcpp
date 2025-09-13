##### Gauss elimination to "row"-echelon forms #####

# M: a matrix, c(m, n) matrix

Gauss.row <- function(M)
{
  m <- nrow(M)
  n <- ncol(M)
  k.max <- min(m, n)
  n.pivot <- 0

  for (k in 1:n)
  {
    i.max <- which.max(
      abs(
        M[(n.pivot + 1):m, k]
      )
    )[1] + n.pivot

    if (M[i.max, k] == 0)
      next

    n.pivot <- n.pivot + 1
    M[c(n.pivot, i.max), ] <- M[c(i.max, n.pivot), ]


    for (i in 1:m)
    {
      if (i == n.pivot)
        next

      M[i, ] <- M[i, ] - M[n.pivot, ] * M[i, k] / M[n.pivot, k]
    }

    M[n.pivot, ] <- M[n.pivot, ] / M[n.pivot, k]

    if (n.pivot == k.max)
      break
  }

  return(M)
}

##### Inverse of a matrix

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
