K2_Ep_R <- function(u)
{
  3 / 4 * (1 - u ^ 2) * (abs(u) <= 1)
}

K2_G_R <- function(u)
{
  dnorm(u)
}

K2_Bw_R <- function(u)
{
  15 / 16 * (1 - u ^ 2) ^ 2 * (abs(u) <= 1)
}

K4_Bw_R <- function(u)
{
  (105 / 64) * (1 - 3 * (u ^ 2)) * ((1 - u ^ 2) ^ 2) * (abs(u) <= 1)
}

KDE_R <- function(X, x, K, h, w = NULL) {

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  if (is.null(w)) {

    w <- 1
  }

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]

  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kik * w) / number_n

  return(Dhat)
}

KNW_R <- function(Y, X, x, K, h, w = NULL) {

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  if (is.null(w)) {

    w <- 1
  }

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)

  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kik * w)
  Nhat <- colSums(Kik * Y * w)
  mhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  return(mhat)
}

KNWcdf_R <- function(Y, y, X, x, K, h, w = NULL) {

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_l <- length(y)

  if (is.null(w)) {

    w <- 1
  }

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)


  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kik * w)
  Nhat <- t(Kik) %*% (outer_leq_rcpp(Y, y) * w)
  mhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))
  mhat <- mhat * (mhat > 0) * (mhat < 1) + (mhat >= 1)

  return(mhat)
}

CVKNW_R <- function(Y, X, K, h, w = NULL) {

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(w)) {

    w <- rep(1, number_n)
  }

  Xij <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xij) <- c(number_n, number_n, number_p)

  Kij <- apply(K(aperm(Xij, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  diag(Kij) <- 0
  Dhat <- colSums(Kij * w)
  Nhat <- colSums(Kij * Y * w)
  mhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  cv_value <- sum(((Y - mhat) ^ 2) * w) / sum(w)

  return(cv_value)
}

CVKNWcdf_R <- function(Y, X, K, h, w = NULL){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(w)) {

    w <- rep(1, number_n)
  }

  Yik <- outer_leq_rcpp(Y, Y)

  Xij <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xij) <- c(number_n, number_n, number_p)
  Kij <- apply(K(aperm(Xij, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  diag(Kij) <- 0
  Dhat <- colSums(Kij * w)
  Nhat <- Kij %*% (Yik * w)
  mhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))
  mhat <- mhat * (mhat > 0) * (mhat < 1) + (mhat >= 1)
  squared_error <- (Yik - mhat) ^ 2
  cv_value <- sum(colSums(squared_error * w) * w) / (sum(w) ^ 2)

  return(cv_value)
}

LSKNW_R <- function(Y, X, K, h, w = NULL){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(w)) {

    w <- rep(1, number_n)
  }

  Xij <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xij) <- c(number_n, number_n, number_p)
  Kij <- apply(K(aperm(Xij, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kij * w)
  Nhat <- colSums(Kij * Y * w)
  mhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))
  ls_value <- sum(((Y - mhat) ^ 2) * w) / sum(w)

  return(ls_value)
}

LSKNWcdf_R <- function(Y, X, K, h, w = NULL){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(w)) {

    w <- rep(1, number_n)
  }

  Yik <- outer_leq_rcpp(Y, Y)
  Xij <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xij) <- c(number_n, number_n, number_p)
  Kij <- apply(K(aperm(Xij, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kij * w)
  Nhat <- Kij %*% (Yik * w)
  mhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))
  mhat <- mhat * (mhat > 0) * (mhat < 1) + (mhat >= 1)
  ls_value <- sum(colSums(((Yik - mhat) ^ 2) * w) * w) / (sum(w) ^ 2)

  return(ls_value)
}








