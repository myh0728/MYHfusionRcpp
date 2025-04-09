K2_Ep <- function(u)
{
  3 / 4 * (1 - u ^ 2) * (abs(u) <= 1)
}

K2_Bw <- function(u)
{
  15 / 16 * (1 - u ^ 2) ^ 2 * (abs(u) <= 1)
}

K4_Bw <- function(u)
{
  (105 / 64) * (1 - 3 * (u ^ 2)) * ((1 - u ^ 2) ^ 2) * (abs(u) <= 1)
}

KDE_R <- function(X, x, K, h){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kik) / number_n

  return(Dhat)
}

KDE_w_R <- function(X, x, K, h, w){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kik * w) / number_n

  return(Dhat)
}
