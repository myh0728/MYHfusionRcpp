eXsq_R <- function(data_X){

  number_n <- dim(data_X)[1]
  number_p <- dim(data_X)[2]

  eXsq <- colMeans(
    data_X[, rep(1:number_p, times = number_p)] *
      data_X[, rep(1:number_p, each = number_p)]
  )

  return(eXsq)
}

eXsq_w_R <- function(data_X,
                     weight){

  number_n <- dim(data_X)[1]
  number_p <- dim(data_X)[2]

  eXsq <- colMeans(
    data_X[, rep(1:number_p, times = number_p)] *
      data_X[, rep(1:number_p, each = number_p)] * weight
  )

  return(eXsq)
}

Xsq_lowtri_R <- function(data_X){

  number_n <- dim(data_X)[1]
  number_p <- dim(data_X)[2]

  Xsq <- data_X[, rep(1:number_p, times = number_p)] *
    data_X[, rep(1:number_p, each = number_p)]

  A <- matrix(1:(number_p ^ 2), number_p, number_p)
  index.lowtri <- A[lower.tri(A, diag = TRUE)]
  Xsq_lowtri <- Xsq[, index.lowtri]

  return(Xsq_lowtri)
}

twoXYsym_lowtri_R <- function(data_X,
                              data_Y){

  number_n <- dim(data_X)[1]
  number_p <- dim(data_X)[2]

  XY <- data_X[, rep(1:number_p, times = number_p)] *
    data_Y[, rep(1:number_p, each = number_p)]
  YX <- data_Y[, rep(1:number_p, times = number_p)] *
    data_X[, rep(1:number_p, each = number_p)]
  twoXYsym <- XY + YX

  A <- matrix(1:(number_p ^ 2), number_p, number_p)
  index.lowtri <- A[lower.tri(A, diag = TRUE)]
  twoXYsym_lowtri <- twoXYsym[, index.lowtri]

  return(twoXYsym_lowtri)
}

ctingP_R <- function(Y, y){

  number_n <- dim(Y)[1]
  number_p <- dim(Y)[2]
  number_k <- dim(y)[1]

  CP <- matrix(apply(as.matrix(
    Y[rep(1:number_n, times = number_k), ] <=
      y[rep(1:number_k, each = number_n), ]
  ), 1, prod), nrow = number_n, ncol = number_k)

  return(CP)
}

ctingP_uni_R <- function(Y, y){

  number_n <- length(Y)
  number_k <- length(y)

  CP <- matrix(Y, nrow = number_n, ncol = number_k) <=
    matrix(y, nrow = number_n, ncol = number_k, byrow = TRUE)

  return(CP)
}

ctingP_uni_R_outer <- function(Y, y){

  CP <- outer(Y, y, FUN = "<=")

  return(CP)
}

countAinB_outer <- function(A, B){

  cts <- rowSums(outer(A, B, FUN = "=="))

  return(cts)
}

countAinB_W_outer <- function(A, B, W){

  cts <- colSums(outer(B, A, FUN = "==") * W)

  return(cts)
}

countAinB_for <- function(A, B){

  cts <- rep(0, length(A))

  for (i in 1:length(A)){

    for (j in 1:length(B)){

      if (A[i] == B[j]){

        cts[i] <- cts[i] + 1
      }
    }
  }

  return(cts)
}

countAinB_W_for <- function(A, B, W){

  cts <- rep(0, length(A))

  for (i in 1:length(A)){

    for (j in 1:length(B)){

      if (A[i] == B[j]){

        cts[i] <- cts[i] + W[j]
      }
    }
  }

  return(cts)
}

rankAinB_outer <- function(A, B){

  rk <- rowSums(outer(A, B, FUN = ">="))

  return(rk)
}

rankAinB_W_outer <- function(A, B, W){

  rk <- colSums(outer(B, A, "<=") * W)
}

rankAinB_for <- function(A, B){

  rk <- rep(0, length(A))

  for (i in 1:length(A)){

    for (j in 1:length(B)){

      if (A[i] >= B[j]){

        rk[i] <- rk[i] + 1

      }else

        break
    }
  }

  return(rk)
}

rankAinB_W_for <- function(A, B, W){

  rk <- rep(0, length(A))

  for (i in 1:length(A)){

    for (j in 1:length(B)){

      if (A[i] >= B[j]){

        rk[i] <- rk[i] + W[j]

      }else

        break
    }
  }

  return(rk)
}

