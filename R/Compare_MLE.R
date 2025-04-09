lL_normal <- function(X, Y, alpha, beta, sigma, wi.boot = NULL)
{
  SI <- alpha + X %*% beta
  f.cond.y.x <- dnorm(Y - SI, mean = 0, sd = sigma)

  if (is.null(wi.boot))
  {
    lL <- sum(log(f.cond.y.x))
  }else
  {
    lL <- sum(log(f.cond.y.x) * wi.boot)
  }

  return(lL)
}

diff_lL_normal <- function(X, Y, alpha, beta, sigma)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  gradient <- rep(0, number_p + 2)
  hessian <- matrix(0, nrow = number_p + 2, ncol = number_p + 2)

  res <- as.vector(Y - alpha - X %*% beta)
  rss <- sum(res ^ 2)
  Xrss <- colSums(cbind(1, X) * res)
  gradient[1:(number_p + 1)] <- Xrss / (sigma ^ 2) / number_n
  gradient[number_p + 2] <- rss / (sigma ^ 3) / number_n - 1 / sigma
  hessian[1:(number_p + 1), 1:(number_p + 1)] <- -eXsq_rcpp(cbind(1, X)) / (sigma ^ 2)
  hessian[number_p + 2, number_p + 2] <- -rss * 3 / (sigma ^ 4) / number_n + 1 / (sigma ^ 2)
  hessian[1:(number_p + 1), number_p + 2] <- -Xrss * 2 / (sigma ^ 3) / number_n
  hessian[number_p + 2, 1:(number_p + 1)] <- t(hessian[1:(number_p + 1), number_p + 2])

  results <- list(gradient = gradient,
                  hessian = hessian)

  return(results)
}

lL_logistic <- function(X, Y, alpha, beta, wi.boot = NULL)
{
  SI <- alpha + X %*% beta
  P1 <- 1 / (1 + exp(-SI))
  f.cond.y.x <- P1
  f.cond.y.x[Y == 0] <- (1 - P1[Y == 0])

  if (is.null(wi.boot))
  {
    lL <- sum(log(f.cond.y.x))
  }else
  {
    lL <- sum(log(f.cond.y.x) * wi.boot)
  }

  return(lL)
}

diff_lL_logistic <- function(X, Y, alpha, beta)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]

  enSI <- as.vector(exp(-(alpha + X %*% beta)))
  numerator <- 1 + enSI
  X_intercept <- cbind(1, X)
  gradient <- colSums(X_intercept * (Y - 1 / numerator)) / number_n
  hessian <- -eXsq_w_rcpp(X_intercept / numerator, enSI)

  results <- list(gradient = gradient,
                  hessian = hessian)

  return(results)
}

lL_Gamma <- function(X, Y, alpha, beta, nu, wi.boot = NULL)
{
  SI <- alpha + X %*% beta
  lambda <- nu * exp(-SI)
  f.cond.y.x <- (lambda ^ nu) * (Y ^ (nu - 1)) * exp(-lambda * Y) / gamma(nu)

  if (is.null(wi.boot))
  {
    lL <- sum(log(f.cond.y.x))
  }else
  {
    lL <- sum(log(f.cond.y.x) * wi.boot)
  }

  return(lL)
}

diff_lL_Gamma <- function(X, Y, alpha, beta, nu)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  gradient <- rep(0, number_p + 2)
  hessian <- matrix(0, nrow = number_p + 2, ncol = number_p + 2)

  eSI_i <- as.vector(exp(alpha + X %*% beta))
  lambda_i <- nu / eSI_i

  gradient[1:(number_p + 1)] <- colMeans(cbind(1, X) * (Y * lambda_i - nu))
  gradient[number_p + 2] <- mean(log(lambda_i * Y) - Y / eSI_i) + 1 - digamma(nu)
  hessian[1:(number_p + 1), 1:(number_p + 1)] <- -eXsq_w_R(cbind(1, X), Y * lambda_i)
  hessian[number_p + 2, number_p + 2] <- 1 / nu - trigamma(nu)

  results <- list(gradient = gradient,
                  hessian = hessian)

  return(results)
}

lL_binomial <- function(X, Y, alpha, beta, N, wi.boot = NULL)
{
  SI <- alpha + X %*% beta
  P1 <- 1 / (1 + exp(-SI))

  if (is.null(wi.boot))
  {
    lL <- sum(Y * log(P1) + (N - Y) * log(1 - P1))
  }else
  {
    lL <- sum((Y * log(P1) + (N - Y) * log(1 - P1)) * wi.boot)
  }

  return(lL)
}

lL_Poisson <- function(X, Y, alpha, beta, wi.boot = NULL)
{
  SI <- alpha + X %*% beta
  eSI <- exp(SI)

  if (is.null(wi.boot))
  {
    lL <- sum(Y * SI - eSI)
  }else
  {
    lL <- sum((Y * SI - eSI) * wi.boot)
  }

  return(lL)
}
