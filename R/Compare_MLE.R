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

diff_lL_normal <- function(X, Y, alpha, beta, sigma, wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  # -------------------------------------------------------------------------
  # 1. 處理權重 (Weight Handling)
  # -------------------------------------------------------------------------
  if (is.null(wi.boot)) {
    # 無權重：權重全為 1
    w <- rep(1, number_n)
    sum_w <- number_n
  } else {
    # 有權重：使用傳入的權重
    if (length(wi.boot) != number_n) {
      stop("Length of wi.boot must match number of rows in X")
    }
    w <- as.vector(wi.boot)
    sum_w <- sum(w)
  }

  gradient <- rep(0, number_p + 2)
  hessian <- matrix(0, nrow = number_p + 2, ncol = number_p + 2)

  # -------------------------------------------------------------------------
  # 2. 計算基礎統計量 (Basic Statistics)
  # -------------------------------------------------------------------------
  # 擴充 X 矩陣 (加入截距項的 Column)
  X_aug <- cbind(1, X)

  # 殘差
  res <- as.vector(Y - alpha - X %*% beta)

  # 加權 RSS (Weighted Residual Sum of Squares)
  # formula: sum(w * res^2)
  rss <- sum(w * res^2)

  # 加權 X * res
  # formula: sum(w * X_aug * res) -> X_aug^T * (w * res)
  # colSums(w * X_aug * res) 也可以，但在 R 中 crossprod 效率較好
  Xrss <- as.vector(crossprod(X_aug, w * res))

  # -------------------------------------------------------------------------
  # 3. 計算 Gradient
  # -------------------------------------------------------------------------
  # 對 alpha 和 beta 的微分
  gradient[1:(number_p + 1)] <- Xrss / (sigma ^ 2) / number_n

  # 對 sigma 的微分
  # dL/d_sigma = -sum(w)/sigma + sum(w * res^2)/sigma^3
  # 最後除以 number_n (平均梯度)
  gradient[number_p + 2] <- rss / (sigma ^ 3) / number_n - sum_w / sigma / number_n

  # -------------------------------------------------------------------------
  # 4. 計算 Hessian
  # -------------------------------------------------------------------------
  # Alpha & Beta Block (左上角)
  # d^2L / d_beta^2 = - X^T * W * X / sigma^2
  # 這裡使用 crossprod(X, w*X) 來計算加權的 X^T X
  # 並除以 number_n 以符合期望值尺度 (取代原本的 eXsq_rcpp)
  XtWX <- crossprod(X_aug, w * X_aug)
  hessian[1:(number_p + 1), 1:(number_p + 1)] <- -XtWX / (sigma ^ 2) / number_n

  # Sigma Block (右下角)
  # d^2L / d_sigma^2 = sum(w)/sigma^2 - 3 * sum(w * res^2) / sigma^4
  hessian[number_p + 2, number_p + 2] <- -rss * 3 / (sigma ^ 4) / number_n + sum_w / (sigma ^ 2) / number_n

  # Interaction Block (Alpha/Beta 對 Sigma 的微分)
  # d^2L / d_beta d_sigma = -2 * X^T * W * res / sigma^3
  hessian[1:(number_p + 1), number_p + 2] <- -Xrss * 2 / (sigma ^ 3) / number_n

  # 對稱填補
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

diff_lL_logistic <- function(X, Y, alpha, beta, wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]

  # -------------------------------------------------------------------------
  # 1. 處理權重 (Weight Handling)
  # -------------------------------------------------------------------------
  if (is.null(wi.boot)) {
    # 無權重：全為 1
    w <- rep(1, number_n)
  } else {
    # 有權重
    if (length(wi.boot) != number_n) {
      stop("Length of wi.boot must match number of rows in X")
    }
    w <- as.vector(wi.boot)
  }

  # -------------------------------------------------------------------------
  # 2. 計算機率與中間項
  # -------------------------------------------------------------------------
  # P(Y=1) = 1 / (1 + exp(-(alpha + X beta)))
  # 使用原本的邏輯計算
  eta <- as.vector(alpha + X %*% beta)
  enSI <- exp(-eta)
  numerator <- 1 + enSI

  # p = 1 / (1 + exp(-eta))
  p <- 1 / numerator

  # 準備擴充矩陣 [1, X]
  X_intercept <- cbind(1, X)

  # -------------------------------------------------------------------------
  # 3. 計算 Gradient
  # -------------------------------------------------------------------------
  # Gradient = X^T * W * (Y - p) / n
  res <- Y - p

  # crossprod(A, b) 等同於 t(A) %*% b，但在 R 中效率更高
  gradient <- as.vector(crossprod(X_intercept, w * res)) / number_n

  # -------------------------------------------------------------------------
  # 4. 計算 Hessian
  # -------------------------------------------------------------------------
  # Hessian = - X^T * W_hess * X / n
  # 其中 Logistic 的權重 W_hess = w * p * (1-p)
  # p * (1-p) = (1/num) * (enSI/num) = enSI / (num^2)

  p_times_1_minus_p <- enSI / (numerator^2)
  hess_weights <- w * p_times_1_minus_p

  # 計算加權的 X^T X
  # crossprod(X, weights * X)
  hessian <- -crossprod(X_intercept, hess_weights * X_intercept) / number_n

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

diff_lL_Gamma <- function(X, Y, alpha, beta, nu, wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  # -------------------------------------------------------------------------
  # 1. 處理權重 (Weight Handling)
  # -------------------------------------------------------------------------
  if (is.null(wi.boot)) {
    # 無權重：全為 1
    w <- rep(1, number_n)
    sum_w <- number_n
  } else {
    # 有權重
    if (length(wi.boot) != number_n) {
      stop("Length of wi.boot must match number of rows in X")
    }
    w <- as.vector(wi.boot)
    sum_w <- sum(w)
  }

  gradient <- rep(0, number_p + 2)
  hessian <- matrix(0, nrow = number_p + 2, ncol = number_p + 2)

  # -------------------------------------------------------------------------
  # 2. 計算中間變數
  # -------------------------------------------------------------------------
  # Linear predictor: eta = alpha + X %*% beta
  eta <- as.vector(alpha + X %*% beta)

  # eSI_i = exp(eta)
  eSI_i <- exp(eta)

  # lambda = nu / exp(eta)
  lambda_i <- nu / eSI_i

  # 擴充矩陣 [1, X] 用於計算截距項與 beta
  X_intercept <- cbind(1, X)

  # -------------------------------------------------------------------------
  # 3. 計算 Gradient
  # -------------------------------------------------------------------------
  # dL/d_beta (包含 alpha): sum( w * X * (Y * lambda - nu) ) / n
  res_mean <- Y * lambda_i - nu
  gradient[1:(number_p + 1)] <- as.vector(crossprod(X_intercept, w * res_mean)) / number_n

  # dL/d_nu: sum( w * (log(lambda * Y) - Y/eSI + 1 - digamma(nu)) ) / n
  term_nu <- log(lambda_i * Y) - Y / eSI_i
  gradient[number_p + 2] <- (sum(w * term_nu) + sum_w * (1 - digamma(nu))) / number_n

  # -------------------------------------------------------------------------
  # 4. 計算 Hessian
  # -------------------------------------------------------------------------
  # Alpha & Beta Block (左上角)
  # d^2L / d_beta^2 = - sum( w * Y * lambda * X * X^T ) / n
  # Hessian 權重: w * Y * lambda
  hess_weights <- w * Y * lambda_i

  # 計算加權 X^T X
  hessian[1:(number_p + 1), 1:(number_p + 1)] <- -crossprod(X_intercept, hess_weights * X_intercept) / number_n

  # Nu Block (右下角)
  # d^2L / d_nu^2 = sum( w * (1/nu - trigamma(nu)) ) / n
  hessian[number_p + 2, number_p + 2] <- (sum_w / number_n) * (1 / nu - trigamma(nu))

  # 根據原函數邏輯，Alpha/Beta 與 Nu 之間的交互項 (Cross terms) 設為 0

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

diff_lL_binomial <- function(X, Y, alpha, beta, N, wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]

  # -------------------------------------------------------------------------
  # 1. 處理權重 (Weight Handling)
  # -------------------------------------------------------------------------
  if (is.null(wi.boot)) {
    # 無權重：全為 1
    w <- rep(1, number_n)
  } else {
    # 有權重
    if (length(wi.boot) != number_n) {
      stop("Length of wi.boot must match number of rows in X")
    }
    w <- as.vector(wi.boot)
  }

  # -------------------------------------------------------------------------
  # 2. 計算機率與中間變數
  # -------------------------------------------------------------------------
  # Linear Predictor: SI (Single Index) = alpha + X beta
  SI <- as.vector(alpha + X %*% beta)

  # Probability: P1 = 1 / (1 + exp(-SI))
  # 為了數值穩定，也可以寫成 plogis(SI)
  P1 <- 1 / (1 + exp(-SI))

  # 準備擴充矩陣 [1, X]
  X_intercept <- cbind(1, X)

  # -------------------------------------------------------------------------
  # 3. 計算 Gradient
  # -------------------------------------------------------------------------
  # Score function for Binomial: Y - N * p
  # Gradient = X^T * W * (Y - N*p) / n
  res <- Y - N * P1

  # 使用 crossprod 加速矩陣運算
  gradient <- as.vector(crossprod(X_intercept, w * res)) / number_n

  # -------------------------------------------------------------------------
  # 4. 計算 Hessian
  # -------------------------------------------------------------------------
  # Hessian = - X^T * W_hess * X / n
  # Variance of Binomial: Var(Y) = N * p * (1-p)
  # Hessian weights = w * N * p * (1-p)

  hess_weights <- w * N * P1 * (1 - P1)

  # 計算加權的 X^T X
  hessian <- -crossprod(X_intercept, hess_weights * X_intercept) / number_n

  results <- list(gradient = gradient,
                  hessian = hessian)

  return(results)
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

diff_lL_Poisson <- function(X, Y, alpha, beta, wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]

  # -------------------------------------------------------------------------
  # 1. 處理權重 (Weight Handling)
  # -------------------------------------------------------------------------
  if (is.null(wi.boot)) {
    # 無權重：全為 1
    w <- rep(1, number_n)
  } else {
    # 有權重
    if (length(wi.boot) != number_n) {
      stop("Length of wi.boot must match number of rows in X")
    }
    w <- as.vector(wi.boot)
  }

  # -------------------------------------------------------------------------
  # 2. 計算 Lambda 與中間變數
  # -------------------------------------------------------------------------
  # Linear Predictor: SI = alpha + X beta
  SI <- as.vector(alpha + X %*% beta)

  # Poisson Mean: lambda = exp(SI)
  lambda <- exp(SI)

  # 準備擴充矩陣 [1, X]
  X_intercept <- cbind(1, X)

  # -------------------------------------------------------------------------
  # 3. 計算 Gradient
  # -------------------------------------------------------------------------
  # Gradient = X^T * W * (Y - lambda) / n
  res <- Y - lambda

  # 使用 crossprod 加速矩陣運算
  gradient <- as.vector(crossprod(X_intercept, w * res)) / number_n

  # -------------------------------------------------------------------------
  # 4. 計算 Hessian
  # -------------------------------------------------------------------------
  # Hessian = - X^T * W_hess * X / n
  # Hessian weights = w * lambda

  hess_weights <- w * lambda

  # 計算加權的 X^T X
  # crossprod(X, weights * X)
  hessian <- -crossprod(X_intercept, hess_weights * X_intercept) / number_n

  results <- list(gradient = gradient,
                  hessian = hessian)

  return(results)
}
