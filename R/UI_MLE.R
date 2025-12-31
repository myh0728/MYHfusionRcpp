GLM.MLE <- function(data = NULL, X.name = NULL, Y.name = NULL,
                    X = NULL, Y = NULL,
                    distribution = "normal", N.binomial = NULL,
                    initial = NULL, wi.boot = NULL, do.SE = TRUE,
                    X.future = NULL)
{
  # 1. 資料前處理
  if (!is.null(data)) {

    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])

  } else {

    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- nrow(Y)
  number_p <- ncol(X)

  # 檢查分佈輸入是否正確
  dist_type <- match.arg(distribution, c("normal", "Bernoulli", "Gamma", "binomial", "Poisson"))

  # -------------------------------------------------------------------------
  # 2. 定義策略 (Strategy Definition)
  # 根據分佈設定初始值、目標函數(obj_fun)、海森函數(hess_fun)與預測邏輯
  # -------------------------------------------------------------------------

  # 預設變數
  param_names <- c("alpha", paste0("beta", 1:number_p))
  start_params <- if(is.null(initial)) c(0, rep(0, number_p)) else initial

  # 用於存放該分佈特有的邏輯
  runner <- list()

  if (dist_type == "normal") {
    # 參數：alpha, beta, sigma (需 log 轉換優化)
    if (is.null(initial)) start_params <- c(start_params, 0) # sigma initial (log scale 0 => sigma=1)
    else start_params[number_p + 2] <- log(start_params[number_p + 2])

    param_names <- c(param_names, "sigma")

    runner$obj_fun <- function(theta) {

      -lL_normal_rcpp(X, Y, theta[1], theta[2:(number_p+1)], exp(theta[number_p+2]), wi.boot)
    }
    runner$hess_fun <- function(theta) {

      -diff_lL_normal_rcpp(X, Y, theta[1], theta[2:(number_p+1)], theta[number_p+2], wi.boot)$hessian
    }
    runner$trans_back <- function(theta) {

      theta[number_p+2] <- exp(theta[number_p+2]) # log-sigma -> sigma
      theta
    }
    runner$predict <- function(theta, X_new) {

      as.vector(theta[1] + X_new %*% theta[2:(number_p+1)])
    }

  } else if (dist_type == "Bernoulli") {
    # 參數：alpha, beta
    runner$obj_fun <- function(theta) {

      -lL_logistic_rcpp(X, Y, theta[1], theta[2:(number_p+1)], wi.boot)
    }
    runner$hess_fun <- function(theta) {

      -diff_lL_logistic_rcpp(X, Y, theta[1], theta[2:(number_p+1)], wi.boot)$hessian
    }
    runner$trans_back <- function(theta) theta # 無需轉換
    runner$predict <- function(theta, X_new) {

      eta <- as.vector(theta[1] + X_new %*% theta[2:(number_p+1)])
      p <- 1 / (1 + exp(-eta))
      list(Y.predict = as.vector((p >= 0.5) * 1), Y.posterior = p)
    }

  } else if (dist_type == "Gamma") {
    # 參數：alpha, beta, nu (需 log 轉換優化)
    if (is.null(initial)) start_params <- c(start_params, 0)
    else start_params[number_p + 2] <- log(start_params[number_p + 2])

    param_names <- c(param_names, "shape") # nu

    runner$obj_fun <- function(theta) {

      -lL_Gamma_rcpp(X, Y, theta[1], theta[2:(number_p+1)], exp(theta[number_p+2]), wi.boot)
    }
    runner$hess_fun <- function(theta) {

      -diff_lL_Gamma_rcpp(X, Y, theta[1], theta[2:(number_p+1)], exp(theta[number_p+2]), wi.boot)$hessian
    }
    runner$trans_back <- function(theta) {

      theta[number_p+2] <- exp(theta[number_p+2])
      theta
    }
    runner$predict <- function(theta, X_new) {

      as.vector(exp(theta[1] + X_new %*% theta[2:(number_p+1)]))
    }

  } else if (dist_type == "binomial") {
    # 參數：alpha, beta
    if (is.null(N.binomial)) N.binomial <- max(Y)

    runner$obj_fun <- function(theta) {

      -lL_binomial_rcpp(X, Y, theta[1], theta[2:(number_p+1)], N.binomial, wi.boot)
    }
    runner$hess_fun <- function(theta) {

      -diff_lL_binomial_rcpp(X, Y, theta[1], theta[2:(number_p+1)], N.binomial, wi.boot)$hessian
    }
    runner$trans_back <- function(theta) theta
    runner$predict <- function(theta, X_new) {

      eta <- as.vector(theta[1] + X_new %*% theta[2:(number_p+1)])
      as.vector(N.binomial / (1 + exp(-eta)))
    }

  } else if (dist_type == "Poisson") {
    # 參數：alpha, beta
    runner$obj_fun <- function(theta) {

      -lL_Poisson_rcpp(X, Y, theta[1], theta[2:(number_p+1)], wi.boot)
    }
    runner$hess_fun <- function(theta) {

      -diff_lL_Poisson_rcpp(X, Y, theta[1], theta[2:(number_p+1)], wi.boot)$hessian
    }
    runner$trans_back <- function(theta) theta
    runner$predict <- function(theta, X_new) {

      as.vector(exp(theta[1] + X_new %*% theta[2:(number_p+1)]))
    }
  }

  # -------------------------------------------------------------------------
  # 3. 執行優化 (Execution)
  # 這段程式碼現在對所有分佈通用
  # -------------------------------------------------------------------------

  # A. MLE Optimization
  opt_res <- nlminb(start = start_params, objective = runner$obj_fun)

  # B. 參數處理 (Back-transformation & Naming)
  theta.hat.trans <- opt_res$par
  theta.hat <- runner$trans_back(theta.hat.trans)
  names(theta.hat) <- param_names

  # 儲存基本結果
  results <- list(
    alpha = theta.hat[1],
    beta  = theta.hat[2:(number_p + 1)]
  )

  # 針對有額外參數的分佈 (Normal, Gamma) 加入結果
  if (dist_type == "normal") results$sigma <- theta.hat["sigma"]
  if (dist_type == "Gamma")  results$nu    <- theta.hat["shape"]

  results$parameter <- theta.hat

  # -------------------------------------------------------------------------
  # 4. 計算標準誤 (SE Calculation)
  # 這段程式碼現在對所有分佈通用
  # -------------------------------------------------------------------------
  if (do.SE) {

    hessian_mat <- runner$hess_fun(theta.hat) # 傳入已還原的參數 (sigma, nu)

    # 計算 Covariance
    asy.Cov <- inv_sympd_rcpp(hessian_mat)

    # 設定名稱
    dimnames(asy.Cov) <- list(param_names, param_names)

    results$Cov.coef <- asy.Cov / number_n
  }

  # -------------------------------------------------------------------------
  # 5. 預測 (Prediction)
  # -------------------------------------------------------------------------
  if (!is.null(X.future)) {

    X.future <- as.matrix(X.future)
    pred_res <- runner$predict(theta.hat, X.future)

    if (is.list(pred_res)) {

      results <- c(results, pred_res) # 用於 Bernoulli 同時回傳類別和機率

    } else {

      results$Y.predict <- pred_res
    }
  }

  return(results)
}




