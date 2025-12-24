SemiSI <- function(data = NULL, X.name = NULL, Y.name = NULL,
                   X = NULL, Y = NULL,
                   type = "mean", kernel = "K2.Biweight",
                   bandwidth = NULL, bandwidth.initial = 1,
                   beta.initial = NULL,
                   wi.boot = NULL)
{
  # ============================================================================
  # 1. 資料前處理 (Data Preprocessing)
  # ============================================================================
  if (!is.null(data)) {

    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])

  } else {

    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (number_p < 2) {
    warning("There should be at least two covariates.")
  }

  # 處理初始 Beta (標準化：設第一個分量為 1)
  if (is.null(beta.initial)) {

    beta.initial <- rep(0, number_p - 1)

  } else {

    beta.initial <- as.vector(beta.initial)
    # 若使用者傳入完整 beta，將其標準化為 beta[-1]/beta[1]
    # 假設 beta[1] 不為 0
    beta.initial <- beta.initial[-1] / beta.initial[1]
  }

  # ============================================================================
  # 2. 設定參數與選擇函數 (Configuration)
  # ============================================================================

  # 2.1 Kernel 名稱對應 (R -> C++)
  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  # 2.2 根據 type 決定要使用的 C++ 函數 (CV 與 LS)
  if (type == "mean") {

    cv_func <- CVKNW_rcpp
    ls_func <- LSKNW_rcpp

  } else if (type == "distribution") {

    cv_func <- CVKNWcdf_rcpp
    ls_func <- LSKNWcdf_rcpp

  } else {

    stop("Invalid 'type'. Must be 'mean' or 'distribution'.")
  }

  results <- list()

  # ============================================================================
  # 3. 步驟一：頻寬選擇 (Bandwidth Selection via LOOCV) - Optional
  #    如果沒有提供 bandwidth，則同時優化 beta 和 bandwidth
  # ============================================================================
  if (is.null(bandwidth)) {

    # 定義 CV 目標函數 (theta 包含 beta 和 log(h))
    LOOCV_Obj <- function(theta) {
      # 重組 beta: 固定第一個係數為 1
      beta_curr <- c(1, theta[1:(number_p - 1)])
      # 取出 bandwidth (指數還原，確保為正)
      h_curr <- rep(exp(theta[number_p]), 1) # Single Index 是 1 維，h 長度為 1

      # 計算 Single Index (X * beta)
      # 注意：C++ 函數接受矩陣 X，這裡傳入 (n x 1) 的矩陣
      index_val <- X %*% beta_curr

      # 呼叫對應的 C++ CV 函數
      # wi.boot 若為 NULL，C++ 端會自動處理
      cv_val <- cv_func(Y = Y, X = index_val, h = h_curr, kernel = cpp_kernel, w = wi.boot)

      return(cv_val)
    }

    # 執行優化
    # 參數：前 (p-1) 個是 beta，最後 1 個是 log(bandwidth)
    start_params <- c(beta.initial, log(bandwidth.initial))
    LOOCV.esti <- nlminb(start = start_params, objective = LOOCV_Obj)

    # 提取結果
    beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
    bandwidth <- exp(LOOCV.esti$par[number_p])

    results <- list(coef.CV = c(1, beta.CV),
                    bandwidth = bandwidth,
                    CV.value = LOOCV.esti$objective)

    # 更新 beta.initial 供下一步使用
    beta.initial <- beta.CV

  } else {
    # 若 bandwidth 已固定，確保格式正確
    bandwidth <- rep(bandwidth, length.out = 1)
  }

  # ============================================================================
  # 4. 步驟二：係數估計 (Coefficient Estimation via Least Squares)
  #    使用確定的 bandwidth，優化 beta
  # ============================================================================

  # 定義 LS 目標函數 (theta 僅包含 beta)
  SemiLS_Obj <- function(theta) {

    beta_curr <- c(1, theta)
    index_val <- X %*% beta_curr

    # 呼叫對應的 C++ LS 函數
    ls_val <- ls_func(Y = Y, X = index_val, h = bandwidth, kernel = cpp_kernel, w = wi.boot)

    return(ls_val)
  }

  # 執行優化
  SemiLS.esti <- nlminb(start = beta.initial, objective = SemiLS_Obj)

  # 儲存最終結果
  results$coef <- c(1, SemiLS.esti$par)
  results$LS.value <- SemiLS.esti$objective

  return(results)
}

SemiMI <- function(data = NULL, X.name = NULL, Y.name = NULL,
                   X = NULL, Y = NULL,
                   type = "mean", dimension = 2, kernel = "K2.Biweight",
                   bandwidth = NULL, bandwidth.initial = 1,
                   beta.initial = NULL,
                   wi.boot = NULL)
{
  # ============================================================================
  # 1. 資料前處理 (Data Preprocessing)
  # ============================================================================
  if (!is.null(data)) {

    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])

  } else {

    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  # 基本檢查
  if (number_p <= dimension) {

    warning("The number of covariates should be larger than the number of indices (dimension).")
  }

  # 提早檢查 type
  if (!type %in% c("mean", "distribution")) {

    stop("Invalid 'type'. Must be 'mean' or 'distribution'.")
  }

  # 計算自由參數的數量 (排除掉固定為 Identity 的部分)
  # beta 是 p x d 矩陣，前 d x d 固定為 I，剩下 (p-d) x d 為自由參數
  number_k <- (number_p - dimension) * dimension

  # ============================================================================
  # 2. 初始參數設定 (Initialization)
  # ============================================================================

  # 處理 beta.initial，將其轉換為優化器需要的向量格式 (長度 number_k)
  if (is.null(beta.initial)) {
    # 預設為 0
    theta_start <- rep(0, number_k)

  } else {
    # 確保輸入是矩陣 p x d
    beta_mat <- matrix(beta.initial, nrow = number_p, ncol = dimension)
    # 取出下半部 (p-d) x d 的部分作為起始值
    lower_block <- beta_mat[(dimension + 1):number_p, , drop = FALSE]
    theta_start <- as.vector(lower_block)
  }

  # ============================================================================
  # 3. 函數與 Kernel 選擇 (Function Selection)
  # ============================================================================

  # Kernel 名稱對應
  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  # 選擇 C++ 函數
  if (type == "mean") {

    cv_func <- CVKNW_rcpp
    ls_func <- LSKNW_rcpp

  } else {

    cv_func <- CVKNWcdf_rcpp
    ls_func <- LSKNWcdf_rcpp
  }

  # ============================================================================
  # 4. 輔助函數：重建 Beta 矩陣 (Helper: Reconstruct Beta)
  #    將優化向量 theta 轉換回完整的 p x d 矩陣
  # ============================================================================
  reconstruct_beta <- function(theta_vec) {

    lower_mat <- matrix(theta_vec, nrow = number_p - dimension, ncol = dimension)
    # 上半部固定為單位矩陣 (Identity)，下半部為優化參數
    rbind(diag(dimension), lower_mat)
  }

  results <- list()

  # ============================================================================
  # 5. 步驟一：頻寬選擇 (Bandwidth Selection via LOOCV)
  # ============================================================================
  if (is.null(bandwidth)) {

    # 定義 CV 目標函數 (theta 包含 beta 參數 和 log(h))
    LOOCV_Obj <- function(theta) {
      # 1. 拆解參數
      theta_beta <- theta[1:number_k]
      h_scalar <- exp(theta[number_k + 1]) # 頻寬縮放因子 (scalar)

      # 2. 重組 beta 並投影 X
      beta_curr <- reconstruct_beta(theta_beta)
      X_beta <- X %*% beta_curr

      # 3. 計算實際頻寬向量 (Scalar * SD of each index)
      # apply(..., 2, sd) 會回傳長度為 dimension 的向量
      h_curr <- h_scalar * apply(X_beta, 2, sd)

      # 4. 呼叫 C++
      cv_val <- cv_func(Y = Y, X = X_beta, h = h_curr, kernel = cpp_kernel, w = wi.boot)

      return(cv_val)
    }

    # 執行優化
    # 參數：beta 參數 + 1 個 log(bandwidth)
    start_params <- c(theta_start, log(bandwidth.initial))
    LOOCV.esti <- nlminb(start = start_params, objective = LOOCV_Obj)

    # 提取結果
    theta_beta_opt <- LOOCV.esti$par[1:number_k]
    h_scalar_opt <- exp(LOOCV.esti$par[number_k + 1])

    beta.CV <- reconstruct_beta(theta_beta_opt)

    # 計算最終頻寬向量
    bandwidth <- h_scalar_opt * apply(X %*% beta.CV, 2, sd)
    results <- list(coef.CV = beta.CV,
                    bandwidth = bandwidth,
                    CV.value = LOOCV.esti$objective)

    # 更新 theta_start 供下一步使用
    theta_start <- theta_beta_opt

  } else {
    # 若頻寬已固定，確保長度正確 (應為 dimension 長度)
    bandwidth <- rep(bandwidth, length.out = dimension)
  }

  # ============================================================================
  # 6. 步驟二：係數估計 (Coefficient Estimation via Least Squares)
  # ============================================================================

  # 定義 LS 目標函數 (theta 僅包含 beta 參數)
  SemiLS_Obj <- function(theta) {

    beta_curr <- reconstruct_beta(theta)
    X_beta <- X %*% beta_curr

    # 使用確定的 bandwidth (向量)
    ls_val <- ls_func(Y = Y, X = X_beta, h = bandwidth, kernel = cpp_kernel, w = wi.boot)

    return(ls_val)
  }

  # 執行優化
  SemiLS.esti <- nlminb(start = theta_start, objective = SemiLS_Obj)

  # 儲存最終結果
  results$coef <- reconstruct_beta(SemiLS.esti$par)
  results$LS.value <- SemiLS.esti$objective

  return(results)
}

SemiSI.adtLASSO.ADMM <- function(data = NULL, X.name = NULL, Y.name = NULL,
                                 X = NULL, Y = NULL,
                                 type = "mean", kernel = "K2.Biweight",
                                 beta.initial, bandwidth,
                                 weight.beta, lambda.L1.penalty,
                                 rho.augLagrange.penalty = 1,
                                 iter.max = 50, tol = 1e-7)
{
  # ============================================================================
  # 1. 資料前處理 (Data Preprocessing)
  # ============================================================================
  if (!is.null(data)) {

    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])

  } else {

    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (number_p < 2) {

    stop("There should be at least two covariates.")
  }

  if (!type %in% c("mean", "distribution")) {

    stop("Invalid type. Must be 'mean' or 'distribution'.")
  }

  # ============================================================================
  # 2. 參數初始化與設定 (Initialization & Configuration)
  # ============================================================================

  # 2.1 轉換 Kernel 名稱 (R -> C++)
  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  # 2.2 選擇目標函數 (根據 type)
  # 使用 Unified C++ 函數
  if (type == "mean") {

    ls_func <- LSKNW_rcpp

  } else {

    ls_func <- LSKNWcdf_rcpp
  }

  # 2.3 初始化 ADMM 變數
  # 假設 beta[1] 固定為 1，我們只優化剩下的 p-1 個參數
  beta.free.initial <- beta.initial[-1] / beta.initial[1]

  # 計算 Adaptive Weights (對應 free parameters)
  weight.free.beta <- lambda.L1.penalty * weight.beta[-1]

  # ADMM 變數
  beta_free_iter  <- beta.free.initial
  gamma_free_iter <- beta.free.initial
  u_free_iter     <- rep(0, length(beta.free.initial))

  # 確保 bandwidth 格式正確
  h_curr <- rep(bandwidth, length.out = 1)

  # ============================================================================
  # 3. ADMM 主迴圈 (Main Loop)
  # ============================================================================
  converged <- FALSE
  final_iter <- 0

  for (iter in 1:iter.max) {

    final_iter <- iter

    # --------------------------------------------------------------------------
    # Step 1: Beta Update (Primal Variable Update)
    # 目標：最小化 (Loss + Augmented Lagrangian Penalty)
    # --------------------------------------------------------------------------
    SemiLS.penalized <- function(theta) {
      # 重組 beta (beta[1] = 1)
      beta_full <- c(1, theta)

      # 計算 Single Index
      index_val <- X %*% beta_full

      # 計算 Kernel Smoothing Loss (Scaled by n)
      # 注意：C++ 函數回傳 Mean Squared Error，這裡乘回 n 變成 Sum of Squares
      # 以便與 Penalty 項的量級匹配
      loss_val <- ls_func(Y = Y, X = index_val, h = h_curr, kernel = cpp_kernel) * number_n

      # 計算 Augmented Lagrangian Penalty
      # (rho / 2) * || beta - gamma + u ||^2
      penalty_val <- (rho.augLagrange.penalty / 2) * sum((theta - gamma_free_iter + u_free_iter) ^ 2)

      return(loss_val + penalty_val)
    }

    # 執行數值優化 (nlminb)
    # start 使用上一次迭代的 beta
    SemiLS.penalized.esti <- nlminb(start = beta_free_iter, objective = SemiLS.penalized)
    beta_free_iter_new <- SemiLS.penalized.esti$par

    # --------------------------------------------------------------------------
    # Step 2: Gamma Update (Dual Variable Update)
    # 目標：L1 Proximal Operator (Soft Thresholding)
    # --------------------------------------------------------------------------
    # 變數變換： v = beta + u
    v <- beta_free_iter_new + u_free_iter

    # 門檻值： lambda_j / rho
    threshold <- weight.free.beta / rho.augLagrange.penalty

    # Soft Thresholding: sign(v) * max(|v| - threshold, 0)
    gamma_free_iter_new <- sign(v) * pmax(abs(v) - threshold, 0)

    # --------------------------------------------------------------------------
    # Step 3: U Update (Multiplier Update)
    # --------------------------------------------------------------------------
    u_free_iter_new <- u_free_iter + beta_free_iter_new - gamma_free_iter_new

    # --------------------------------------------------------------------------
    # Step 4: 收斂檢查 (Convergence Check)
    # --------------------------------------------------------------------------
    # Primal Residual: || beta - gamma ||
    res_primal <- sqrt(sum((beta_free_iter_new - gamma_free_iter_new) ^ 2))

    # Dual Residual: || rho * (gamma_new - gamma_old) ||
    res_dual <- rho.augLagrange.penalty * sqrt(sum((gamma_free_iter_new - gamma_free_iter) ^ 2))

    # 更新迭代變數
    beta_free_iter  <- beta_free_iter_new
    gamma_free_iter <- gamma_free_iter_new
    u_free_iter     <- u_free_iter_new

    # 檢查是否收斂
    if (res_primal < tol && res_dual < tol) {

      converged <- TRUE
      break
    }

    # Optional: 顯示進度
    # if (iter %% 10 == 0) cat(sprintf("Iter %d: Primal Res = %.2e, Dual Res = %.2e\n", iter, res_primal, res_dual))
  }

  # ============================================================================
  # 4. 結果整理 (Results)
  # ============================================================================

  if (!converged) {

    warning("ADMM did not converge within the maximum number of iterations.")
  }

  # 最終係數：使用 gamma (因為它是稀疏的)
  final_coef <- c(1, gamma_free_iter)

  # 原始優化變數 (非稀疏)
  beta_coef <- c(1, beta_free_iter)

  results <- list(coef = final_coef,           # 稀疏係數 (Lasso 結果)
                  coef.beta = beta_coef,       # 原始係數 (Ridge-like 結果)
                  iteration = final_iter,
                  converged = converged,
                  residuals = c(primal = res_primal, dual = res_dual))

  return(results)
}

SemiMI.adtLASSO.ADMM <- function(data = NULL, X.name = NULL, Y.name = NULL,
                                 X = NULL, Y = NULL,
                                 type = "mean", dimension = 2, kernel = "K2.Biweight",
                                 beta.initial, bandwidth,
                                 weight.beta, lambda.L1.penalty,
                                 rho.augLagrange.penalty = 1,
                                 iter.max = 50, tol = 1e-7)
{
  # ============================================================================
  # 1. 資料前處理 (Data Preprocessing)
  # ============================================================================
  if (!is.null(data)) {

    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])

  } else {

    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (number_p <= dimension) {

    stop("The number of covariates must be larger than the dimension (indices) for variable selection.")
  }

  if (!type %in% c("mean", "distribution")) {

    stop("Invalid type. Must be 'mean' or 'distribution'.")
  }

  # ============================================================================
  # 2. 參數初始化與設定 (Configuration)
  # ============================================================================

  # 2.1 Kernel 選擇
  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  # 2.2 目標函數選擇
  if (type == "mean") {

    ls_func <- LSKNW_rcpp

  } else {

    ls_func <- LSKNWcdf_rcpp
  }

  # 2.3 提取自由參數 (Free Parameters)
  # beta 是 p x d，前 d x d 為 Identity (固定)，只優化後 (p-d) x d
  # beta.initial 應為 p x d 矩陣
  if (is.vector(beta.initial)) {

    beta.initial <- matrix(beta.initial, nrow = number_p, ncol = dimension)
  }

  # 自由參數起始值 (矩陣格式: (p-d) x d)
  beta_free_mat <- beta.initial[(dimension + 1):number_p, , drop = FALSE]

  # 轉換為向量格式供 nlminb 使用
  theta_start <- as.vector(beta_free_mat)

  # 2.4 ADMM 變數初始化 (皆為向量格式，長度 = (p-d)*d)
  beta_vec_iter  <- theta_start
  gamma_vec_iter <- theta_start
  u_vec_iter     <- rep(0, length(theta_start))

  # 2.5 處理權重 (Adaptive Weights)
  # weight.beta 預期長度為 p (對應每個變數)
  # 我們只取出需要被懲罰的變數權重 (第 d+1 到 p 個變數)
  if (length(weight.beta) != number_p) {

    stop("weight.beta length must equal number of covariates (p).")
  }
  weight_free <- weight.beta[(dimension + 1):number_p] # 長度為 p-d

  # 確保 bandwidth 格式
  h_curr <- rep(bandwidth, length.out = dimension)

  # 輔助函數：將向量還原為 (p-d) x d 矩陣
  vec_to_mat <- function(v) {

    matrix(v, nrow = number_p - dimension, ncol = dimension)
  }

  # 輔助函數：重建完整 beta (p x d)
  reconstruct_full_beta <- function(free_mat) {

    rbind(diag(dimension), free_mat)
  }

  # ============================================================================
  # 3. ADMM 主迴圈 (Group LASSO)
  # ============================================================================
  converged <- FALSE
  final_iter <- 0

  for (iter in 1:iter.max) {

    final_iter <- iter

    # --------------------------------------------------------------------------
    # Step 1: Beta Update (Smooth Optimization)
    # 目標：最小化 Kernel Loss + Augmented Lagrangian Quadratic Penalty
    # --------------------------------------------------------------------------
    SemiLS.penalized <- function(theta) {
      # 1. 重建完整 Beta 用於計算 Kernel Loss
      beta_free_curr <- vec_to_mat(theta)
      beta_full <- reconstruct_full_beta(beta_free_curr)

      # 投影數據
      index_val <- X %*% beta_full

      # 計算 Loss (MSE * n = Sum of Squares)
      loss_val <- ls_func(Y = Y, X = index_val, h = h_curr, kernel = cpp_kernel) * number_n

      # 計算 Penalty (Ridge-like part of ADMM)
      # (rho / 2) * || beta - gamma + u ||^2
      diff_vec <- theta - gamma_vec_iter + u_vec_iter
      penalty_val <- (rho.augLagrange.penalty / 2) * sum(diff_vec ^ 2)

      return(loss_val + penalty_val)
    }

    # 執行優化
    # start 使用上一次迭代的 beta
    SemiLS.esti <- nlminb(start = beta_vec_iter, objective = SemiLS.penalized)
    beta_vec_iter_new <- SemiLS.esti$par

    # --------------------------------------------------------------------------
    # Step 2: Gamma Update (Group LASSO Proximal Operator)
    # --------------------------------------------------------------------------
    # 準備變數 v = beta + u
    v_vec <- beta_vec_iter_new + u_vec_iter

    # 將 V 還原為矩陣形式 (Rows = Variables, Cols = Indices)
    # 矩陣大小：(p-d) x d
    V_mat <- vec_to_mat(v_vec)
    Gamma_mat <- matrix(0, nrow = nrow(V_mat), ncol = ncol(V_mat))

    # 針對每一個變數 (Row) 進行 Group Soft Thresholding
    for (j in 1:nrow(V_mat)) {
      # 取出第 j 個變數對應的所有 Index 係數向量
      v_j <- V_mat[j, ]

      # 計算該向量的 L2 Norm
      norm_v_j <- sqrt(sum(v_j^2))

      # 計算閾值 (lambda * w_j / rho)
      threshold <- (lambda.L1.penalty * weight_free[j]) / rho.augLagrange.penalty

      # Group Soft Thresholding 公式:
      # S(v, t) = max(0, 1 - t / ||v||) * v
      if (norm_v_j > threshold) {

        scaling_factor <- 1 - (threshold / norm_v_j)
        Gamma_mat[j, ] <- scaling_factor * v_j

      } else {

        Gamma_mat[j, ] <- 0 # 整個群組被壓縮為 0 (變數被剔除)
      }
    }

    # 將矩陣轉回向量
    gamma_vec_iter_new <- as.vector(Gamma_mat)

    # --------------------------------------------------------------------------
    # Step 3: U Update
    # --------------------------------------------------------------------------
    u_vec_iter_new <- u_vec_iter + beta_vec_iter_new - gamma_vec_iter_new

    # --------------------------------------------------------------------------
    # Step 4: 收斂檢查
    # --------------------------------------------------------------------------
    res_primal <- sqrt(sum((beta_vec_iter_new - gamma_vec_iter_new)^2))
    res_dual   <- rho.augLagrange.penalty * sqrt(sum((gamma_vec_iter_new - gamma_vec_iter)^2))

    # 更新迭代變數
    beta_vec_iter  <- beta_vec_iter_new
    gamma_vec_iter <- gamma_vec_iter_new
    u_vec_iter     <- u_vec_iter_new

    if (res_primal < tol && res_dual < tol) {

      converged <- TRUE
      break
    }
  }

  # ============================================================================
  # 4. 結果整理
  # ============================================================================
  if (!converged) {

    warning("ADMM did not converge within the maximum number of iterations.")
  }

  # 還原矩陣形式
  beta_free_final <- vec_to_mat(beta_vec_iter)
  gamma_free_final <- vec_to_mat(gamma_vec_iter) # 稀疏解

  # 組合最終 Beta (包含固定的 Identity part)
  final_coef <- reconstruct_full_beta(gamma_free_final)
  beta_coef_smooth <- reconstruct_full_beta(beta_free_final)

  # 設定 row/col names 以利閱讀
  rownames(final_coef) <- colnames(X)
  colnames(final_coef) <- paste0("Index.", 1:dimension)

  results <- list(coef = final_coef,           # Group LASSO 結果 (稀疏)
                  coef.beta = beta_coef_smooth,# 原始優化結果 (平滑)
                  iteration = final_iter,
                  converged = converged,
                  residuals = c(primal = res_primal, dual = res_dual))

  return(results)
}



