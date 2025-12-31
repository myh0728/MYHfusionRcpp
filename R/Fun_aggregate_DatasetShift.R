GLMcombineADavar <- function(name.par,
                             MLE.hessian,
                             Psi.square,
                             Psi.diff.beta,
                             Psi.diff.phi = NULL,
                             Psi.diff.theta = NULL,
                             kappa = NULL,
                             avar.phi = NULL,
                             avar.working = NULL,
                             eps.inv = 1e-7) {

  # --- 1. 確定維度與參數名稱 ---

  n_beta <- length(name.par)      # MLE 參數個數
  n_eta  <- nrow(Psi.square)      # 限制式/Lagrange 參數個數 (number_m)

  # 建立所有參數的名稱向量
  names_all <- name.par

  # 處理 Phi (Auxiliary parameters)
  has_phi <- !is.null(kappa)
  n_phi <- 0
  if (has_phi) {
    n_phi <- n_eta # 通常 phi 的維度對應 eta
    names_all <- c(names_all, paste0("phi", 1:n_phi))
  }

  # 處理 Theta (其他參數)
  has_theta <- !is.null(Psi.diff.theta)
  n_theta <- 0
  if (has_theta) {
    # 假設輸入是向量，需計算列數
    n_theta <- length(Psi.diff.theta) / n_eta
    names_all <- c(names_all, paste0("theta", 1:n_theta))
  }

  # 加入 Eta
  names_all <- c(names_all, paste0("eta", 1:n_eta))

  n_total <- length(names_all)

  # --- 2. 定義區塊索引 (Block Indices) ---
  # 使用數值索引比字串索引更高效且不易出錯

  idx_beta  <- 1:n_beta
  current_idx <- n_beta

  idx_phi <- NULL
  if (has_phi) {
    idx_phi <- (current_idx + 1):(current_idx + n_phi)
    current_idx <- current_idx + n_phi
  }

  idx_theta <- NULL
  if (has_theta) {
    idx_theta <- (current_idx + 1):(current_idx + n_theta)
    current_idx <- current_idx + n_theta
  }

  idx_eta <- (current_idx + 1):(current_idx + n_eta)

  # --- 3. 初始化矩陣 (Meat S & Bread J) ---

  avar.S <- matrix(0, nrow = n_total, ncol = n_total, dimnames = list(names_all, names_all))
  J.V    <- matrix(0, nrow = n_total, ncol = n_total, dimnames = list(names_all, names_all))

  # --- 4. 填充矩陣區塊 ---

  # Block: beta (MLE 部分)
  # Hessian 通常是負的 Information Matrix，Variance 計算需要 Information
  avar.S[idx_beta, idx_beta] <- -MLE.hessian
  J.V[idx_beta, idx_beta]    <- -MLE.hessian

  # Block: eta (Constraint Variance 部分)
  avar.S[idx_eta, idx_eta] <- Psi.square
  J.V[idx_eta, idx_eta]    <- -Psi.square

  # Block: Cross terms eta-beta (Psi 對 beta 的微分)
  # 確保轉為矩陣格式 (n_eta x n_beta)
  mat_psi_diff_beta <- matrix(as.vector(Psi.diff.beta), nrow = n_eta, ncol = n_beta)
  J.V[idx_eta, idx_beta] <- mat_psi_diff_beta
  J.V[idx_beta, idx_eta] <- t(mat_psi_diff_beta)

  # Block: Phi (如果存在)
  if (has_phi) {
    if (is.null(avar.phi)) stop("kappa provided but avar.phi is NULL")

    # 來自 C++ 的 rcond_rcpp 檢查條件數
    if (is.null(avar.working)) {
      # 情境 1: 沒有 working variance
      cond_num <- tryCatch(rcond_rcpp(avar.phi), error = function(e) 0)

      avar.inv <- if (cond_num > eps.inv) {
        inv_sympd_rcpp(avar.phi)
      } else {
        diag(n_phi) # 若奇異則退化為單位矩陣
      }

      block_val <- kappa * avar.inv
      avar.S[idx_phi, idx_phi] <- block_val
      J.V[idx_phi, idx_phi]    <- block_val

    } else {
      # 情境 2: 有 working variance
      cond_num <- tryCatch(rcond_rcpp(avar.working), error = function(e) 0)

      avar.working.inv <- if (cond_num > eps.inv) {
        inv_sympd_rcpp(avar.working)
      } else {
        diag(n_phi)
      }

      # Sandwich construction for phi block
      avar.S[idx_phi, idx_phi] <- kappa * (avar.working.inv %*% avar.phi %*% avar.working.inv)
      J.V[idx_phi, idx_phi]    <- kappa * avar.working.inv
    }

    # Cross terms eta-Phi
    mat_psi_diff_phi <- matrix(as.vector(Psi.diff.phi), nrow = n_eta, ncol = n_phi)
    J.V[idx_eta, idx_phi] <- mat_psi_diff_phi
    J.V[idx_phi, idx_eta] <- t(mat_psi_diff_phi)
  }

  # Block: Theta (如果存在)
  if (has_theta) {
    mat_psi_diff_theta <- matrix(as.vector(Psi.diff.theta), nrow = n_eta, ncol = n_theta)
    J.V[idx_eta, idx_theta] <- mat_psi_diff_theta
    J.V[idx_theta, idx_eta] <- t(mat_psi_diff_theta)
  }

  # --- 5. 計算漸近變異數 (Sandwich Estimator) ---

  # 計算 J 的反矩陣
  cond_J <- tryCatch(rcond_rcpp(J.V), error = function(e) 0)

  if (cond_J > eps.inv) {
    # 矩陣條件良好，直接求解
    J.V.inv <- solve_rcpp(J.V, diag(n_total))
  } else {
    # 矩陣接近奇異，使用您定義的廣義逆矩陣函數
    J.V.inv <- ginv.sym.eigen(M = J.V, eps.inv = eps.inv)
  }
  dimnames(J.V.inv) <- dimnames(J.V)

  # Sandwich Formula: Cov = J^-1 * S * J^-T
  asy.Cov <- J.V.inv %*% avar.S %*% t(J.V.inv)
  dimnames(asy.Cov) <- list(names_all, names_all)

  # 提取主要參數 (Beta) 的變異數
  asy.Cov.par <- asy.Cov[idx_beta, idx_beta, drop = FALSE]

  # --- 6. 回傳結果 ---
  return(list(
    avar.S      = avar.S,
    J.V         = J.V,
    J.V.inv     = J.V.inv,
    asy.Cov     = asy.Cov,
    asy.Cov.par = asy.Cov.par
  ))
}

