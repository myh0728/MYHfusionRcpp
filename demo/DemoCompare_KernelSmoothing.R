library(dplyr)

# ==============================================================================
# 1. 資料生成 (數值驗證用小樣本即可)
# ==============================================================================
set.seed(999)
n <- 50
p <- 2
X <- matrix(rnorm(n * p), n, p)
x_eval <- matrix(rnorm(5 * p), 5, p) # 5個評估點
Y <- sin(X[,1]) + rnorm(n, 0, 0.1)
y_grid <- seq(min(Y), max(Y), length.out = 5)
h <- rep(1.0, p)
w_vec <- runif(n, 0.5, 1.5)

cat("=== 開始數值一致性檢查 (n =", n, ") ===\n\n")

# ==============================================================================
# 2. 執行檢查
# ==============================================================================

# 設定 Kernel 列表 (包含 Gaussian)
kernel_config <- list(
  list(name="Epanechnikov", r_func=K2_Ep_R),
  list(name="Gaussian",     r_func=K2_G_R),
  list(name="K2_Bw",        r_func=K2_Bw_R),
  list(name="K4_Bw",        r_func=K4_Bw_R)
)

# 用來儲存結果的 Data Frame
check_report <- data.frame(
  Kernel = character(),
  Function = character(),
  Weighted = character(),
  Status = character(),     # PASS / FAIL
  Message = character(),    # 差異訊息
  stringsAsFactors = FALSE
)

# 輔助函數：執行比較並回傳結果列表
compare_res <- function(res_rcpp, res_r) {
  # 將兩者都攤平為純向量，避免 attribute mismatch (例如 dim)
  # R 的 matrix 攤平是 column-major，Rcpp arma::mat 轉回 R 也是 column-major，順序會一致
  val1 <- as.vector(res_rcpp)
  val2 <- as.vector(res_r)

  is_eq <- all.equal(val1, val2)

  if (isTRUE(is_eq)) {
    return(list(status = "PASS", msg = ""))
  } else {
    return(list(status = "FAIL", msg = is_eq[1]))
  }
}

for (cfg in kernel_config) {
  k_name <- cfg$name
  k_r    <- cfg$r_func

  for (weighted in c(FALSE, TRUE)) {
    w_arg <- if(weighted) w_vec else NULL
    w_str <- if(weighted) "Weighted" else "Unweighted"

    # --- 1. KDE ---
    val_R    <- KDE_R(X, x_eval, k_r, h, w_arg)
    val_Rcpp <- KDE_rcpp(X, x_eval, h, k_name, w_arg)
    res <- compare_res(val_Rcpp, val_R)

    check_report[nrow(check_report)+1, ] <- list(k_name, "KDE", w_str, res$status, res$msg)

    # --- 2. KNW ---
    val_R    <- KNW_R(Y, X, x_eval, k_r, h, w_arg)
    val_Rcpp <- KNW_rcpp(Y, X, x_eval, h, k_name, w_arg)
    res <- compare_res(val_Rcpp, val_R)

    check_report[nrow(check_report)+1, ] <- list(k_name, "KNW", w_str, res$status, res$msg)

    # --- 3. KNWcdf ---
    val_R    <- KNWcdf_R(Y, y_grid, X, x_eval, k_r, h, w_arg)
    val_Rcpp <- KNWcdf_rcpp(Y, y_grid, X, x_eval, h, k_name, w_arg)
    res <- compare_res(val_Rcpp, val_R)

    check_report[nrow(check_report)+1, ] <- list(k_name, "KNWcdf", w_str, res$status, res$msg)

    # --- 4. CV ---
    val_R    <- CVKNW_R(Y, X, k_r, h, w_arg)
    val_Rcpp <- CVKNW_rcpp(Y, X, h, k_name, w_arg)
    res <- compare_res(val_Rcpp, val_R)

    check_report[nrow(check_report)+1, ] <- list(k_name, "CV", w_str, res$status, res$msg)

    # --- 5. LS ---
    val_R    <- LSKNW_R(Y, X, k_r, h, w_arg)
    val_Rcpp <- LSKNW_rcpp(Y, X, h, k_name, w_arg)
    res <- compare_res(val_Rcpp, val_R)

    check_report[nrow(check_report)+1, ] <- list(k_name, "LS", w_str, res$status, res$msg)
  }
}

# ==============================================================================
# 3. 輸出報告
# ==============================================================================

# 簡化顯示：如果 Message 為空，顯示 "-"
check_report$Message[check_report$Message == ""] <- "-"

if (all(check_report$Status == "PASS")) {
  cat("🎉 完美！所有 Kernel (包含 Gaussian) 與所有函數在 Rcpp 與 Pure R 的結果皆完全一致！\n\n")
} else {
  cat("⚠️ 注意：發現部分結果不一致，請檢查下方報表。\n\n")
}

print(check_report)
































