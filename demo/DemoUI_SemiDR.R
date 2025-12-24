# ==============================================================================
# Comprehensive Test Script for Semi-parametric Index Models
#
# 涵蓋函數:
# 1. SemiSI (Single Index Model)
# 2. SemiMI (Multi-Index Model)
# 3. SemiSI.adtLASSO.ADMM (Adaptive LASSO Variable Selection via ADMM)
# ==============================================================================

# ==============================================================================
# 1. 資料生成 (Data Generation)
# ==============================================================================
set.seed(2025)

n <- 200    # 樣本數
p <- 5      # Covariate 維度
k <- 10     # 評估點數量

# 生成 X 與 評估點 x
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- matrix(rnorm(k * p), nrow = k, ncol = p)

# 生成隨機權重
w <- rexp(n = n, rate = 1)

# 設定真實係數 beta0 (稀疏結構: 後面幾個係數為 0)
# beta0: [1, 0.5, 0.5, 0, 0]
beta0 <- c(1, 0.5, 0.5, 0, 0)

# 生成 Y (Single Index 結構)
# Y = sin(X * beta) + noise
Y <- as.vector(sin(X %*% beta0)) + rnorm(n, mean = 0, sd = 0.2)

# 生成 Distribution 用的 Grid
y_grid <- seq(min(Y), max(Y), length.out = 10)

cat("=== 資料準備完成: n =", n, ", p =", p, "===\n\n")

# ==============================================================================
# 2. 測試 SemiSI (Single Index Model)
# ==============================================================================
cat("--- [Test Group 1] SemiSI ---\n")

kernels <- c("K2.Epanechnikov", "K2.Biweight", "K4.Biweight")
types <- c("mean", "distribution")

for (kern in kernels) {
  for (typ in types) {
    cat(sprintf("Testing SemiSI | Type: %s | Kernel: %s\n", typ, kern))

    # 2.1 無權重 (Unweighted)
    # y 參數只在 distribution type 時需要，但在 R 中多傳通常沒問題(若函數有處理)
    # 這裡明確分開呼叫比較嚴謹
    if (typ == "mean") {
      res <- SemiSI(X = X, Y = Y, type = typ, kernel = kern)
    } else {
      res <- SemiSI(X = X, Y = Y, type = typ, kernel = kern) # SemiSI 內部預設會處理 y grid
    }

    # 2.2 有權重 (Weighted)
    if (typ == "mean") {
      res_w <- SemiSI(X = X, Y = Y, type = typ, kernel = kern, wi.boot = w)
    } else {
      res_w <- SemiSI(X = X, Y = Y, type = typ, kernel = kern, wi.boot = w)
    }
  }
}
cat("SemiSI 測試完成。\n\n")

# ==============================================================================
# 3. 測試 SemiMI (Multi-Index Model)
# ==============================================================================
cat("--- [Test Group 2] SemiMI (Dimension = 2) ---\n")

dim_mi <- 2

for (kern in kernels) {
  for (typ in types) {
    cat(sprintf("Testing SemiMI | Type: %s | Kernel: %s\n", typ, kern))

    # 3.1 無權重
    res <- SemiMI(X = X, Y = Y, type = typ, dimension = dim_mi, kernel = kern)

    # 3.2 有權重
    res_w <- SemiMI(X = X, Y = Y, type = typ, dimension = dim_mi, kernel = kern, wi.boot = w)
  }
}
cat("SemiMI 測試完成。\n\n")

# ==============================================================================
# 4. 測試 SemiSI.adtLASSO.ADMM (Variable Selection)
# ==============================================================================
cat("--- [Test Group 3] SemiSI.adtLASSO.ADMM ---\n")

# 4.0 準備工作：先跑一次標準 SemiSI 取得初始值與 Adaptive Weights
# -------------------------------------------------------------
cat("Step 0: Pre-calculation using standard SemiSI...\n")
pre_fit <- SemiSI(X = X, Y = Y, type = "mean", kernel = "K2.Biweight")

# 提取初始值
beta_init_vec <- pre_fit$coef        # 這是 SemiSI 估計出來的 beta
bandwidth_val <- pre_fit$bandwidth   # 這是 CV 選出來的頻寬

# 設定 Adaptive Weights (通常是 1 / |beta_initial|)
# 避免分母為 0
weight_vec <- 1 / (abs(beta_init_vec) + 1e-5)

lambda_val <- 0.1 # 懲罰力度 (範例值)

cat("Step 1: Running ADMM for Mean Regression...\n")

# 4.1 測試 Mean Regression ADMM
admm_mean <- SemiSI.adtLASSO.ADMM(
  X = X, Y = Y,
  type = "mean",
  kernel = "K2.Biweight",
  beta.initial = beta_init_vec,
  bandwidth = bandwidth_val,
  weight.beta = weight_vec,
  lambda.L1.penalty = lambda_val,
  rho.augLagrange.penalty = 1,
  iter.max = 50
)

cat("ADMM (Mean) Result Beta:", round(admm_mean$coef, 4), "\n")
cat("Converged:", admm_mean$converged, "\n\n")


cat("Step 2: Running ADMM for Distribution Estimation...\n")

# 4.2 測試 Distribution ADMM
# 注意：通常 distribution 的 bandwidth 會跟 mean 不一樣，
# 這裡為了演示方便沿用，實際使用建議先跑一次 type="distribution" 的 SemiSI
admm_dist <- SemiSI.adtLASSO.ADMM(
  X = X, Y = Y,
  type = "distribution",
  kernel = "K2.Biweight",
  beta.initial = beta_init_vec,
  bandwidth = bandwidth_val,
  weight.beta = weight_vec,
  lambda.L1.penalty = lambda_val,
  rho.augLagrange.penalty = 1,
  iter.max = 50
)

cat("ADMM (Dist) Result Beta:", round(admm_dist$coef, 4), "\n")
cat("Converged:", admm_dist$converged, "\n")

cat("\n所有測試執行完畢。\n")






# ==============================================================================
# 1. 資料生成 (Data Generation for Multi-Index Model)
# ==============================================================================
set.seed(2025)

n <- 300    # 樣本數
p <- 6      # Covariate 總維度
d <- 2      # Indices 維度 (Dimension)

# 生成 X
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# 設定真實係數 Beta (p x d)
# 識別限制：前 d x d (前 2x2) 必須是單位矩陣
# 變數選取測試：設定第 3 個變數有作用，第 4, 5, 6 個變數為 0 (無作用)
true_beta <- matrix(0, nrow = p, ncol = d)

# 固定部分 (Identity)
true_beta[1:d, 1:d] <- diag(d)

# 自由部分 (Free parameters)
true_beta[3, ] <- c(0.5, -0.5) # 第3個變數在兩個 Index 都有係數
# true_beta[4:6, ] 保持為 0 (稀疏)

cat("真實 Beta 結構:\n")
print(true_beta)

# 生成 Y (Multi-Index 結構)
# Index 1: X * beta[,1]
# Index 2: X * beta[,2]
idx1 <- X %*% true_beta[, 1]
idx2 <- X %*% true_beta[, 2]

# Y 與兩個 Index 非線性相關
Y <- as.vector(sin(idx1) + (idx2^2 - 1)) + rnorm(n, mean = 0, sd = 0.2)

# Distribution 用的 Grid
y_grid <- seq(min(Y), max(Y), length.out = 10)

cat("\n=== 資料準備完成: n =", n, ", p =", p, ", d =", d, "===\n")

# ==============================================================================
# 4. 測試 SemiMI.adtLASSO.ADMM
# ==============================================================================
cat("\n--- [Test Group 4] SemiMI.adtLASSO.ADMM ---\n")

# 4.0 準備工作：先跑一次標準 SemiMI 取得初始值與 Adaptive Weights
# -------------------------------------------------------------
cat("Step 0: Pre-calculation using standard SemiMI...\n")

# 注意：這裡使用 dimension = 2
pre_fit <- SemiMI(
  X = X, Y = Y,
  type = "mean",
  dimension = d,
  kernel = "K2.Biweight"
)

# 提取初始值 (p x d 矩陣)
beta_init_mat <- pre_fit$coef
bandwidth_val <- pre_fit$bandwidth

cat("Initial Beta (from SemiMI):\n")
print(round(beta_init_mat, 4))

# 設定 Adaptive Weights (Group Weights)
# 對於每個變數 j (每一列)，計算其 L2 Norm: ||beta_j||_2
# weight_j = 1 / ||beta_j||_2
beta_row_norms <- apply(beta_init_mat, 1, function(v) sqrt(sum(v^2)))

# 避免分母為 0
weight_vec <- 1 / (beta_row_norms + 1e-5)

# (可選) 強制將前 d 個變數的權重設為 0 (雖然 ADMM 函數內部會自動忽略它們)
weight_vec[1:d] <- 0

cat("Adaptive Weights:\n")
print(round(weight_vec, 4))

lambda_val <- 0.5 # Group LASSO 懲罰力度

# -------------------------------------------------------------
cat("\nStep 1: Running ADMM for Mean Regression (Group LASSO)...\n")

# 4.1 測試 Mean Regression ADMM
admm_mean <- SemiMI.adtLASSO.ADMM(
  X = X, Y = Y,
  type = "mean",
  dimension = d,
  kernel = "K2.Biweight",
  beta.initial = beta_init_mat,   # 傳入矩陣
  bandwidth = bandwidth_val,      # 傳入向量
  weight.beta = weight_vec,       # 傳入向量 (長度 p)
  lambda.L1.penalty = lambda_val,
  rho.augLagrange.penalty = 1,
  iter.max = 50
)

cat("ADMM (Mean) Result Beta (Sparse):\n")
print(round(admm_mean$coef, 4))
cat("Converged:", admm_mean$converged, "\n")

# 檢查變數是否被正確剔除 (第 4-6 列應全為 0)
is_zero <- apply(admm_mean$coef, 1, function(x) all(abs(x) < 1e-4))
cat("Variables selected (False means zero/removed):\n")
print(!is_zero)


# -------------------------------------------------------------
cat("\nStep 2: Running ADMM for Distribution Estimation (Group LASSO)...\n")

# 4.2 測試 Distribution ADMM
admm_dist <- SemiMI.adtLASSO.ADMM(
  X = X, Y = Y,
  type = "distribution",
  dimension = d,
  kernel = "K2.Biweight",
  beta.initial = beta_init_mat,
  bandwidth = bandwidth_val,
  weight.beta = weight_vec,
  lambda.L1.penalty = lambda_val,
  rho.augLagrange.penalty = 1,
  iter.max = 100
)

cat("ADMM (Dist) Result Beta (Sparse):\n")
print(round(admm_dist$coef, 4))
cat("Converged:", admm_dist$converged, "\n")

cat("\n所有測試執行完畢。\n")

