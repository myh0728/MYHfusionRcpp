# ==============================================================================
# KernelSmooth Wrapper Function Test Script
#
# 說明：
# 1. 測試 KernelSmooth 主函數對各種參數組合的支援度
# 2. 涵蓋 Mean Regression 與 Distribution Estimation
# 3. 測試自動頻寬選擇 (Auto Bandwidth) 與固定頻寬 (Fixed Bandwidth)
# ==============================================================================

# ==============================================================================
# 1. 資料生成 (Data Generation)
# ==============================================================================
set.seed(2025)

n <- 100   # 樣本數
p <- 3     # Covariate 維度
k <- 10    # 評估點數量 (Evaluation Points)

# 生成資料
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- matrix(rnorm(k * p), nrow = k, ncol = p) # 評估點
w <- rexp(n = n, rate = 1)                    # 隨機權重
beta0 <- c(1, rep(0.5, p - 1))                # 真實係數

# 生成 Y (帶有雜訊)
Y <- as.vector(sin(X %*% beta0)) + rnorm(n, mean = 0, sd = 0.2)

# 生成 Distribution 用的 Grid
y <- seq(min(Y), max(Y), length.out = 10)

cat("=== 資料準備完成: n =", n, ", p =", p, ", k =", k, "===\n\n")

# ==============================================================================
# 2. 測試：均值迴歸 (Mean Regression)
# type = "mean"
# ==============================================================================
cat("--- [Test Group 1] Mean Regression ---\n")

# Case 1.1: 固定頻寬 (Fixed Bandwidth), 無權重
cat("\n1.1 [Mean] K2.Epanechnikov, Fixed Bandwidth (0.5)\n")
res_mean_1 <- KernelSmooth(X = X, Y = Y, x = x,
                           type = "mean",
                           kernel = "K2.Epanechnikov",
                           bandwidth = 0.5)
print(str(res_mean_1)) # 檢查結構

# Case 1.2: 自動頻寬 (Auto Bandwidth), 無權重
# (不傳入 bandwidth，程式應自動執行 CV)
cat("\n1.2 [Mean] K2.Epanechnikov, Auto Bandwidth\n")
res_mean_2 <- KernelSmooth(X = X, Y = Y, x = x,
                           type = "mean",
                           kernel = "K2.Epanechnikov")
cat("Selected Bandwidth:", paste(round(res_mean_2$bandwidth, 4), collapse=", "), "\n")

# Case 1.3: 固定頻寬, 有權重 (Weighted)
cat("\n1.3 [Mean] K2.Epanechnikov, Fixed BW, Weighted\n")
res_mean_3 <- KernelSmooth(X = X, Y = Y, x = x,
                           type = "mean",
                           kernel = "K2.Epanechnikov",
                           bandwidth = 0.5,
                           wi.boot = w)
# 簡單檢查是否有算出數值
cat("Output head:", head(res_mean_3$mean, 3), "\n")

# Case 1.4: 自動頻寬, 有權重 (Weighted Auto)
cat("\n1.4 [Mean] K2.Epanechnikov, Auto BW, Weighted\n")
res_mean_4 <- KernelSmooth(X = X, Y = Y, x = x,
                           type = "mean",
                           kernel = "K2.Epanechnikov",
                           wi.boot = w)
cat("Selected Bandwidth (Weighted):", paste(round(res_mean_4$bandwidth, 4), collapse=", "), "\n")


# ==============================================================================
# 3. 測試：分佈估計 (Distribution Estimation)
# type = "distribution"
# ==============================================================================
cat("\n--- [Test Group 2] Distribution Estimation ---\n")

# Case 2.1: K2.Epanechnikov, 固定頻寬
cat("\n2.1 [Dist] K2.Epanechnikov, Fixed BW (0.5)\n")
res_dist_1 <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                           type = "distribution",
                           kernel = "K2.Epanechnikov",
                           bandwidth = 0.5)
# 檢查 CDF 輸出維度 (應為 k x length(y))
cat("Distribution Matrix Dim:", dim(res_dist_1$distribution), "\n")

# Case 2.2: K2.Epanechnikov, 自動頻寬, 有權重
cat("\n2.2 [Dist] K2.Epanechnikov, Auto BW, Weighted\n")
res_dist_2 <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                           type = "distribution",
                           kernel = "K2.Epanechnikov",
                           wi.boot = w)
cat("Selected Bandwidth:", paste(round(res_dist_2$bandwidth, 4), collapse=", "), "\n")

# Case 2.3: K2.Biweight, 自動頻寬 (Implied)
cat("\n2.3 [Dist] K2.Biweight, Auto BW\n")
res_dist_3 <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                           type = "distribution",
                           kernel = "K2.Biweight")
cat("K2.Bw Result check (first row sum approx 1?):", sum(diff(y) * head(res_dist_3$distribution[1,], -1)), "(Rough check)\n")

# Case 2.4: K4.Biweight, 自動頻寬 (Implied)
cat("\n2.4 [Dist] K4.Biweight, Auto BW\n")
res_dist_4 <- KernelSmooth(X = X, Y = Y, x = x, y = y,
                           type = "distribution",
                           kernel = "K4.Biweight")
cat("Selected Bandwidth:", paste(round(res_dist_4$bandwidth, 4), collapse=", "), "\n")


# ==============================================================================
# 4. (Optional) Gaussian 測試
# 如果您的 map 有加入 Gaussian，可以取消註解測試
# ==============================================================================
cat("\n--- [Test Group 3] Gaussian Kernel ---\n")
res_gauss <- KernelSmooth(X = X, Y = Y, x = x, type = "mean", kernel = "Gaussian")
cat("Gaussian Auto Bandwidth:", paste(round(res_gauss$bandwidth, 4), collapse=", "), "\n")

cat("\n所有測試執行完畢。\n")





