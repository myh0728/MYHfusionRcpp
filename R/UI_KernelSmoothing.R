KernelSmooth <- function(data = NULL, X.name = NULL, Y.name = NULL,
                         X = NULL, Y = NULL, x = NULL, y = NULL,
                         type = "mean", kernel = "K2.Biweight",
                         bandwidth = NULL, bandwidth.initial = 1,
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

  if (is.null(x)) {

    x <- X

  } else {

    x <- as.matrix(x)
  }

  if (is.null(y)) {

    y <- sort(unique(as.vector(Y)))

  } else {

    y <- as.vector(y)
  }

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_l <- length(y)

  # ============================================================================
  # 1.5 參數檢查 (Parameter Validation)
  #     提早檢查 type，避免進入耗時運算後才報錯
  # ============================================================================
  if (!type %in% c("mean", "distribution")) {

    stop("Invalid 'type'. Must be 'mean' or 'distribution'.")
  }

  # ============================================================================
  # 2. Kernel 名稱對應 (Kernel Name Mapping)
  #    將 R 的命名習慣轉換為 C++ get_kernel_func 能識別的字串
  # ============================================================================
  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  # 如果使用者輸入的名稱在對應表中，就轉換；否則維持原樣
  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  # ============================================================================
  # 3. 頻寬選擇 (Bandwidth Selection)
  #    使用 Unified CV 函數
  # ============================================================================
  if (is.null(bandwidth)) {

    # 定義優化目標函數 (Objective Function)
    LOOCV_Obj <- function(h.log) {

      h_curr <- rep(exp(h.log), number_p)

      if (type == "mean") {

        cv_val <- CVKNW_rcpp(Y = Y, X = X, h = h_curr, kernel = cpp_kernel, w = wi.boot)

      } else {
        # 因為前面已經檢查過 type，這裡 else 一定是 "distribution"
        cv_val <- CVKNWcdf_rcpp(Y = Y, X = X, h = h_curr, kernel = cpp_kernel, w = wi.boot)
      }

      return(cv_val)
    }

    # 執行數值優化
    LOOCV.esti <- nlminb(start = log(bandwidth.initial), objective = LOOCV_Obj)
    bandwidth <- rep(exp(LOOCV.esti$par), number_p)

  } else {

    bandwidth <- rep(bandwidth, length.out = number_p)
  }

  # ============================================================================
  # 4. 估計 (Estimation)
  #    使用 Unified Estimation 函數
  # ============================================================================
  results <- list()

  if (type == "mean") {

    mhat <- KNW_rcpp(Y = Y, X = X, x = x, h = bandwidth, kernel = cpp_kernel, w = wi.boot)

    # 確保輸出格式正確 (column vector)
    if (is.null(dim(mhat))) mhat <- matrix(mhat, ncol = 1)

    rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
    colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")

    results <- list(mean = mhat,
                    eval.covariate = x,
                    bandwidth = bandwidth)

  } else {
    # type == "distribution"
    Fhat <- KNWcdf_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, kernel = cpp_kernel, w = wi.boot)

    rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
    colnames(Fhat) <- paste("y", signif(y, 4), sep = "=")
    colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")

    results <- list(distribution = Fhat,
                    eval.covariate = x,
                    bandwidth = bandwidth)
  }

  return(results)
}






