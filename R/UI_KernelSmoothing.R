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
  # 2. Kernel 名稱對應 (Kernel Name Mapping)
  #    將 R 的命名習慣轉換為 C++ get_kernel_func 能識別的字串
  # ============================================================================
  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  # 如果使用者輸入的名稱在對應表中，就轉換；否則維持原樣 (假設使用者輸入了正確的 C++ 名稱)
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
    if (is.null(dim(mhat))) mhat <- matrix(mhat, ncol = 1)
    rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
    colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")

    results <- list(mean = mhat,
                    eval.covariate = x,
                    bandwidth = bandwidth)

  } else if (type == "distribution") {

    Fhat <- KNWcdf_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, kernel = cpp_kernel, w = wi.boot)
    rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
    colnames(Fhat) <- paste("y", signif(y, 4), sep = "=")
    colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")

    results <- list(distribution = Fhat,
                    eval.covariate = x,
                    bandwidth = bandwidth)
  } else {

    stop("Invalid 'type'. Must be 'mean' or 'distribution'.")
  }

  return(results)
}




KernelSmooth.old <- function(data = NULL, X.name = NULL, Y.name = NULL,
                             X = NULL, Y = NULL, x = NULL, y = NULL,
                             type = "mean", kernel = "K2.Biweight",
                             bandwidth = NULL, bandwidth.initial = 1,
                             wi.boot = NULL)
{
  if (!is.null(data))
  {
    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])
  }else
  {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  if (is.null(x))
  {
    x <- X
  }else
  {
    x <- as.matrix(x)
  }

  if (is.null(y))
  {
    y <- sort(unique(as.vector(Y)))
  }else
  {
    y <- as.vector(y)
  }

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_l <- length(y)

  if (type == "mean")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNW_K2Ep_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        mhat <- KNW_K2Ep_rcpp(Y = Y, X = X, x = x, h = bandwidth)
      }else
      {
        mhat <- KNW_K2Ep_w_rcpp(Y = Y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(mean = mhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNW_K2Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        mhat <- KNW_K2Bw_rcpp(Y = Y, X = X, x = x, h = bandwidth)
      }else
      {
        mhat <- KNW_K2Bw_w_rcpp(Y = Y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(mean = mhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNW_K4Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        mhat <- KNW_K4Bw_rcpp(Y = Y, X = X, x = x, h = bandwidth)
      }else
      {
        mhat <- KNW_K4Bw_w_rcpp(Y = Y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(mean = mhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }
  }

  if (type == "distribution")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        Fhat <- KNWcdf_K2Ep_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth)
      }else
      {
        Fhat <- KNWcdf_K2Ep_w_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(Fhat) <- y
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(distribution = Fhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNWcdf_K2Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        Fhat <- KNWcdf_K2Bw_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth)
      }else
      {
        Fhat <- KNWcdf_K2Bw_w_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(Fhat) <- y
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(distribution = Fhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNWcdf_K4Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        Fhat <- KNWcdf_K4Bw_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth)
      }else
      {
        Fhat <- KNWcdf_K4Bw_w_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(Fhat) <- y
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(distribution = Fhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }
  }

  return(results)
}


