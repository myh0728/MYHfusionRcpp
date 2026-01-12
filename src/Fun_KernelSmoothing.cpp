#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include "head_KernelFunctions.h"

// kernel functions

// =========================================================
//  定義 Kernel 函數 (Kernels available in R)
// =========================================================

// [[Rcpp::export]]
double K2_Ep_rcpp(const double & u) {

  return K2_Ep_inline(u);
}

// [[Rcpp::export]]
double K2_Ep_diff1_rcpp(const double & u) {

  return K2_Ep_diff1_inline(u);
}

// [[Rcpp::export]]
double K2_Ep_diff2_rcpp(const double & u) {

  return K2_Ep_diff2_inline(u);
}

// [[Rcpp::export]]
double K2_Bw_rcpp(const double & u) {

  return K2_Bw_inline(u);
}

// [[Rcpp::export]]
double K2_Bw_diff1_rcpp(const double & u) {

  return K2_Bw_diff1_inline(u);
}

// [[Rcpp::export]]
double K2_Bw_diff2_rcpp(const double & u) {

  return K2_Bw_diff2_inline(u);
}

// [[Rcpp::export]]
double K4_Bw_rcpp(const double & u) {

  return K4_Bw_inline(u);
}

// [[Rcpp::export]]
double K4_Bw_diff1_rcpp(const double & u) {

  return K4_Bw_diff1_inline(u);
}

// [[Rcpp::export]]
double K4_Bw_diff2_rcpp(const double & u) {

  return K4_Bw_diff2_inline(u);
}

// [[Rcpp::export]]
double K2_G_rcpp(const double & u) {

  return K2_G_inline(u);
}

// [[Rcpp::export]]
double K2_G_diff1_rcpp(const double & u) {

  return K2_G_diff1_inline(u);
}

// [[Rcpp::export]]
double K2_G_diff2_rcpp(const double & u) {

  return K2_G_diff2_inline(u);
}

// =========================================================
//  製作 kernel 選擇器
// =========================================================

// 輔助函數：根據字串名稱回傳對應的函數指標
KernelFunc get_kernel_func(std::string kernel_name) {

  if (kernel_name == "K2_Ep") {

    return K2_Ep_inline;

  } else if (kernel_name == "K2_G") {

    return K2_G_inline;

  } else if (kernel_name == "K2_Bw") {

    return K2_Bw_inline;

  } else if (kernel_name == "K4_Bw") {

    return K4_Bw_inline;

  } else {

    stop("Unknown kernel name."); // Rcpp 的報錯機制
    return NULL;
  }
}

// =========================================================
//  整合型 KDE 函數 (Unified KDE)
//  Kernel density estimation
// =========================================================

// [[Rcpp::export]]
arma::vec KDE_rcpp(const arma::mat & X,
                   const arma::mat & x,
                   const arma::vec & h,
                   std::string kernel = "K2_Ep",
                   Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows; // 樣本數
  const arma::uword n_p = X.n_cols; // 維度
  const arma::uword n_k = x.n_rows; // 評估點數

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 初始化結果向量
  arma::vec Dhat(n_k, arma::fill::zeros);

  // 檢查是否有權重
  bool has_weights = w.isNotNull();
  arma::vec weights;
  if (has_weights) {
    weights = as<arma::vec>(w);
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");
  }

  // =========================================================
  // OpenMP 平行化策略：
  // 針對每一個評估點 k 進行平行運算
  // =========================================================

  if (has_weights) {
    // === 情境 A：有權重 ===
#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < n_k; ++k) {

      // 每個執行緒獨立計算該 k 點的密度總和
      double density_sum = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        double Kik_h = 1.0;

        // 內層維度迴圈 (Product Kernel)
        for (arma::uword l = 0; l < n_p; ++l) {
          // [優化] 直接讀取元素，不建立 Xrow_i 暫存向量
          // Armadillo: X(row, col)
          double u = (X(i, l) - x(k, l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        density_sum += Kik_h * weights(i);
      }
      Dhat(k) = density_sum;
    }

  } else {
    // === 情境 B：無權重 ===
#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < n_k; ++k) {

      double density_sum = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {
          // [優化] 直接讀取元素
          double u = (X(i, l) - x(k, l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        density_sum += Kik_h;
      }
      Dhat(k) = density_sum;
    }
  }

  // 最後平均
  Dhat /= n_n;

  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_rcpp0(const arma::mat & X,
                    const arma::mat & x,
                    const arma::vec & h,
                    std::string kernel = "K2_Ep", // 使用者選擇 Kernel
                    Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) { // 選填參數

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;

  // 取得使用者選擇的 Kernel 函數指標
  KernelFunc k_func = get_kernel_func(kernel);

  // 初始化結果向量
  arma::vec Dhat(n_k, arma::fill::zeros);

  // 預先宣告變數
  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);

  // 判斷是否有傳入權重 w
  if (w.isNotNull()) {

    // ============================
    //  情境 A：有權重 (Weighted)
    // ============================

    arma::vec weights = as<arma::vec>(w);

    // 防呆：檢查權重長度
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");

    // 主迴圈

    for (arma::uword k = 0; k < n_k; ++k) {

      xrow_k = x.row(k).t();

      for (arma::uword i = 0; i < n_n; ++i) {

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - xrow_k(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        Dhat(k) += Kik_h * weights(i);
      }
    }

  } else {

    // ============================
    //  情境 B：無權重 (Unweighted)
    // ============================

    for (arma::uword k = 0; k < n_k; ++k) {

      xrow_k = x.row(k).t();

      for (arma::uword i = 0; i < n_n; ++i) {

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - xrow_k(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        Dhat(k) += Kik_h;
      }
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// =========================================================
//  整合型 KNW 函數 (Unified KNW)
//  Nadaraya-Watson Estimator
// =========================================================

// [[Rcpp::export]]
arma::vec KNW_rcpp(const arma::vec & Y,
                   const arma::mat & X,
                   const arma::mat & x,
                   const arma::vec & h,
                   std::string kernel = "K2_Ep",
                   Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows; // 樣本數
  const arma::uword n_p = X.n_cols; // 維度
  const arma::uword n_k = x.n_rows; // 評估點數

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 初始化結果向量
  arma::vec mhat(n_k, arma::fill::zeros);

  // 檢查是否有權重
  bool has_weights = w.isNotNull();
  arma::vec weights;
  if (has_weights) {
    weights = as<arma::vec>(w);
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");
  }

  // =========================================================
  // OpenMP 平行化策略：針對每一個評估點 k 進行平行運算
  // =========================================================

  if (has_weights) {
    // === 情境 A：有權重 ===
#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < n_k; ++k) {

      double num_k = 0.0; // 分子累加 (Thread-private)
      double den_k = 0.0; // 分母累加 (Thread-private)

      for (arma::uword i = 0; i < n_n; ++i) {
        // 優化：權重為 0 直接跳過
        if (weights(i) == 0) continue;

        double Kik_h = 1.0;

        // 內層維度迴圈
        for (arma::uword l = 0; l < n_p; ++l) {
          // [優化] 直接讀取元素，避免建立暫存向量
          double u = (X(i, l) - x(k, l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        double weighted_K = Kik_h * weights(i);
        den_k += weighted_K;
        num_k += weighted_K * Y(i);
      }

      // 計算該點的最終估計值
      if (std::abs(den_k) > 1e-10) {
        mhat(k) = num_k / den_k;
      }
    }

  } else {
    // === 情境 B：無權重 ===
#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < n_k; ++k) {

      double num_k = 0.0;
      double den_k = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {
          // [優化] 直接讀取元素
          double u = (X(i, l) - x(k, l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        den_k += Kik_h;
        num_k += Kik_h * Y(i);
      }

      if (std::abs(den_k) > 1e-10) {
        mhat(k) = num_k / den_k;
      }
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_rcpp0(const arma::vec & Y,
                    const arma::mat & X,
                    const arma::mat & x,
                    const arma::vec & h,
                    std::string kernel = "K2_Ep",
                    Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;

  // 取得使用者選擇的 Kernel 函數指標
  KernelFunc k_func = get_kernel_func(kernel);

  // 初始化為 0
  arma::vec Dhat(n_k, arma::fill::zeros); // 分母 (Denominator)
  arma::vec Nhat(n_k, arma::fill::zeros); // 分子 (Numerator)
  arma::vec mhat(n_k, arma::fill::zeros); // 結果

  // 預先宣告變數
  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);

  // 判斷是否有傳入權重 w
  if (w.isNotNull()) {

    // ============================
    //  情境 A：有權重 (Weighted)
    // ============================

    arma::vec weights = as<arma::vec>(w);

    // 防呆：檢查權重長度
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");

    // 主迴圈
    for (arma::uword k = 0; k < n_k; ++k) {

      xrow_k = x.row(k).t();

      for (arma::uword i = 0; i < n_n; ++i) {
        // 如果權重是 0，直接跳過計算，省時間 (Optional optimization)
        if (weights(i) == 0) continue;

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - xrow_k(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        double weighted_K = Kik_h * weights(i);
        Dhat(k) += weighted_K;
        Nhat(k) += weighted_K * Y(i);
      }
    }
  } else {

    // ============================
    //  情境 B：無權重 (Unweighted)
    // ============================

    for (arma::uword k = 0; k < n_k; ++k) {

      xrow_k = x.row(k).t();

      for (arma::uword i = 0; i < n_n; ++i) {

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - xrow_k(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        Dhat(k) += Kik_h;
        Nhat(k) += Kik_h * Y(i);
      }
    }
  }

  // -------------------------------------------------------
  // 計算最終估計值 (Division)
  // -------------------------------------------------------

  for (arma::uword k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);

    } else {
      mhat(k) = 0;
    }
  }

  return mhat;
}

// ===========================================================================
//  整合型 KNW CDF 函數 (Unified KNWcdf)
//  Nadaraya-Watson Estimator for conditional cumulative distribution function
// ===========================================================================

// [[Rcpp::export]]
arma::mat KNWcdf_rcpp(const arma::vec & Y,
                      const arma::vec & y,  // CDF 的評估點 (grid for Y)
                      const arma::mat & X,
                      const arma::mat & x,  // Covariates 的評估點 (grid for X)
                      const arma::vec & h,
                      std::string kernel = "K2_Ep",
                      Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows; // 樣本數
  const arma::uword n_p = X.n_cols; // X 維度
  const arma::uword n_k = x.n_rows; // X grid 數 (Rows of result)
  const arma::uword n_l = y.n_elem; // Y grid 數 (Cols of result)

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 初始化結果矩陣
  arma::mat mhat(n_k, n_l, arma::fill::zeros);

  // 檢查是否有權重
  bool has_weights = w.isNotNull();
  arma::vec weights;
  if (has_weights) {
    weights = as<arma::vec>(w);
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");
  }

  // =========================================================
  // OpenMP 平行化策略：針對每一個 Covariate Grid k 進行平行運算
  // =========================================================

  if (has_weights) {
    // === 情境 A：有權重 ===
#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < n_k; ++k) {

      // [記憶體優化]：不需要巨大的 Nhat 矩陣
      // 每個執行緒只需維護當前 k 的分子向量 (Thread-private)
      arma::rowvec num_row_k(n_l, arma::fill::zeros);
      double den_k = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {
        if (weights(i) == 0) continue;

        double Kik_h = 1.0;

        // Kernel 計算 (直接讀取元素，不配置向量)
        for (arma::uword j = 0; j < n_p; ++j) {
          double u = (X(i, j) - x(k, j)) / h(j);
          Kik_h *= k_func(u) / h(j);
        }

        double val = Kik_h * weights(i);
        den_k += val;

        // 針對 CDF Grid y 進行累加
        // 這裡無法輕易優化，因為 Y(i) 和 y(l) 的關係是不確定的
        for (arma::uword l = 0; l < n_l; ++l) {
          if (Y(i) <= y(l)) {
            num_row_k(l) += val;
          }
        }
      }

      // 計算完該 row (k) 的所有樣本後，直接做除法寫入結果
      if (std::abs(den_k) > 1e-10) {
        mhat.row(k) = num_row_k / den_k;
      }
    }

  } else {
    // === 情境 B：無權重 ===
#pragma omp parallel for schedule(static)
    for (arma::uword k = 0; k < n_k; ++k) {

      arma::rowvec num_row_k(n_l, arma::fill::zeros);
      double den_k = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        double Kik_h = 1.0;

        for (arma::uword j = 0; j < n_p; ++j) {
          double u = (X(i, j) - x(k, j)) / h(j);
          Kik_h *= k_func(u) / h(j);
        }

        den_k += Kik_h;

        for (arma::uword l = 0; l < n_l; ++l) {
          if (Y(i) <= y(l)) {
            num_row_k(l) += Kik_h;
          }
        }
      }

      if (std::abs(den_k) > 1e-10) {
        mhat.row(k) = num_row_k / den_k;
      }
    }
  }

  // 限制範圍 [0, 1] (Armadillo 的 clamp 是平行化友善的)
  mhat.clamp(0.0, 1.0);

  return mhat;
}

// [[Rcpp::export]]
arma::mat KNWcdf_rcpp0(const arma::vec & Y,
                       const arma::vec & y,  // CDF 的評估點 (grid for Y)
                       const arma::mat & X,
                       const arma::mat & x,  // Covariates 的評估點 (grid for X)
                       const arma::vec & h,
                       std::string kernel = "K2_Ep",
                       Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows; // 樣本數
  const arma::uword n_p = X.n_cols; // X 維度
  const arma::uword n_k = x.n_rows; // X grid 數
  const arma::uword n_l = y.n_elem; // Y grid 數 (CDF 評估點數)

  // 取得使用者選擇的 Kernel 函數指標
  KernelFunc k_func = get_kernel_func(kernel);

  // 初始化為 0
  arma::vec Dhat(n_k, arma::fill::zeros);       // 分母
  arma::mat Nhat(n_k, n_l, arma::fill::zeros);  // 分子 (矩陣)
  arma::mat mhat(n_k, n_l, arma::fill::zeros);  // 結果

  // 預先宣告變數
  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);

  // 判斷是否有傳入權重 w
  if (w.isNotNull()) {

    // ============================
    //  情境 A：有權重 (Weighted)
    // ============================

    arma::vec weights = as<arma::vec>(w);

    // 防呆：檢查權重長度
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");

    // 主迴圈
    for (arma::uword k = 0; k < n_k; ++k) {

      xrow_k = x.row(k).t();

      for (arma::uword i = 0; i < n_n; ++i) {
        // 如果權重是 0，直接跳過計算，省時間 (Optional optimization)
        if (weights(i) == 0) continue;

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        // 計算 Product Kernel
        for (arma::uword j = 0; j < n_p; ++j) {

          double u = (Xrow_i(j) - xrow_k(j)) / h(j);
          Kik_h *= k_func(u) / h(j);
        }

        // 乘上權重
        double val = Kik_h * weights(i);

        // 更新分母
        Dhat(k) += val;

        // 更新分子 (針對每個 y grid 計算 Indicator)
        for (arma::uword l = 0; l < n_l; ++l) {

          if (Y(i) <= y(l)) {

            Nhat(k, l) += val;
          }
        }
      }
    }
  } else {

    // ============================
    //  情境 B：無權重 (Unweighted)
    // ============================

    for (arma::uword k = 0; k < n_k; ++k) {

      xrow_k = x.row(k).t();

      for (arma::uword i = 0; i < n_n; ++i) {

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        for (arma::uword j = 0; j < n_p; ++j) {

          double u = (Xrow_i(j) - xrow_k(j)) / h(j);
          Kik_h *= k_func(u) / h(j);
        }

        Dhat(k) += Kik_h;

        for (arma::uword l = 0; l < n_l; ++l) {

          if (Y(i) <= y(l)) {

            Nhat(k, l) += Kik_h;
          }
        }
      }
    }
  }

  // -------------------------------------------------------
  // 計算比率與邊界限制
  // -------------------------------------------------------

  for (arma::uword k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat.row(k) = Nhat.row(k) / Dhat(k);
    }
  }

  mhat.clamp(0.0, 1.0);

  return mhat;
}

// =============================================================
//  整合型 LOOCV 函數
//  Leave-one-out cross-validation for Nadaraya-Watson estimator
// =============================================================

// [[Rcpp::export]]
double CVKNW_rcpp(const arma::vec & Y,
                  const arma::mat & X,
                  const arma::vec & h,
                  std::string kernel = "K2_Ep",
                  Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 檢查是否有權重
  bool has_weights = w.isNotNull();
  arma::vec weights;
  double sum_w = 0.0;

  if (has_weights) {
    weights = as<arma::vec>(w);
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");
    sum_w = arma::sum(weights);
  }

  // =========================================================
  // OpenMP 平行化策略：針對每一個驗證樣本 j 進行平行運算
  // =========================================================

  if (has_weights) {
    // === 情境 A：有權重 (Weighted) ===

    // 使用 reduction 累加 cv_value
    // 使用 dynamic 排程，因為有 weights(j)==0 的跳過機制，負載可能不均
#pragma omp parallel for schedule(dynamic) reduction(+:cv_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      // 優化：若驗證點權重為 0，跳過
      if (weights(j) == 0) continue;

      double Dhat = 0.0;
      double Nhat = 0.0;

      // 內層迴圈：訓練樣本 i
      for (arma::uword i = 0; i < n_n; ++i) {

        // 條件：不是自己 且 權重不為 0
        if (i != j && weights(i) > 0) {

          double Kij_h = 1.0;

          // 計算 Product Kernel (直接讀取 X 元素，不建立暫存向量)
          for (arma::uword l = 0; l < n_p; ++l) {
            // [優化] 直接存取記憶體，大幅減少 Overhead
            double u = (X(i, l) - X(j, l)) / h(l);
            Kij_h *= k_func(u) / h(l);
          }

          double weighted_K = Kij_h * weights(i);
          Dhat += weighted_K;
          Nhat += weighted_K * Y(i);
        }
      }

      double mhat = 0.0;
      if (std::abs(Dhat) > 1e-10) { // 增加數值保護
        mhat = Nhat / Dhat;
      }

      double err = Y(j) - mhat;

      // 累加結果 (OpenMP 會自動處理 reduction)
      cv_value += (err * err) * weights(j);
    }

    if (sum_w != 0) {
      cv_value /= sum_w;
    }

  } else {
    // === 情境 B：無權重 (Unweighted) ===

    // 無權重情境下負載固定，使用 static 排程效率最高
#pragma omp parallel for schedule(static) reduction(+:cv_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      double Dhat = 0.0;
      double Nhat = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        if (i != j) {

          double Kij_h = 1.0;

          for (arma::uword l = 0; l < n_p; ++l) {
            // [優化] 直接讀取 X
            double u = (X(i, l) - X(j, l)) / h(l);
            Kij_h *= k_func(u) / h(l);
          }

          Dhat += Kij_h;
          Nhat += Kij_h * Y(i);
        }
      }

      double mhat = 0.0;
      if (std::abs(Dhat) > 1e-10) {
        mhat = Nhat / Dhat;
      }

      double err = Y(j) - mhat;
      cv_value += err * err;
    }

    cv_value /= n_n;
  }

  return cv_value;
}

// =========================================================
//  整合型 LOOCV-CDF 函數
// =========================================================

// [[Rcpp::export]]
double CVKNWcdf_rcpp(const arma::vec & Y,
                     const arma::mat & X,
                     const arma::vec & h,
                     std::string kernel = "K2_Ep",
                     Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows; // 樣本數
  const arma::uword n_p = X.n_cols; // 變數維度
  double cv_value = 0.0;

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 檢查是否有權重
  bool has_weights = w.isNotNull();
  arma::vec weights;
  double sum_w = 0.0;

  if (has_weights) {
    weights = as<arma::vec>(w);
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");
    sum_w = arma::sum(weights);
  }

  // =========================================================
  // OpenMP 平行化策略：針對 Leave-One-Out 的驗證點 j 進行平行
  // =========================================================

  if (has_weights) {
    // === 情境 A：有權重 (Weighted) ===

    // 使用 dynamic 排程，因為 weights(j)==0 會導致工作量不均
#pragma omp parallel for schedule(dynamic) reduction(+:cv_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      // 優化：權重為 0 直接跳過
      if (weights(j) == 0) continue;

      // [變數私有化] 每個執行緒獨立擁有 Nhat, mhat, Dhat
      arma::vec Nhat(n_n, arma::fill::zeros);
      double Dhat = 0.0;

      // Loop 2: 訓練點 i
      for (arma::uword i = 0; i < n_n; ++i) {

        if (i != j && weights(i) > 0) {

          double Kik_h = 1.0;

          // Kernel 計算 (直接讀取，不建立 Xrow 向量)
          for (arma::uword l = 0; l < n_p; ++l) {
            double u = (X(i, l) - X(j, l)) / h(l);
            Kik_h *= k_func(u) / h(l);
          }

          // 乘上訓練樣本權重
          double val = Kik_h * weights(i);
          Dhat += val;

          // Loop 3: Grid 評估點 k (累積 Indicator)
          // 這是 O(N) 的操作，內嵌在 O(N) 迴圈內，總共 O(N^2)
          for (arma::uword k = 0; k < n_n; ++k) {
            if (Y(i) <= Y(k)) {
              Nhat(k) += val;
            }
          }
        }
      }

      // 計算預測值 mhat 與 CV 誤差
      // 這部分計算量相對小，直接在執行緒內完成
      double local_err_sum = 0.0;

      if (std::abs(Dhat) > 1e-10) {
        // 這裡不需要建立 mhat 向量，直接算直接加總，更省記憶體
        for (arma::uword k = 0; k < n_n; ++k) {
          double m_k = Nhat(k) / Dhat;
          // clamp logic
          if (m_k < 0.0) m_k = 0.0;
          else if (m_k > 1.0) m_k = 1.0;

          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          double err = indicator - m_k;
          local_err_sum += (err * err) * weights(k);
        }
      } else {
        // Dhat 為 0，mhat 為 0
        for (arma::uword k = 0; k < n_n; ++k) {
          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          // err = indicator - 0
          local_err_sum += (indicator * indicator) * weights(k);
        }
      }

      // 加權並累加到總誤差
      cv_value += local_err_sum * weights(j);
    }

    if (sum_w != 0) {
      cv_value /= (sum_w * sum_w);
    }

  } else {
    // === 情境 B：無權重 (Unweighted) ===

    // 使用 static 排程，負載均衡
#pragma omp parallel for schedule(static) reduction(+:cv_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      // [變數私有化] 每個執行緒獨立擁有
      arma::vec Nhat(n_n, arma::fill::zeros);
      double Dhat = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        if (i != j) {

          double Kik_h = 1.0;

          for (arma::uword l = 0; l < n_p; ++l) {
            // 直接讀取，不產生暫存向量
            double u = (X(i, l) - X(j, l)) / h(l);
            Kik_h *= k_func(u) / h(l);
          }

          Dhat += Kik_h;

          for (arma::uword k = 0; k < n_n; ++k) {
            if (Y(i) <= Y(k)) {
              Nhat(k) += Kik_h;
            }
          }
        }
      }

      // 計算誤差
      double local_err_sum = 0.0;

      if (std::abs(Dhat) > 1e-10) {
        for (arma::uword k = 0; k < n_n; ++k) {
          double m_k = Nhat(k) / Dhat;
          if (m_k < 0.0) m_k = 0.0;
          else if (m_k > 1.0) m_k = 1.0;

          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          double err = indicator - m_k;
          local_err_sum += err * err;
        }
      } else {
        for (arma::uword k = 0; k < n_n; ++k) {
          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          local_err_sum += indicator * indicator;
        }
      }

      cv_value += local_err_sum;
    }

    cv_value /= static_cast<double>(n_n * n_n);
  }

  return cv_value;
}

// =========================================================
//  整合型 LSKNW 函數
//  Least squares criterion
// =========================================================

// [[Rcpp::export]]
double LSKNW_rcpp(const arma::vec & Y,
                  const arma::mat & X,
                  const arma::vec & h,
                  std::string kernel = "K2_Ep",
                  Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 檢查是否有權重
  bool has_weights = w.isNotNull();
  arma::vec weights;
  double sum_w = 0.0;

  if (has_weights) {
    weights = as<arma::vec>(w);
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");
    sum_w = arma::sum(weights);
  }

  // =========================================================
  // OpenMP 平行化策略：針對每一個樣本點 j 進行平行運算
  // =========================================================

  if (has_weights) {
    // === 情境 A：有權重 (Weighted) ===

    // 使用 dynamic 排程，因為 weights(j)==0 的跳過機制會導致負載不均
#pragma omp parallel for schedule(dynamic) reduction(+:ls_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      // 優化：權重為 0 直接跳過
      if (weights(j) == 0) continue;

      double Dhat = 0.0;
      double Nhat = 0.0;

      // Loop 2: 訓練樣本 i (包含自己)
      for (arma::uword i = 0; i < n_n; ++i) {

        if (weights(i) > 0) { // 優化：跳過無效訓練點

          double Kik_h = 1.0;

          // Loop 3: Kernel 計算 (直接讀取，不建立暫存向量)
          for (arma::uword l = 0; l < n_p; ++l) {
            // [優化] 直接存取記憶體
            double u = (X(i, l) - X(j, l)) / h(l);
            Kik_h *= k_func(u) / h(l);
          }

          double val = Kik_h * weights(i);
          Dhat += val;
          Nhat += val * Y(i);
        }
      }

      double mhat = 0.0;
      if (std::abs(Dhat) > 1e-10) { // 增加數值保護
        mhat = Nhat / Dhat;
      }

      double err = Y(j) - mhat;
      // 誤差加權累積 (Reduction)
      ls_value += (err * err) * weights(j);
    }

    if (sum_w != 0) ls_value /= sum_w;

  } else {
    // === 情境 B：無權重 (Unweighted) ===

    // 負載均衡，使用 static 排程
#pragma omp parallel for schedule(static) reduction(+:ls_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      double Dhat = 0.0;
      double Nhat = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {
          // [優化] 直接讀取 X
          double u = (X(i, l) - X(j, l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        Dhat += Kik_h;
        Nhat += Kik_h * Y(i);
      }

      double mhat = 0.0;
      if (std::abs(Dhat) > 1e-10) {
        mhat = Nhat / Dhat;
      }

      double err = Y(j) - mhat;
      ls_value += err * err;
    }

    ls_value /= static_cast<double>(n_n);
  }

  return ls_value;
}

// [[Rcpp::export]]
double LSKNW_rcpp0(const arma::vec & Y,
                   const arma::mat & X,
                   const arma::vec & h,
                   std::string kernel = "K2_Ep",
                   Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先宣告變數
  arma::vec Xrow_j(n_p);
  arma::vec Xrow_i(n_p);

  // 判斷是否有傳入權重 w
  if (w.isNotNull()) {

    // ============================
    //  情境 A：有權重 (Weighted)
    // ============================

    arma::vec weights = as<arma::vec>(w);

    // 防呆：檢查權重長度
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");

    double sum_w = arma::sum(weights);

    // Loop 1: 針對每一個樣本點 j (計算誤差)
    for (arma::uword j = 0; j < n_n; ++j) {

      if (weights(j) == 0) continue; // 優化：跳過權重為 0 的點

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      double Nhat = 0.0;

      // Loop 2: 針對每一個樣本點 i (訓練) - 包含自己
      for (arma::uword i = 0; i < n_n; ++i) {

        if (weights(i) == 0) continue; // 優化：跳過權重為 0 的訓練點

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        // Loop 3: Kernel 計算
        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - Xrow_j(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        // 加權累加
        double val = Kik_h * weights(i);
        Dhat += val;
        Nhat += val * Y(i);
      }

      double mhat = 0.0;
      if (Dhat != 0) {

        mhat = Nhat / Dhat;
      }

      double err = Y(j) - mhat;
      // 誤差加權累積
      ls_value += (err * err) * weights(j);
    }

    if (sum_w != 0) ls_value /= sum_w;

  } else {

    // ============================
    //  情境 B：無權重 (Unweighted)
    // ============================

    for (arma::uword j = 0; j < n_n; ++j) {

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      double Nhat = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - Xrow_j(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        Dhat += Kik_h;
        Nhat += Kik_h * Y(i);
      }

      double mhat = 0.0;
      if (Dhat != 0) {

        mhat = Nhat / Dhat;
      }

      double err = Y(j) - mhat;
      ls_value += err * err;
    }

    ls_value /= static_cast<double>(n_n);
  }

  return ls_value;
}

// =========================================================
//  整合型 LSKNW-CDF 函數
// =========================================================

// [[Rcpp::export]]
double LSKNWcdf_rcpp(const arma::vec & Y,
                     const arma::mat & X,
                     const arma::vec & h,
                     std::string kernel = "K2_Ep",
                     Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows; // 樣本數
  const arma::uword n_p = X.n_cols; // 維度
  double ls_value = 0.0;

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 檢查是否有權重
  bool has_weights = w.isNotNull();
  arma::vec weights;
  double sum_w = 0.0;

  if (has_weights) {
    weights = as<arma::vec>(w);
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");
    sum_w = arma::sum(weights);
  }

  // =========================================================
  // OpenMP 平行化策略：針對 Loop j 進行平行
  // =========================================================

  if (has_weights) {
    // === 情境 A：有權重 (Weighted) ===

    // 使用 dynamic 排程，處理 weights(j)==0 跳過的情況
#pragma omp parallel for schedule(dynamic) reduction(+:ls_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      // 優化：權重為 0 直接跳過
      if (weights(j) == 0) continue;

      // [變數私有化] 每個執行緒獨立擁有 Nhat (分子)
      // 不需要 mhat 向量，節省記憶體
      arma::vec Nhat(n_n, arma::fill::zeros);
      double Dhat = 0.0;

      // Loop 2: 訓練點 i
      for (arma::uword i = 0; i < n_n; ++i) {

        if (weights(i) > 0) { // 優化：跳過無效訓練點

          double Kik_h = 1.0;

          // Loop 3: Kernel (直接讀取，不建立 Xrow 向量)
          for (arma::uword l = 0; l < n_p; ++l) {
            double u = (X(i, l) - X(j, l)) / h(l);
            Kik_h *= k_func(u) / h(l);
          }

          // 權重處理
          double val = Kik_h * weights(i);
          Dhat += val;

          // Loop 4: 累積 CDF 分子 (Grid k)
          for (arma::uword k = 0; k < n_n; ++k) {
            if (Y(i) <= Y(k)) {
              Nhat(k) += val;
            }
          }
        }
      }

      // 計算誤差累積 (Loop 5)
      // 直接在迴圈中計算 scalar m_k，不需儲存整條 mhat 向量
      double local_err_sum = 0.0;

      if (std::abs(Dhat) > 1e-10) {
        for (arma::uword k = 0; k < n_n; ++k) {
          double m_k = Nhat(k) / Dhat;

          // Clamp
          if (m_k < 0.0) m_k = 0.0;
          else if (m_k > 1.0) m_k = 1.0;

          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          double err = indicator - m_k;

          // LS 誤差加權：乘上 w(k) 和 w(j)
          local_err_sum += (err * err) * weights(k);
        }
      } else {
        // Dhat 為 0，預測值視為 0
        for (arma::uword k = 0; k < n_n; ++k) {
          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          local_err_sum += (indicator * indicator) * weights(k);
        }
      }

      ls_value += local_err_sum * weights(j);
    }

    // 正規化
    if (sum_w != 0) ls_value /= (sum_w * sum_w);

  } else {
    // === 情境 B：無權重 (Unweighted) ===

    // 使用 static 排程，負載均衡
#pragma omp parallel for schedule(static) reduction(+:ls_value)
    for (arma::uword j = 0; j < n_n; ++j) {

      // [變數私有化]
      arma::vec Nhat(n_n, arma::fill::zeros);
      double Dhat = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {
          // [優化] 直接讀取 X
          double u = (X(i, l) - X(j, l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        Dhat += Kik_h;

        for (arma::uword k = 0; k < n_n; ++k) {
          if (Y(i) <= Y(k)) {
            Nhat(k) += Kik_h;
          }
        }
      }

      double local_err_sum = 0.0;

      if (std::abs(Dhat) > 1e-10) {
        for (arma::uword k = 0; k < n_n; ++k) {
          double m_k = Nhat(k) / Dhat;

          if (m_k < 0.0) m_k = 0.0;
          else if (m_k > 1.0) m_k = 1.0;

          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          double err = indicator - m_k;
          local_err_sum += err * err;
        }
      } else {
        for (arma::uword k = 0; k < n_n; ++k) {
          double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
          local_err_sum += indicator * indicator;
        }
      }

      ls_value += local_err_sum;
    }

    ls_value /= static_cast<double>(n_n * n_n);
  }

  return ls_value;
}

// [[Rcpp::export]]
double LSKNWcdf_rcpp0(const arma::vec & Y,
                      const arma::mat & X,
                      const arma::vec & h,
                      std::string kernel = "K2_Ep",
                      Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先宣告變數
  arma::vec Xrow_j(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  // 判斷是否有傳入權重 w
  if (w.isNotNull()) {

    // ============================
    //  情境 A：有權重 (Weighted)
    // ============================

    arma::vec weights = as<arma::vec>(w);

    // 防呆：檢查權重長度
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");

    double sum_w = arma::sum(weights);

    // Loop 1: 評估點 j
    for (arma::uword j = 0; j < n_n; ++j) {

      if (weights(j) == 0) continue;

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      Nhat.zeros(); // 歸零

      // Loop 2: 訓練點 i (包含自己)
      for (arma::uword i = 0; i < n_n; ++i) {

        if (weights(i) == 0) continue;

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        // Loop 3: Kernel
        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - Xrow_j(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        // 權重處理
        double val = Kik_h * weights(i);
        Dhat += val;

        // Loop 4: 累積 CDF 分子 (Grid k)
        for (arma::uword k = 0; k < n_n; ++k) {

          if (Y(i) <= Y(k)) {
            Nhat(k) += val;
          }
        }
      }

      // 計算預測值 mhat
      if (Dhat != 0) {

        mhat = Nhat / Dhat;
        mhat.clamp(0.0, 1.0); // 限制範圍

      } else {

        mhat.zeros();
      }

      // Loop 5: 計算誤差累積
      for (arma::uword k = 0; k < n_n; ++k) {

        double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
        double err = indicator - mhat(k);
        ls_value += (err * err) * weights(k) * weights(j);
      }
    }

    // 正規化
    if (sum_w != 0) ls_value /= (sum_w * sum_w);
  } else {

    // ============================
    //  情境 B：無權重 (Unweighted)
    // ============================

    for (arma::uword j = 0; j < n_n; ++j) {

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      Nhat.zeros();

      for (arma::uword i = 0; i < n_n; ++i) {

        Xrow_i = X.row(i).t();
        double Kik_h = 1.0;

        for (arma::uword l = 0; l < n_p; ++l) {

          double u = (Xrow_i(l) - Xrow_j(l)) / h(l);
          Kik_h *= k_func(u) / h(l);
        }

        Dhat += Kik_h;

        for (arma::uword k = 0; k < n_n; ++k) {

          if (Y(i) <= Y(k)) {
            Nhat(k) += Kik_h;
          }
        }
      }

      if (Dhat != 0) {

        mhat = Nhat / Dhat;
        mhat.clamp(0.0, 1.0);

      } else {

        mhat.zeros();
      }

      for (arma::uword k = 0; k < n_n; ++k) {

        double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
        double err = indicator - mhat(k);
        ls_value += err * err;
      }
    }

    // 正規化
    ls_value /= static_cast<double>(n_n * n_n);
  }

  return ls_value;
}











