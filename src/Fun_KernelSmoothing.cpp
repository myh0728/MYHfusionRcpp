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
//  定義函數指標與選擇器
// =========================================================

// 定義一個「函數指標型別」，它指向一個「輸入 double&，回傳 double」的函數
typedef double (*KernelFunc)(const double &);

// 輔助函數：根據字串名稱回傳對應的函數指標
KernelFunc get_kernel_func(std::string kernel_name) {

  if (kernel_name == "Epanechnikov" || kernel_name == "K2_Ep") {

    return K2_Ep_inline;

  } else if (kernel_name == "Gaussian" || kernel_name == "K2_G") {

    return K2_G_inline;

  } else if (kernel_name == "K2_Biweight" || kernel_name == "K2_Bw") {

    return K2_Bw_inline;

  } else if (kernel_name == "K4_Biweight" || kernel_name == "K4_Bw") {

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
                   std::string kernel = "Epanechnikov", // 使用者選擇 Kernel
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
                   std::string kernel = "Epanechnikov",
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
                      std::string kernel = "Epanechnikov",
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
                  std::string kernel = "Epanechnikov",
                  Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  // 取得使用者選擇的 Kernel 函數指標
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

    double sum_w = arma::sum(weights); // 計算總權重作為分母

    // 外層迴圈 j: 被留下來驗證的那一個 (Leave-One-Out)
    for (arma::uword j = 0; j < n_n; ++j) {

      // 若驗證點權重為 0，對 CV 無貢獻，直接跳過 (優化)
      if (weights(j) == 0) continue;

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      double Nhat = 0.0;

      // 內層迴圈 i: 用來訓練的樣本
      for (arma::uword i = 0; i < n_n; ++i) {

        // 1. 排除自己 (i != j)
        // 2. 排除權重為 0 的訓練樣本 (優化)
        if (i != j && weights(i) > 0) {

          Xrow_i = X.row(i).t();
          double Kij_h = 1.0;

          // 計算 Product Kernel
          for (arma::uword l = 0; l < n_p; ++l) {

            double u = (Xrow_i(l) - Xrow_j(l)) / h(l);
            Kij_h *= k_func(u) / h(l);
          }

          // 累加時乘上訓練樣本權重 weights(i)
          double weighted_K = Kij_h * weights(i);
          Dhat += weighted_K;
          Nhat += weighted_K * Y(i);
        }
      }

      double mhat = 0.0;

      if (Dhat != 0) {

        mhat = Nhat / Dhat;
      }

      // 計算誤差平方
      double err = Y(j) - mhat;

      // 累加 CV 值時，乘上驗證樣本權重 weights(j)
      cv_value += (err * err) * weights(j);
    }

    // Weighted Mean Squared Error
    if (sum_w != 0) {

      cv_value /= sum_w;
    }

  } else {

    // ============================
    //  情境 B：無權重 (Unweighted)
    // ============================

    for (arma::uword j = 0; j < n_n; ++j) {

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      double Nhat = 0.0;

      for (arma::uword i = 0; i < n_n; ++i) {

        if (i != j) {

          Xrow_i = X.row(i).t();
          double Kij_h = 1.0;

          for (arma::uword l = 0; l < n_p; ++l) {

            double u = (Xrow_i(l) - Xrow_j(l)) / h(l);
            Kij_h *= k_func(u) / h(l);
          }

          Dhat += Kij_h;
          Nhat += Kij_h * Y(i);
        }
      }

      double mhat = 0.0;

      if (Dhat != 0) {

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
                     std::string kernel = "Epanechnikov",
                     Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows; // 樣本數
  const arma::uword n_p = X.n_cols; // 變數維度
  double cv_value = 0.0;

  // 取得使用者選擇的 Kernel 函數指標
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先宣告變數 (避免在迴圈內重複宣告)
  arma::vec Xrow_j(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Nhat(n_n); // 分子
  arma::vec mhat(n_n); // 估計結果

  // 判斷是否有傳入權重 w
  if (w.isNotNull()) {

    // ============================
    //  情境 A：有權重 (Weighted)
    // ============================

    arma::vec weights = as<arma::vec>(w);

    // 防呆：檢查權重長度
    if (weights.n_elem != n_n) stop("Weights length must match X rows.");

    double sum_w = arma::sum(weights); // 計算總權重作為分母

    // Loop 1: 驗證點 j (Leave-One-Out)
    for (arma::uword j = 0; j < n_n; ++j) {

      if (weights(j) == 0) continue; // 優化：權重為 0 直接跳過

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      Nhat.zeros(); // 必須歸零

      // Loop 2: 訓練點 i
      for (arma::uword i = 0; i < n_n; ++i) {

        if (i != j && weights(i) > 0) { // 排除自己 & 排除無效權重

          Xrow_i = X.row(i).t();
          double Kik_h = 1.0;

          // Kernel 計算 (Product Kernel)
          for (arma::uword l = 0; l < n_p; ++l) {

            double u = (Xrow_i(l) - Xrow_j(l)) / h(l);
            Kik_h *= k_func(u) / h(l);
          }

          // 乘上訓練樣本權重
          double val = Kik_h * weights(i);
          Dhat += val;

          // Loop 3: Grid 評估點 k (累積 Indicator)
          // 這裡計算所有 Y(k) 的 CDF 分子
          for (arma::uword k = 0; k < n_n; ++k) {

            if (Y(i) <= Y(k)) {

              Nhat(k) += val;
            }
          }
        }
      }

      // 計算預測值 mhat
      if (Dhat != 0) {

        mhat = Nhat / Dhat;
        mhat.clamp(0.0, 1.0); // 限制在 [0, 1] 之間

      } else {

        mhat.zeros();
      }

      // 計算誤差累積
      for (arma::uword k = 0; k < n_n; ++k) {

        double indicator = (Y(j) <= Y(k)) ? 1.0 : 0.0;
        double err = indicator - mhat(k);
        // CV 誤差加權：乘上 w(k) (grid權重) 和 w(j) (驗證點權重)
        cv_value += (err * err) * weights(k) * weights(j);
      }
    }

    // 正規化
    if (sum_w != 0) {

      cv_value /= (sum_w * sum_w);
    }

  } else {

    // ============================
    //  情境 B：無權重 (Unweighted)
    // ============================

    for (arma::uword j = 0; j < n_n; ++j) {

      Xrow_j = X.row(j).t();
      double Dhat = 0.0;
      Nhat.zeros();

      for (arma::uword i = 0; i < n_n; ++i) {

        if (i != j) {

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
        cv_value += err * err;
      }
    }

    // 正規化
    cv_value /= static_cast<double>(n_n * n_n);
  }

  return cv_value;
}




//==========================================================================

// [[Rcpp::export]]
double LSKNW_K2Ep_rcpp(const arma::vec & Y,
                       const arma::mat & X,
                       const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h;
      Nhat += Kij_h * Y(i);
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    ls_value += pow(Y(j) - mhat, 2);
  }

  ls_value /= n_n;

  return ls_value;
}

// [[Rcpp::export]]
double LSKNW_K2Ep_w_rcpp(const arma::vec & Y,
                         const arma::mat & X,
                         const arma::vec & h,
                         const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h * w(i);
      Nhat += Kij_h * Y(i) * w(i);
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    ls_value += pow(Y(j) - mhat, 2) * w(j);
  }

  ls_value /= sum(w);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNW_K2Bw_rcpp(const arma::vec & Y,
                       const arma::mat & X,
                       const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h;
      Nhat += Kij_h * Y(i);
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    ls_value += pow(Y(j) - mhat, 2);
  }

  ls_value /= n_n;

  return ls_value;
}

// [[Rcpp::export]]
double LSKNW_K2Bw_w_rcpp(const arma::vec & Y,
                         const arma::mat & X,
                         const arma::vec & h,
                         const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h * w(i);
      Nhat += Kij_h * Y(i) * w(i);
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    ls_value += pow(Y(j) - mhat, 2) * w(j);
  }

  ls_value /= sum(w);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNW_K4Bw_rcpp(const arma::vec & Y,
                       const arma::mat & X,
                       const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h;
      Nhat += Kij_h * Y(i);
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    ls_value += pow(Y(j) - mhat, 2);
  }

  ls_value /= n_n;

  return ls_value;
}

// [[Rcpp::export]]
double LSKNW_K4Bw_w_rcpp(const arma::vec & Y,
                         const arma::mat & X,
                         const arma::vec & h,
                         const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h * w(i);
      Nhat += Kij_h * Y(i) * w(i);
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    ls_value += pow(Y(j) - mhat, 2) * w(j);
  }

  ls_value /= sum(w);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNWcdf_K2Ep_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h;

      for (size_t k = 0; k < n_n; ++k) {

        Nhat(k) += Kij_h * (Y(i) <= Y(k));
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      ls_value += pow((Y(j) <= Y(k)) - mhat(k), 2);
    }
  }

  ls_value /= pow(n_n, 2);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNWcdf_K2Ep_w_rcpp(const arma::vec & Y,
                            const arma::mat & X,
                            const arma::vec & h,
                            const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h * w(i);

      for (size_t k = 0; k < n_n; ++k) {

        Nhat(k) += Kij_h * (Y(i) <= Y(k)) * w(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      ls_value += pow((Y(j) <= Y(k)) - mhat(k), 2) * w(k) * w(j);
    }
  }

  ls_value /= pow(sum(w), 2);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNWcdf_K2Bw_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h;

      for (size_t k = 0; k < n_n; ++k) {

        Nhat(k) += Kij_h * (Y(i) <= Y(k));
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      ls_value += pow((Y(j) <= Y(k)) - mhat(k), 2);
    }
  }

  ls_value /= pow(n_n, 2);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNWcdf_K2Bw_w_rcpp(const arma::vec & Y,
                            const arma::mat & X,
                            const arma::vec & h,
                            const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h * w(i);

      for (size_t k = 0; k < n_n; ++k) {

        Nhat(k) += Kij_h * (Y(i) <= Y(k)) * w(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      ls_value += pow((Y(j) <= Y(k)) - mhat(k), 2) * w(k) * w(j);
    }
  }

  ls_value /= pow(sum(w), 2);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNWcdf_K4Bw_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h;

      for (size_t k = 0; k < n_n; ++k) {

        Nhat(k) += Kij_h * (Y(i) <= Y(k));
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      ls_value += pow((Y(j) <= Y(k)) - mhat(k), 2);
    }
  }

  ls_value /= pow(n_n, 2);

  return ls_value;
}

// [[Rcpp::export]]
double LSKNWcdf_K4Bw_w_rcpp(const arma::vec & Y,
                            const arma::mat & X,
                            const arma::vec & h,
                            const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double ls_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dij_h = (Xrow_i - Xrow_j) / h;
      Kij_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
      }

      Dhat += Kij_h * w(i);

      for (size_t k = 0; k < n_n; ++k) {

        Nhat(k) += Kij_h * (Y(i) <= Y(k)) * w(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      ls_value += pow((Y(j) <= Y(k)) - mhat(k), 2) * w(k) * w(j);
    }
  }

  ls_value /= pow(sum(w), 2);

  return ls_value;
}








// For comparison

// [[Rcpp::export]]
double K2_Ep_rcpp_v1(const double & u) {

  double K_val = 0.0;

  K_val = 0.75 * (1.0 - pow(u, 2)) * (abs(u) < 1.0);

  return K_val;
}

// [[Rcpp::export]]
arma::vec KDE_K2Ep_rcpp(const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K2Ep_w_rcpp(const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K2Bw_rcpp(const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K2Bw_w_rcpp(const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K4Bw_rcpp(const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K4Bw_w_rcpp(const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Ep_rcpp(const arma::vec & Y,
                        const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
      Nhat(k) += Kik_h * Y(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Ep_w_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
      Nhat(k) += Kik_h * Y(i) * w(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Bw_rcpp(const arma::vec & Y,
                        const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
      Nhat(k) += Kik_h * Y(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Bw_w_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
      Nhat(k) += Kik_h * Y(i) * w(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K4Bw_rcpp(const arma::vec & Y,
                        const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
      Nhat(k) += Kik_h * Y(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K4Bw_w_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
      Nhat(k) += Kik_h * Y(i) * w(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::mat KNWcdf_K2Ep_rcpp(const arma::vec & Y,
                           const arma::vec & y,
                           const arma::mat & X,
                           const arma::mat & x,
                           const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  const arma::uword n_l = y.n_elem;
  arma::vec Dhat(n_k);
  arma::mat Nhat(n_k, n_l);
  arma::mat mhat(n_k, n_l);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;

      for (size_t l = 0; l < n_l; ++l) {

        Nhat(k, l) += Kik_h * (Y(i) <= y(l));
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      for (size_t l = 0; l < n_l; ++l) {

        mhat(k, l) = Nhat(k, l) / Dhat(k);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
arma::mat KNWcdf_K2Ep_w_rcpp(const arma::vec & Y,
                             const arma::vec & y,
                             const arma::mat & X,
                             const arma::mat & x,
                             const arma::vec & h,
                             const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  const arma::uword n_l = y.n_elem;
  arma::vec Dhat(n_k);
  arma::mat Nhat(n_k, n_l);
  arma::mat mhat(n_k, n_l);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);

      for (size_t l = 0; l < n_l; ++l) {

        Nhat(k, l) += Kik_h * (Y(i) <= y(l)) * w(i);
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      for (size_t l = 0; l < n_l; ++l) {

        mhat(k, l) = Nhat(k, l) / Dhat(k);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
arma::mat KNWcdf_K2Bw_rcpp(const arma::vec & Y,
                           const arma::vec & y,
                           const arma::mat & X,
                           const arma::mat & x,
                           const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  const arma::uword n_l = y.n_elem;
  arma::vec Dhat(n_k);
  arma::mat Nhat(n_k, n_l);
  arma::mat mhat(n_k, n_l);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;

      for (size_t l = 0; l < n_l; ++l) {

        Nhat(k, l) += Kik_h * (Y(i) <= y(l));
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      for (size_t l = 0; l < n_l; ++l) {

        mhat(k, l) = Nhat(k, l) / Dhat(k);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
arma::mat KNWcdf_K2Bw_w_rcpp(const arma::vec & Y,
                             const arma::vec & y,
                             const arma::mat & X,
                             const arma::mat & x,
                             const arma::vec & h,
                             const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  const arma::uword n_l = y.n_elem;
  arma::vec Dhat(n_k);
  arma::mat Nhat(n_k, n_l);
  arma::mat mhat(n_k, n_l);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);

      for (size_t l = 0; l < n_l; ++l) {

        Nhat(k, l) += Kik_h * (Y(i) <= y(l)) * w(i);
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      for (size_t l = 0; l < n_l; ++l) {

        mhat(k, l) = Nhat(k, l) / Dhat(k);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
arma::mat KNWcdf_K4Bw_rcpp(const arma::vec & Y,
                           const arma::vec & y,
                           const arma::mat & X,
                           const arma::mat & x,
                           const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  const arma::uword n_l = y.n_elem;
  arma::vec Dhat(n_k);
  arma::mat Nhat(n_k, n_l);
  arma::mat mhat(n_k, n_l);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;

      for (size_t l = 0; l < n_l; ++l) {

        Nhat(k, l) += Kik_h * (Y(i) <= y(l));
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      for (size_t l = 0; l < n_l; ++l) {

        mhat(k, l) = Nhat(k, l) / Dhat(k);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
arma::mat KNWcdf_K4Bw_w_rcpp(const arma::vec & Y,
                             const arma::vec & y,
                             const arma::mat & X,
                             const arma::mat & x,
                             const arma::vec & h,
                             const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  const arma::uword n_l = y.n_elem;
  arma::vec Dhat(n_k);
  arma::mat Nhat(n_k, n_l);
  arma::mat mhat(n_k, n_l);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);

      for (size_t l = 0; l < n_l; ++l) {

        Nhat(k, l) += Kik_h * (Y(i) <= y(l)) * w(i);
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      for (size_t l = 0; l < n_l; ++l) {

        mhat(k, l) = Nhat(k, l) / Dhat(k);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
double CVKNW_K2Ep_rcpp(const arma::vec & Y,
                       const arma::mat & X,
                       const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h;
        Nhat += Kij_h * Y(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    cv_value += pow(Y(j) - mhat, 2);
  }

  cv_value /= n_n;

  return cv_value;
}

// [[Rcpp::export]]
double CVKNW_K2Ep_w_rcpp(const arma::vec & Y,
                         const arma::mat & X,
                         const arma::vec & h,
                         const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h * w(i);
        Nhat += Kij_h * Y(i) * w(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    cv_value += pow(Y(j) - mhat, 2) * w(j);
  }

  cv_value /= sum(w);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNW_K2Bw_rcpp(const arma::vec & Y,
                       const arma::mat & X,
                       const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h;
        Nhat += Kij_h * Y(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    cv_value += pow(Y(j) - mhat, 2);
  }

  cv_value /= n_n;

  return cv_value;
}

// [[Rcpp::export]]
double CVKNW_K2Bw_w_rcpp(const arma::vec & Y,
                         const arma::mat & X,
                         const arma::vec & h,
                         const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h * w(i);
        Nhat += Kij_h * Y(i) * w(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    cv_value += pow(Y(j) - mhat, 2) * w(j);
  }

  cv_value /= sum(w);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNW_K4Bw_rcpp(const arma::vec & Y,
                       const arma::mat & X,
                       const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h;
        Nhat += Kij_h * Y(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    cv_value += pow(Y(j) - mhat, 2);
  }

  cv_value /= n_n;

  return cv_value;
}

// [[Rcpp::export]]
double CVKNW_K4Bw_w_rcpp(const arma::vec & Y,
                         const arma::mat & X,
                         const arma::vec & h,
                         const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  double Nhat = 0.0;
  double mhat = 0.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat = 0.0;
    mhat = 0.0;

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h * w(i);
        Nhat += Kij_h * Y(i) * w(i);
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
    }

    cv_value += pow(Y(j) - mhat, 2) * w(j);
  }

  cv_value /= sum(w);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNWcdf_K2Ep_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h;

        for (size_t k = 0; k < n_n; ++k) {

          Nhat(k) += Kij_h * (Y(i) <= Y(k));
        }
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      cv_value += pow((Y(j) <= Y(k)) - mhat(k), 2);
    }
  }

  cv_value /= pow(n_n, 2);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNWcdf_K2Ep_w_rcpp(const arma::vec & Y,
                            const arma::mat & X,
                            const arma::vec & h,
                            const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h * w(i);

        for (size_t k = 0; k < n_n; ++k) {

          Nhat(k) += Kij_h * (Y(i) <= Y(k)) * w(i);
        }
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      cv_value += pow((Y(j) <= Y(k)) - mhat(k), 2) * w(k) * w(j);
    }
  }

  cv_value /= pow(sum(w), 2);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNWcdf_K2Bw_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h;

        for (size_t k = 0; k < n_n; ++k) {

          Nhat(k) += Kij_h * (Y(i) <= Y(k));
        }
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      cv_value += pow((Y(j) <= Y(k)) - mhat(k), 2);
    }
  }

  cv_value /= pow(n_n, 2);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNWcdf_K2Bw_w_rcpp(const arma::vec & Y,
                            const arma::mat & X,
                            const arma::vec & h,
                            const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h * w(i);

        for (size_t k = 0; k < n_n; ++k) {

          Nhat(k) += Kij_h * (Y(i) <= Y(k)) * w(i);
        }
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      cv_value += pow((Y(j) <= Y(k)) - mhat(k), 2) * w(k) * w(j);
    }
  }

  cv_value /= pow(sum(w), 2);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNWcdf_K4Bw_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h;

        for (size_t k = 0; k < n_n; ++k) {

          Nhat(k) += Kij_h * (Y(i) <= Y(k));
        }
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      cv_value += pow((Y(j) <= Y(k)) - mhat(k), 2);
    }
  }

  cv_value /= pow(n_n, 2);

  return cv_value;
}

// [[Rcpp::export]]
double CVKNWcdf_K4Bw_w_rcpp(const arma::vec & Y,
                            const arma::mat & X,
                            const arma::vec & h,
                            const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;
  double Dhat = 0.0;
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();
    Dhat = 0.0;
    Nhat.zeros();
    mhat.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K4_Bw_rcpp(Dij_h(l)) / h(l);
        }

        Dhat += Kij_h * w(i);

        for (size_t k = 0; k < n_n; ++k) {

          Nhat(k) += Kij_h * (Y(i) <= Y(k)) * w(i);
        }
      }
    }

    if (Dhat != 0) {

      mhat = Nhat / Dhat;
      mhat.elem(arma::find(mhat < 0)).fill(0);
      mhat.elem(arma::find(mhat > 1)).fill(1);
    }

    for (size_t k = 0; k < n_n; ++k) {

      cv_value += pow((Y(j) <= Y(k)) - mhat(k), 2) * w(k) * w(j);
    }
  }

  cv_value /= pow(sum(w), 2);

  return cv_value;
}




