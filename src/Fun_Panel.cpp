#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include "head_KernelFunctions.h"

// [[Rcpp::export]]
double criterion_panel_CV_rcpp(const arma::mat & X,
                               const arma::mat & Y,
                               const double & h,
                               std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先計算 h 的倒數，加速運算
  const double inv_h = 1.0 / h;

  double cv_val = 0.0;

  // Loop 1: Target Individual i
  for (arma::uword i = 0; i < n_n; ++i) {

    // Loop 2: Target Time j
    for (arma::uword j = 0; j < n_t; ++j) {

      double target_x = X(i, j);
      double p_ij_sum = 0.0;
      double valid_t_count = 0.0; // [恢復]：用來記錄有效估計的時間點數量

      // Loop 3: Feature Time j_F
      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {

        double sum_num = 0.0;
        double sum_den = 0.0;

        double u_self = (X(i, j_F) - target_x) * inv_h;
        double k_self = k_func(u_self);
        double term_num_self = Y(i, j_F) * k_self;
        double term_den_self = k_self;

        // Loop 4: Training Individuals (純累加)
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {

          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);

          sum_den += k_val;
          sum_num += Y(i_F, j_F) * k_val;
        }

        // Leave-One-Out：扣除自己的貢獻
        sum_den -= term_den_self;
        sum_num -= term_num_self;

        // 若分母足夠大 (可估計)，則累加並計數
        if (std::abs(sum_den) > 1e-10) {
          p_ij_sum += sum_num / sum_den;
          valid_t_count += 1.0; // 計數器 +1
        }
      }

      // [修改]：除以有效次數 (valid_t_count)，而非總時間 (n_t)
      double p_ij = 0.0;
      if (valid_t_count > 0) {
        p_ij = p_ij_sum / valid_t_count;
      }

      double diff_ij = Y(i, j) - p_ij;
      cv_val += diff_ij * diff_ij;
    }
  }

  return cv_val;
}

// [[Rcpp::export]]
double criterion_panel_SS_rcpp(const arma::mat & X,
                               const arma::mat & Y,
                               const double & h,
                               std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先計算 h 的倒數
  const double inv_h = 1.0 / h;

  double ss_val = 0.0;

  // Loop 1: Target Individual i
  for (arma::uword i = 0; i < n_n; ++i) {

    // Loop 2: Target Time j
    for (arma::uword j = 0; j < n_t; ++j) {

      double target_x = X(i, j);
      double p_ij_sum = 0.0;
      double valid_t_count = 0.0;

      // Loop 3: Feature Time j_F
      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {

        double sum_num = 0.0;
        double sum_den = 0.0;

        // Loop 4: Training Individuals (包含自己)
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {

          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);

          sum_den += k_val;
          sum_num += Y(i_F, j_F) * k_val;
        }

        // 若分母足夠大 (可估計)，則累加並計數
        if (std::abs(sum_den) > 1e-10) {
          p_ij_sum += sum_num / sum_den;
          valid_t_count += 1.0;
        }
      }

      // 計算平均值
      double p_ij = 0.0;
      if (valid_t_count > 0) {
        p_ij = p_ij_sum / valid_t_count;
      }

      double diff_ij = Y(i, j) - p_ij;
      ss_val += diff_ij * diff_ij;
    }
  }

  return ss_val;
}

// [[Rcpp::export]]
double criterion_panel_SSeff_rcpp(const arma::mat & X,
                                  const arma::mat & Y,
                                  const arma::cube & V_inv, // 輸入應為變異數矩陣的反矩陣 (n_t, n_t, n_n)
                                  const double & h,
                                  std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先計算 h 的倒數
  const double inv_h = 1.0 / h;

  // 預先宣告殘差向量 (重複使用以節省記憶體配置時間)
  arma::vec resid_i(n_t);

  double ss_val = 0.0;

  // Loop 1: Target Individual i
  for (arma::uword i = 0; i < n_n; ++i) {

    // Loop 2: Target Time j (計算該個體所有時間點的預測值)
    for (arma::uword j = 0; j < n_t; ++j) {

      double target_x = X(i, j);
      double p_ij_sum = 0.0;
      double valid_t_count = 0.0;

      // Loop 3: Feature Time j_F
      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {

        double sum_num = 0.0;
        double sum_den = 0.0;

        // Loop 4: Training Individuals (SS 版本包含自己，無須扣除)
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {

          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);

          sum_den += k_val;
          sum_num += Y(i_F, j_F) * k_val;
        }

        // 避免分母過小
        if (std::abs(sum_den) > 1e-10) {
          p_ij_sum += sum_num / sum_den;
          valid_t_count += 1.0;
        }
      }

      // 計算平均預測值
      double p_ij = 0.0;
      if (valid_t_count > 0) {
        p_ij = p_ij_sum / valid_t_count;
      }

      // 儲存殘差
      resid_i(j) = Y(i, j) - p_ij;
    }

    // 計算加權誤差平方和: (Y_i - p_i)' * V_inv_i * (Y_i - p_i)
    double term_i = arma::as_scalar(resid_i.t() * V_inv.slice(i) * resid_i);

    ss_val += term_i;
  }

  return ss_val;
}

// [[Rcpp::export]]
arma::cube get_effVinv_panel_rcpp(const arma::mat & X,
                                  const arma::mat & Y,
                                  const double & h,
                                  std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows;
  const arma::uword n_t = Y.n_cols;

  KernelFunc k_func = get_kernel_func(kernel);
  const double inv_h = 1.0 / h;

  // 必須儲存 P (用於修正對角線)
  arma::mat P(n_n, n_t, arma::fill::zeros);

  // 結果 Cube
  arma::cube V_inv_cube(n_t, n_t, n_n, arma::fill::zeros);

  double rho_sum = 0.0;
  double valid_c_count = 0.0;

  // =========================================================
  // Loop 1: 計算 P, 預填外積, 累加 Rho
  // =========================================================
  for (arma::uword i = 0; i < n_n; ++i) {

    // 1. 計算 p_ij (Smoothing)
    for (arma::uword j = 0; j < n_t; ++j) {
      double target_x = X(i, j);
      double p_num = 0.0;
      double p_den = 0.0;

      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {
          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);
          p_den += k_val;
          p_num += Y(i_F, j_F) * k_val;
        }
      }
      if (std::abs(p_den) > 1e-10) P(i, j) = p_num / p_den;
    }

    // 2. [關鍵] 直接使用外積填滿，這是最快的方法
    // 對角線此時是 p^2 (錯誤值)，但我們稍後會覆蓋它，所以沒關係
    arma::vec p_i = P.row(i).t();
    V_inv_cube.slice(i) = p_i * p_i.t();

    // 3. 累加 Rho
    // 這裡我們可以直接讀取 cube 裡的數值 (即 p_j * p_k)
    for (arma::uword j = 0; j < n_t; ++j) {
      double target_x_j = X(i, j);

      for (arma::uword k = 0; k < n_t; ++k) {

        // [優化] 通常 Rho 只考慮非對角線 (Correlation)
        if (j == k) continue;

        double p_prod = V_inv_cube(j, k, i); // 直接讀取，省去乘法

        // 若 p_j * p_k 太小，計算不穩定，跳過
        if (std::abs(p_prod) <= 1e-10) continue;

        double target_x_k = X(i, k);
        double c_Ntor = 0.0;
        double c_Dtor = 0.0;

        for (arma::uword i_c = 0; i_c < n_n; ++i_c) {
          double u_j = (X(i_c, j) - target_x_j) * inv_h;
          double u_k = (X(i_c, k) - target_x_k) * inv_h;
          double k_val = k_func(u_j) * k_func(u_k);

          c_Dtor += k_val;
          c_Ntor += Y(i_c, j) * Y(i_c, k) * k_val;
        }

        if (std::abs(c_Dtor) > 1e-10) {
          double expected_joint = c_Ntor / c_Dtor;
          rho_sum += (expected_joint - p_prod) / p_prod;
          valid_c_count += 1.0;
        }
      }
    }
  }

  double rho = (valid_c_count > 0) ? (rho_sum / valid_c_count) : 0.0;

  // =========================================================
  // Loop 2: 縮放與對角線覆蓋 (Overwrite Strategy)
  // =========================================================
  for (arma::uword i = 0; i < n_n; ++i) {

    arma::mat & V_i = V_inv_cube.slice(i);

    // 1. 全局縮放: 所有元素 (含對角線) 都乘上 rho
    // 這是 BLAS level 運算，非常快
    V_i *= rho;

    // 2. 對角線覆蓋: 修正為正確的 p(1-p)
    // 這裡不需要除以 rho，因為我們是直接 assign 新值
    arma::vec p_i = P.row(i).t();

    // 為了向量化計算，先算出 Variance 向量
    arma::vec var_vec = p_i % (1.0 - p_i);

    // 加強數值穩定性 (Nugget effect)
    // 這種 element-wise 操作在 Armadillo 中也是優化的
    var_vec.elem( find(var_vec < 1e-6) ).fill(1e-6);
    var_vec += 1e-8;

    // 覆蓋對角線
    V_i.diag() = var_vec;

    // 3. 計算反矩陣 (Cholesky 加速)
    arma::mat V_inv_temp;
    bool success = arma::inv_sympd(V_inv_temp, V_i);

    if (success) {
      V_i = V_inv_temp;
    } else {
      V_i = arma::inv(V_i);
    }
  }

  return V_inv_cube;
}


