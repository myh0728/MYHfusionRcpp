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
