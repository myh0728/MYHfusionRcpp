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

  // 預先計算 h 的倒數
  const double inv_h = 1.0 / h;

  double cv_val = 0.0;

  // Loop 1: Target Individual i (平行化)
  // 使用 reduction 安全累加 cv_val
  // 使用 static 排程，負載均衡
#pragma omp parallel for schedule(static) reduction(+:cv_val)
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

        // 計算自己的貢獻 (稍後要扣除)
        // 這些變數宣告在迴圈內，是執行緒私有的，安全
        double u_self = (X(i, j_F) - target_x) * inv_h;
        double k_self = k_func(u_self);
        double term_num_self = Y(i, j_F) * k_self;
        double term_den_self = k_self;

        // Loop 4: Training Individuals (全體累加)
        // 內層維持單線程，Column-Major 讀取 X(i_F, j_F) 效率最高
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {
          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);

          sum_den += k_val;
          sum_num += Y(i_F, j_F) * k_val;
        }

        // Leave-One-Out：扣除自己的貢獻
        // 這種 "先加總再扣除" 的寫法比 "if (i_F != i)" 更快，因為減少了分支判斷
        sum_den -= term_den_self;
        sum_num -= term_num_self;

        // 若分母足夠大 (可估計)，則累加並計數
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

      double diff_ij = Y(i, j) - p_ij;

      // 累加到 reduction 變數
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

  // Loop 1: Target Individual i (平行化)
  // 使用 reduction(+:ss_val) 確保所有執行緒算出的誤差平方能正確加總
  // 使用 schedule(static) 因為每個 i 的計算量完全相同
#pragma omp parallel for schedule(static) reduction(+:ss_val)
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

        // Loop 4: Training Individuals
        // 內層維持單線程，利用 Column-Major 優勢快速讀取 X, Y
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {

          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);

          sum_den += k_val;
          sum_num += Y(i_F, j_F) * k_val;
        }

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

      double diff_ij = Y(i, j) - p_ij;

      // 直接累加到 reduction 變數
      ss_val += diff_ij * diff_ij;
    }
  }

  return ss_val;
}

// [[Rcpp::export]]
double criterion_panel_SSeff_rcpp(const arma::mat & X,
                                  const arma::mat & Y,
                                  const arma::cube & V_inv, // 輸入為 Cube (n_t, n_t, n_n)
                                  const double & h,
                                  std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先計算 h 的倒數
  const double inv_h = 1.0 / h;

  double ss_val = 0.0;

  // Loop 1: Target Individual i (平行化)
  // 使用 reduction(+:ss_val) 確保加總正確
  // 使用 schedule(static) 因為每個 i 的負載是均衡的
#pragma omp parallel for schedule(static) reduction(+:ss_val)
  for (arma::uword i = 0; i < n_n; ++i) {

    // [關鍵]：將殘差向量移入迴圈，確保每個執行緒獨立擁有
    arma::vec resid_i(n_t);

    // Loop 2: Target Time j
    for (arma::uword j = 0; j < n_t; ++j) {

      double target_x = X(i, j);
      double p_ij_sum = 0.0;
      double valid_t_count = 0.0;

      // Loop 3: Feature Time j_F
      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {

        double sum_num = 0.0;
        double sum_den = 0.0;

        // Loop 4: Training Individuals
        // 內層維持單線程，利用 Column-Major 優勢快速讀取 X, Y
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {

          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);

          sum_den += k_val;
          sum_num += Y(i_F, j_F) * k_val;
        }

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
    // V_inv 是 const reference，讀取 slice(i) 是 Thread-safe 的
    double term_i = arma::as_scalar(resid_i.t() * V_inv.slice(i) * resid_i);

    ss_val += term_i;
  }

  return ss_val;
}

// [[Rcpp::export]]
double criterion_panel_SSw_rcpp(const arma::mat & X,
                                const arma::mat & Y,
                                const arma::mat & W_inv, // 共同變異數矩陣的反矩陣 (n_t, n_t)
                                const double & h,
                                std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先計算 h 的倒數
  const double inv_h = 1.0 / h;

  double ss_val = 0.0;

  // Loop 1: Target Individual i (平行化)
  // 使用 reduction 處理 ss_val 的累加
  // 使用 static 排程，因為每個 i 的負載是均衡的
#pragma omp parallel for schedule(static) reduction(+:ss_val)
  for (arma::uword i = 0; i < n_n; ++i) {

    // [重要修正]：必須在平行迴圈內宣告，確保每個執行緒有獨立的殘差向量
    arma::vec resid_i(n_t);

    // Loop 2: Target Time j
    for (arma::uword j = 0; j < n_t; ++j) {

      double target_x = X(i, j);
      double p_ij_sum = 0.0;
      double valid_t_count = 0.0;

      // Loop 3: Feature Time j_F
      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {

        double sum_num = 0.0;
        double sum_den = 0.0;

        // Loop 4: Training Individuals
        // 內層維持單執行緒掃描，因為 X 是 Column-Major，這樣讀取記憶體效率高
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {

          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);

          sum_den += k_val;
          sum_num += Y(i_F, j_F) * k_val;
        }

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

      // 儲存殘差 (寫入各執行緒私有的 resid_i)
      resid_i(j) = Y(i, j) - p_ij;
    }

    // 計算加權誤差平方和: (Y_i - p_i)' * W_inv * (Y_i - p_i)
    // W_inv 是 const reference，多執行緒同時讀取是安全的
    double term_i = arma::as_scalar(resid_i.t() * W_inv * resid_i);

    ss_val += term_i;
  }

  return ss_val;
}

// [[Rcpp::export]]
arma::cube get_effVinv_panel_rcpp(const arma::mat & X,
                                  const arma::mat & Y,
                                  const double & h,
                                  std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // N
  const arma::uword n_t = Y.n_cols; // T

  KernelFunc k_func = get_kernel_func(kernel);
  const double inv_h = 1.0 / h;

  // 必須儲存 P (用於修正對角線)
  arma::mat P(n_n, n_t, arma::fill::zeros);

  // 結果 Cube
  arma::cube V_inv_cube(n_t, n_t, n_n, arma::fill::zeros);

  double rho_sum = 0.0;
  double valid_c_count = 0.0;

  // =========================================================
  // Loop 1: 計算 P, 預填外積, 累加 Rho (平行化)
  // =========================================================
  // 使用 dynamic 排程，因為內部有 continue 跳過計算，負載不均
#pragma omp parallel for schedule(dynamic) reduction(+:rho_sum, valid_c_count)
  for (arma::uword i = 0; i < n_n; ++i) {

    // 1. 計算 p_ij (Smoothing)
    for (arma::uword j = 0; j < n_t; ++j) {
      double target_x = X(i, j);
      double p_num = 0.0;
      double p_den = 0.0;

      // 內層這兩圈是讀取唯讀的 X, Y，無須擔心
      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {
          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);
          p_den += k_val;
          p_num += Y(i_F, j_F) * k_val;
        }
      }
      // P(i, j) 寫入各自獨立的記憶體位置，安全
      if (std::abs(p_den) > 1e-10) P(i, j) = p_num / p_den;
    }

    // 2. 預填外積
    arma::vec p_i = P.row(i).t();
    V_inv_cube.slice(i) = p_i * p_i.t(); // 寫入獨立 slice，安全

    // 3. 累加 Rho
    for (arma::uword j = 0; j < n_t; ++j) {
      double target_x_j = X(i, j);

      for (arma::uword k = 0; k < n_t; ++k) {
        if (j == k) continue;

        double p_prod = V_inv_cube(j, k, i);

        // 這個 check 會導致不同 i 的運算時間差異很大 -> 所以用 dynamic schedule
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
          // reduction 變數累加，OpenMP 會自動處理
          rho_sum += (expected_joint - p_prod) / p_prod;
          valid_c_count += 1.0;
        }
      }
    }
  }

  // 計算平均 rho (必須在 Loop 1 結束，同步後才能算)
  double rho = (valid_c_count > 0) ? (rho_sum / valid_c_count) : 0.0;

  // =========================================================
  // Loop 2: 縮放與對角線覆蓋 (平行化)
  // =========================================================
  // 矩陣反運算負載相當固定，使用 static 排程即可
#pragma omp parallel for schedule(static)
  for (arma::uword i = 0; i < n_n; ++i) {

    arma::mat & V_i = V_inv_cube.slice(i); // Reference to slice

    // 1. 全局縮放
    V_i *= rho;

    // 2. 對角線覆蓋
    arma::vec p_i = P.row(i).t();
    arma::vec var_vec = p_i % (1.0 - p_i);
    var_vec.elem( find(var_vec < 1e-6) ).fill(1e-6);
    var_vec += 1e-8;
    V_i.diag() = var_vec;

    // 3. 計算反矩陣
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

// [[Rcpp::export]]
arma::vec get_linkF_panel_rcpp(const arma::mat & X,
                               const arma::mat & Y,
                               const arma::vec & u_points, // 評估點向量
                               const double & h,
                               std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)
  const arma::uword n_l = u_points.n_elem; // 評估點數 (L)

  // 取得 Kernel 函數
  KernelFunc k_func = get_kernel_func(kernel);

  // 預先計算 h 的倒數
  const double inv_h = 1.0 / h;

  // 初始化結果向量 (填滿 0)
  arma::vec Fhat(n_l, arma::fill::zeros);

  // Loop 1: 針對每一個評估點 u (平行化)
  // 因為每個 u 的計算互不影響，這裡是平行運算的最佳切入點
#pragma omp parallel for schedule(static)
  for (arma::uword l = 0; l < n_l; ++l) {

    double target_u = u_points(l);
    double p_ij_sum = 0.0;
    double valid_t_count = 0.0;

    // Loop 2: 遍歷每一個時間點 (Time-wise Smoothing)
    for (arma::uword j_F = 0; j_F < n_t; ++j_F) {

      double sum_num = 0.0;
      double sum_den = 0.0;

      // Loop 3: 遍歷每一個個體
      for (arma::uword i_F = 0; i_F < n_n; ++i_F) {

        double u = (X(i_F, j_F) - target_u) * inv_h;
        double k_val = k_func(u);

        sum_den += k_val;
        sum_num += Y(i_F, j_F) * k_val;
      }

      // 避免分母過小 (數值穩定性)
      if (std::abs(sum_den) > 1e-10) {
        p_ij_sum += sum_num / sum_den;
        valid_t_count += 1.0;
      }
    }

    // 計算平均值 (Average of Marginal Smoothers)
    if (valid_t_count > 0) {
      Fhat(l) = p_ij_sum / valid_t_count;
    }
  }

  return Fhat;
}




