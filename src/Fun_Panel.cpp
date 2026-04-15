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

      p_ij = std::max(0.0, std::min(1.0, p_ij));
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

      p_ij = std::max(0.0, std::min(1.0, p_ij));
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

      p_ij = std::max(0.0, std::min(1.0, p_ij));
      // 儲存殘差
      resid_i(j) = Y(i, j) - p_ij;
    }

    // 計算加權誤差平方和: (Y_i - p_i)' * V_inv_i * (Y_i - p_i)
    // V_inv 是 const reference，讀取 slice(i) 是 Thread-safe 的
    ss_val += arma::dot(resid_i, V_inv.slice(i) * resid_i);
  }

  return ss_val;
}

// [[Rcpp::export]]
arma::cube get_Vinv_panel_frailty_rcpp(const arma::mat & X,
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

      if (std::abs(p_den) > 1e-10) {
        double raw_p = p_num / p_den;
        P(i, j) = std::max(0.0, std::min(1.0, raw_p));
      } else {
        // 若分母過小，給予預設值 0 (P 初始化已為 0，這行只是明確寫出邏輯)
        P(i, j) = 0.0;
      }
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
      double raw_f = p_ij_sum / valid_t_count;
      Fhat(l) = std::max(0.0, std::min(1.0, raw_f));
    }
  }

  return Fhat;
}

// [[Rcpp::export]]
arma::cube get_Vinv_panel_exchangable_rcpp(const arma::mat & X,
                                           const arma::mat & Y,
                                           const double & h,
                                           std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // N
  const arma::uword n_t = Y.n_cols; // T

  KernelFunc k_func = get_kernel_func(kernel);
  const double inv_h = 1.0 / h;

  // 用來儲存標準差 (Sigma)，供第二階段使用
  arma::mat Sigma_mat(n_n, n_t, arma::fill::zeros);

  // 結果 Cube
  arma::cube V_inv_cube(n_t, n_t, n_n, arma::fill::zeros);

  double rho_sum = 0.0;

  // =========================================================
  // Stage 1: 平行計算 P (Smoothing), Sigma, 與 Rho 的分子
  // =========================================================
#pragma omp parallel for schedule(dynamic) reduction(+:rho_sum)
  for (arma::uword i = 0; i < n_n; ++i) {

    arma::vec p_i(n_t);
    arma::vec sigma_i(n_t);

    // 1. 計算 P_ij (Smoothing)
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

      if (std::abs(p_den) > 1e-10) {
        double raw_p = p_num / p_den;
        p_i(j) = std::max(0.0, std::min(1.0, raw_p));
      } else {
        p_i(j) = 0.0;
      }
    }

    // 2. 計算標準差 Sigma_i
    arma::vec var_i = p_i % (1.0 - p_i);

    // [重要] 數值穩定性處理 (Nugget effect)
    // 確保後面除以 sigma 時不會爆炸
    for(arma::uword k = 0; k < n_t; ++k) {
      if(var_i(k) < 1e-8) var_i(k) = 1e-8;
    }
    sigma_i = arma::sqrt(var_i);

    // 存入矩陣
    Sigma_mat.row(i) = sigma_i.t();

    // 3. 累加 Rho 分子
    for (arma::uword j = 0; j < n_t - 1; ++j) {
      double resid_j_std = (Y(i, j) - p_i(j)) / sigma_i(j);

      for (arma::uword k = j + 1; k < n_t; ++k) {
        double resid_k_std = (Y(i, k) - p_i(k)) / sigma_i(k);
        rho_sum += resid_j_std * resid_k_std;
      }
    }
  }

  // 計算平均 Rho
  double rho = rho_sum * 2.0 / (double)(n_n * n_t * (n_t - 1));

  // 防呆：限制 rho 的範圍，避免分母為 0
  // 理論極限是 -1/(T-1) < rho < 1
  double lower_bound = -1.0 / (double)(n_t - 1) + 1e-6;
  if (rho > 0.999) rho = 0.999;
  if (rho < lower_bound) rho = lower_bound;

  // =========================================================
  // Stage 2: 使用解析解 (Analytical Solution) 直接填入 V_inv
  // =========================================================

  // 預先計算 R inverse 的係數 (A 和 B)
  // 這些係數對所有個體都是一樣的
  double denom = (1.0 - rho) * (1.0 + (n_t - 1) * rho);
  double r_inv_diag = (1.0 + (n_t - 2) * rho) / denom; // 對角線係數 A
  double r_inv_off  = -rho / denom;                    // 非對角線係數 B

#pragma omp parallel for schedule(static)
  for (arma::uword i = 0; i < n_n; ++i) {

    // 取出該個體的 sigma
    arma::vec sigma_i = Sigma_mat.row(i).t();
    arma::mat & V_inv = V_inv_cube.slice(i); // 引用 Cube 的 slice

    // 雙層迴圈填值
    for (arma::uword j = 0; j < n_t; ++j) {

      // 1. 對角線元素: A / sigma_j^2
      // 注意：sigma_i(j) * sigma_i(j) 就是變異數
      V_inv(j, j) = r_inv_diag / (sigma_i(j) * sigma_i(j));

      // 2. 非對角線元素: B / (sigma_j * sigma_k)
      // 利用對稱性一次填兩個
      for (arma::uword k = j + 1; k < n_t; ++k) {
        double val = r_inv_off / (sigma_i(j) * sigma_i(k));
        V_inv(j, k) = val;
        V_inv(k, j) = val;
      }
    }
  }

  return V_inv_cube;
}

// [[Rcpp::export]]
arma::cube get_Vinv_panel_AR1_rcpp(const arma::mat & X,
                                   const arma::mat & Y,
                                   const double & h,
                                   std::string kernel = "K2_Ep") {

  const arma::uword n_n = Y.n_rows; // N
  const arma::uword n_t = Y.n_cols; // T

  KernelFunc k_func = get_kernel_func(kernel);
  const double inv_h = 1.0 / h;

  // 儲存 Sigma 供 Stage 2 使用
  arma::mat Sigma_mat(n_n, n_t, arma::fill::zeros);

  // 結果 Cube
  arma::cube V_inv_cube(n_t, n_t, n_n, arma::fill::zeros);

  double rho_sum = 0.0;

  // =========================================================
  // Stage 1: 平行計算 P (Smoothing), Sigma, 與 Rho
  // =========================================================
#pragma omp parallel for schedule(dynamic) reduction(+:rho_sum)
  for (arma::uword i = 0; i < n_n; ++i) {

    arma::vec p_i(n_t);
    arma::vec sigma_i(n_t);

    // 1. Smoothing (Kernel Regression)
    for (arma::uword j = 0; j < n_t; ++j) {
      double target_x = X(i, j);
      double p_num = 0.0;
      double p_den = 0.0;

      for (arma::uword j_F = 0; j_F < n_t; ++j_F) {
        for (arma::uword i_F = 0; i_F < n_n; ++i_F) {
          // 優化：直接讀取，不產生暫存向量
          double u = (X(i_F, j_F) - target_x) * inv_h;
          double k_val = k_func(u);
          p_den += k_val;
          p_num += Y(i_F, j_F) * k_val;
        }
      }

      if (std::abs(p_den) > 1e-10) {
        double raw_p = p_num / p_den;
        p_i(j) = std::max(0.0, std::min(1.0, raw_p));
      } else {
        p_i(j) = 0.0;
      }
    }

    // 2. 計算標準差 Sigma_i
    arma::vec var_i = p_i % (1.0 - p_i);

    // 數值穩定性處理
    for(arma::uword k = 0; k < n_t; ++k) {
      if(var_i(k) < 1e-8) var_i(k) = 1e-8;
    }
    sigma_i = arma::sqrt(var_i);
    Sigma_mat.row(i) = sigma_i.t();

    // 3. 累加 Rho 分子 (AR1 只需要看相鄰兩項)
    for (arma::uword j = 0; j < n_t - 1; ++j) {
      double resid_j_std = (Y(i, j) - p_i(j)) / sigma_i(j);
      double resid_jj_std = (Y(i, j + 1) - p_i(j + 1)) / sigma_i(j + 1);

      rho_sum += resid_j_std * resid_jj_std;
    }
  }

  // 計算平均 Rho
  double rho = rho_sum / (double)(n_n * (n_t - 1));

  // [重要] 防呆與邊界限制
  // AR(1) 解析解的分母是 1 - rho^2，若 rho 接近 1 或 -1 會爆炸
  if (rho > 0.999) rho = 0.999;
  if (rho < -0.999) rho = -0.999;

  // =========================================================
  // Stage 2: 使用 AR(1) 解析反矩陣公式填入 V_inv
  // =========================================================

  // 預先計算 AR(1) inverse 的常數係數
  double rho_sq = rho * rho;
  double scale = 1.0 / (1.0 - rho_sq); // 1 / (1 - rho^2)

  // 係數 A: 用於頭尾對角線 (1 * scale)
  // 係數 B: 用於中間對角線 ((1 + rho^2) * scale)
  // 係數 C: 用於次對角線   (-rho * scale)

  double coef_mid_diag = (1.0 + rho_sq) * scale;
  double coef_off_diag = -rho * scale;

#pragma omp parallel for schedule(static)
  for (arma::uword i = 0; i < n_n; ++i) {

    arma::vec sigma_i = Sigma_mat.row(i).t();
    arma::mat & V_inv = V_inv_cube.slice(i); // 引用，直接修改記憶體

    // AR(1) 反矩陣是三對角矩陣 (Tridiagonal)，其餘為 0
    // 我們只遍歷需要填值的地方，大幅提升速度

    // 1. 填入對角線 (Diagonal)
    for (arma::uword j = 0; j < n_t; ++j) {
      double sig_sq = sigma_i(j) * sigma_i(j);

      if (j == 0 || j == n_t - 1) {
        // 頭尾: 1 / (sig^2 * (1-rho^2))
        V_inv(j, j) = scale / sig_sq;
      } else {
        // 中間: (1+rho^2) / (sig^2 * (1-rho^2))
        V_inv(j, j) = coef_mid_diag / sig_sq;
      }
    }

    // 2. 填入次對角線 (Off-diagonal, |j-k|=1)
    for (arma::uword j = 0; j < n_t - 1; ++j) {
      // 公式: -rho / (sig_j * sig_k * (1-rho^2))
      double val = coef_off_diag / (sigma_i(j) * sigma_i(j + 1));

      V_inv(j, j + 1) = val;
      V_inv(j + 1, j) = val; // 對稱
    }

    // 其餘位置保持為 0 (初始化時已經填 0 了)
  }

  return V_inv_cube;
}

// [[Rcpp::export]]
double criterion_panel_logit_SS_rcpp(const arma::mat & X,
                                     const arma::mat & Y,
                                     const double & intercept) {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)

  double ss_val = 0.0;

  // 使用 OpenMP 平行化，負載均衡使用 static 排程
#pragma omp parallel for schedule(static) reduction(+:ss_val)
  for (arma::uword i = 0; i < n_n; ++i) {

    for (arma::uword j = 0; j < n_t; ++j) {

      double SI_ij = intercept + X(i, j);
      double p_ij = 1.0 / (1.0 + std::exp(-SI_ij));
      double diff = Y(i, j) - p_ij;
      ss_val += diff * diff;
    }
  }

  return ss_val;
}

// [[Rcpp::export]]
double criterion_panel_logit_SSeff_rcpp(const arma::mat & X,
                                        const arma::mat & Y,
                                        const double & intercept,
                                        const arma::cube & V_inv) {

  const arma::uword n_n = Y.n_rows; // 個體數 (N)
  const arma::uword n_t = Y.n_cols; // 時間點數 (T)

  double ss_val = 0.0;

  // 使用 OpenMP 平行化，靜態排程
#pragma omp parallel for schedule(static) reduction(+:ss_val)
  for (arma::uword i = 0; i < n_n; ++i) {

    // 宣告在迴圈內，每個執行緒擁有自己獨立的 diff 向量 (長度為 T)
    arma::vec diff(n_t);

    for (arma::uword j = 0; j < n_t; ++j) {

      double SI_ij = intercept + X(i, j);
      double p_ij = 1.0 / (1.0 + std::exp(-SI_ij));

      // 直接算出殘差並存入，省去 vector 相減和轉置 (.t()) 的開銷
      diff(j) = Y(i, j) - p_ij;
    }

    // [效能關鍵]
    // 使用 arma::dot(x, A * x) 來計算二次型 (Quadratic Form: x^T * A * x)
    // 這是 Armadillo 中計算此類矩陣乘法最快且最標準的寫法，直接回傳 double
    ss_val += arma::dot(diff, V_inv.slice(i) * diff);
  }

  return ss_val;
}

// [[Rcpp::export]]
arma::cube get_Vinv_panel_logit_exchangable_rcpp(const arma::mat & X,
                                                 const arma::mat & Y,
                                                 const double & intercept) {

  const arma::uword n_n = Y.n_rows; // N
  const arma::uword n_t = Y.n_cols; // T

  // 儲存標準差
  arma::mat Sigma_mat(n_n, n_t, arma::fill::zeros);
  arma::cube V_inv_cube(n_t, n_t, n_n, arma::fill::zeros);

  double rho_sum = 0.0;

  // =========================================================
  // Stage 1: 平行計算 P, Sigma, 與 Rho 的分子
  // =========================================================
  // 工作量完全對等，使用 static 排程效能最好
#pragma omp parallel for schedule(static) reduction(+:rho_sum)
  for (arma::uword i = 0; i < n_n; ++i) {

    // 1. 高度向量化：一次算出該個體所有時間點的機率 P
    arma::vec SI_i = intercept + X.row(i).t();
    arma::vec p_i = 1.0 / (1.0 + arma::exp(-SI_i));

    // 2. 計算標準差 Sigma_i 與下限保護
    arma::vec var_i = p_i % (1.0 - p_i);
    var_i.elem( arma::find(var_i < 1e-8) ).fill(1e-8);
    arma::vec sigma_i = arma::sqrt(var_i);

    Sigma_mat.row(i) = sigma_i.t();

    // 3. 累加 Rho 分子
    // 先算出標準化殘差向量 Z
    arma::vec resid_std = (Y.row(i).t() - p_i) / sigma_i;

    // (總和)^2
    double sum_resid = arma::sum(resid_std);
    // 平方和
    double sum_sq_resid = arma::accu(arma::square(resid_std));

    // 套用公式: sum(Z_j * Z_k) = 0.5 * [ (sum Z)^2 - sum(Z^2) ]
    rho_sum += 0.5 * (sum_resid * sum_resid - sum_sq_resid);
  }

  // 計算平均 Rho
  double rho = rho_sum * 2.0 / (double)(n_n * n_t * (n_t - 1));

  // 防呆：限制 rho 的範圍
  double lower_bound = -1.0 / (double)(n_t - 1) + 1e-6;
  if (rho > 0.999) rho = 0.999;
  if (rho < lower_bound) rho = lower_bound;

  // =========================================================
  // Stage 2: 使用解析解 (Analytical Solution) 直接填入 V_inv
  // =========================================================

  double denom = (1.0 - rho) * (1.0 + (n_t - 1) * rho);
  double r_inv_diag = (1.0 + (n_t - 2) * rho) / denom; // 係數 A
  double r_inv_off  = -rho / denom;                    // 係數 B

#pragma omp parallel for schedule(static)
  for (arma::uword i = 0; i < n_n; ++i) {

    // 取出 sigma_i 的倒數
    arma::vec inv_sig = 1.0 / Sigma_mat.row(i).t();

    arma::mat & V_inv = V_inv_cube.slice(i);

    // 高度向量化：免除所有迴圈，直接使用矩陣外積生成 V_inv
    // 1. 先把整個矩陣當成非對角線元素計算: B * (inv_sig * inv_sig^T)
    V_inv = r_inv_off * (inv_sig * inv_sig.t());

    // 2. 將對角線元素覆蓋為正確的對角線值: A * (inv_sig^2)
    V_inv.diag() = r_inv_diag * (inv_sig % inv_sig);
  }

  return V_inv_cube;
}

// [[Rcpp::export]]
arma::cube get_Vinv_panel_logit_AR1_rcpp(const arma::mat & X,
                                         const arma::mat & Y,
                                         const double & intercept) {

  const arma::uword n_n = Y.n_rows; // N
  const arma::uword n_t = Y.n_cols; // T

  // 儲存 Sigma 供 Stage 2 使用
  arma::mat Sigma_mat(n_n, n_t, arma::fill::zeros);

  // 結果 Cube
  arma::cube V_inv_cube(n_t, n_t, n_n, arma::fill::zeros);

  double rho_sum = 0.0;

  // =========================================================
  // Stage 1: 高度向量化計算 P, Sigma, 與 Rho
  // =========================================================
  // 工作量完全對等，使用 static 排程效能最好
#pragma omp parallel for schedule(static) reduction(+:rho_sum)
  for (arma::uword i = 0; i < n_n; ++i) {

    // 1. 高度向量化：一次算出該個體所有時間點的機率 P
    arma::vec SI_i = intercept + X.row(i).t();
    arma::vec p_i = 1.0 / (1.0 + arma::exp(-SI_i));

    // 2. 計算標準差 Sigma_i 與下限保護
    arma::vec var_i = p_i % (1.0 - p_i);
    var_i.elem( arma::find(var_i < 1e-8) ).fill(1e-8);
    arma::vec sigma_i = arma::sqrt(var_i);
    Sigma_mat.row(i) = sigma_i.t();

    // 3. 累加 Rho 分子 (免迴圈，使用錯位內積)
    if (n_t > 1) {
      // 先算出整條標準化殘差向量
      arma::vec resid_std = (Y.row(i).t() - p_i) / sigma_i;

      // AR(1) 分子: sum( Z_j * Z_{j+1} )
      // 用 subvec 取出 0~T-2 與 1~T-1，然後做內積
      rho_sum += arma::dot(resid_std.subvec(0, n_t - 2),
                           resid_std.subvec(1, n_t - 1));
    }
  }

  // 計算平均 Rho
  double rho = rho_sum / (double)(n_n * (n_t - 1));

  // [重要] 防呆與邊界限制
  if (rho > 0.999) rho = 0.999;
  if (rho < -0.999) rho = -0.999;

  // =========================================================
  // Stage 2: 使用 AR(1) 解析反矩陣公式快速填值 (免迴圈)
  // =========================================================

  double rho_sq = rho * rho;
  double scale = 1.0 / (1.0 - rho_sq);

  double coef_mid_diag = (1.0 + rho_sq) * scale;
  double coef_off_diag = -rho * scale;

#pragma omp parallel for schedule(static)
  for (arma::uword i = 0; i < n_n; ++i) {

    arma::vec sigma_i = Sigma_mat.row(i).t();
    arma::mat & V_inv = V_inv_cube.slice(i);

    // 預先算出 1 / sigma 以及 1 / sigma^2
    arma::vec inv_sig = 1.0 / sigma_i;
    arma::vec inv_sig_sq = inv_sig % inv_sig;

    // 1. 填入主對角線 (Diagonal)
    // 預設全部填入中間項係數
    arma::vec diag_vals = coef_mid_diag * inv_sig_sq;
    // 修正頭尾兩項
    diag_vals(0) = scale * inv_sig_sq(0);
    diag_vals(n_t - 1) = scale * inv_sig_sq(n_t - 1);

    // 直接貼上主對角線
    V_inv.diag() = diag_vals;

    // 2. 填入次對角線 (Off-diagonal)
    if (n_t > 1) {
      // 次對角線的值: coef * (1/sig_j) * (1/sig_{j+1})
      arma::vec off_diag_vals = coef_off_diag * (inv_sig.subvec(0, n_t - 2) % inv_sig.subvec(1, n_t - 1));

      // 直接貼上上下兩條次對角線 (k=1 代表上對角線, k=-1 代表下對角線)
      V_inv.diag(1) = off_diag_vals;
      V_inv.diag(-1) = off_diag_vals;
    }
  }

  return V_inv_cube;
}



