#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// =========================================================================
// 1. 通用核心 (Generic Optimization Core)
//    這裡集中處理所有的 EL 計算與 Newton-Raphson 演算法，避免重複。
// =========================================================================

namespace EL_Core {

// 定義回傳結構
struct Result {

  double value;
  arma::vec gradient;
  arma::mat hessian;
};

// 通用 Lagrange 計算函數
// Psi_Gen: 一個函數物件 (Lambda)，輸入 index i，回傳該樣本的 Psi 向量
template <typename Psi_Gen>
List compute_lagrange(const arma::uword n,
                      const arma::vec & eta,
                      Psi_Gen get_psi) {

  double value = 0.0;
  const arma::uword dim = eta.n_elem;
  arma::vec gradient = arma::zeros(dim);
  arma::mat hessian = arma::zeros(dim, dim);

  for (size_t i = 0; i < n; ++i) {

    // 呼叫 Lambda 取得該模型的 Psi
    const arma::vec Psi_i = get_psi(i);
    const double denominator = 1.0 + arma::dot(eta, Psi_i);

    if (denominator > 0.0) {

      value += std::log(denominator);
      const arma::vec ratio = Psi_i / denominator;
      gradient += ratio;
      hessian -= ratio * ratio.t(); // 這裡是負的，因為是 Hessian

    } else {

      value -= 10000.0; // Penalty for violation
    }
  }

  // 檢查數值穩定性
  if (!std::isfinite(value) || !gradient.is_finite() || !hessian.is_finite()) {

    value = -10000.0 * n;
    gradient.zeros();
    hessian.zeros();
  }

  return List::create(

    Named("value") = value,
    Named("gradient") = gradient,
    Named("hessian") = hessian
  );
}

// 通用 Newton-Raphson 求解器
// Obj_Func: 計算目標函數的 Lambda (輸入 eta，回傳 List)
template <typename Obj_Func>
List solve_newton(arma::vec eta,
                  Obj_Func calc_obj,
                  const size_t & iter_max,
                  const double & step_rate,
                  const size_t & step_max,
                  const double & tol,
                  const double & eps_inv) {

  // 初始狀態
  List curr = calc_obj(eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double curr_val = curr["value"];
    const arma::vec grad = curr["gradient"];
    const arma::mat hess = curr["hessian"];

    // 計算方向
    arma::vec direction;
    if (arma::rcond(hess) > eps_inv) {

      direction = arma::solve(hess, grad);

    } else {

      direction = -grad;
    }

    // Backtracking Line Search
    double step_size = 1.0;
    arma::vec eta_new = eta - direction;
    List next = calc_obj(eta_new);

    bool improved = false;
    for (size_t s = 0; s < step_max; s++) {

      double next_val = next["value"];

      if (next_val <= curr_val + tol) {

        step_size /= step_rate;
        eta_new = eta - direction * step_size;
        next = calc_obj(eta_new);

      } else {

        improved = true;
        break;
      }
    }

    if (improved) {

      eta = eta_new;
      curr = next;

    } else {
      break; // 無法再優化
    }
  }

  return List::create(
    Named("eta") = eta,
    Named("value") = curr["value"],
    Named("gradient") = curr["gradient"],
    Named("hessian") = curr["hessian"]
  );
}

}

// =========================================================================
// 2. Normal Regression Models
// =========================================================================

// No shift, EY

arma::vec get_psi_AD_EY_normal_rcpp(const arma::vec & Xrow_i,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & phi) {

  // Psi = (alpha + X * beta) - phi
  double SI_i = alpha + arma::dot(Xrow_i, beta);
  return arma::vec({SI_i - phi});
}

// [[Rcpp::export]]
List SolveLagrange_AD_EY_normal_rcpp(const arma::mat & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const double & phi,
                                     const double & eta_initial,
                                     const size_t & iter_max,
                                     const double & step_rate,
                                     const size_t & step_max,
                                     const double & tol,
                                     const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  // 這裡我們定義一個 Lambda，它負責將整體的 eta 傳入，
  // 並利用 EL_Core::compute_lagrange 來對所有樣本求和。
  auto obj_func = [&](const arma::vec & e) {

    // 定義如何取得第 i 個樣本的 Psi
    auto get_psi_lambda = [&](size_t i) -> arma::vec {

      return get_psi_AD_EY_normal_rcpp(X.row(i).t(), alpha, beta, phi);
    };

    // 呼叫 EL Core 計算 Lagrange 值、梯度與 Hessian
    return EL_Core::compute_lagrange(n, arma::vec({e(0)}), get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(arma::vec({eta_initial}),
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_AD_EY_normal_rcpp(const arma::mat & X,
                                const double & alpha,
                                const arma::vec & beta,
                                const double & sigma,
                                const double & phi,
                                const double & eta,
                                const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword n_m = 1; // 明確定義 Psi 的維度

  // 1. 參數打包
  const arma::uword n_params = 1 + p + 1 + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, p) = beta;
  theta(p + 1) = sigma;
  theta(p + 2) = phi;

  // 2. 輔助函數
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {
    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, p);
    double curr_phi = curr_theta(p + 2);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {
      arma::vec psi_vec = get_psi_AD_EY_normal_rcpp(X.row(i).t(),
                                                    curr_alpha,
                                                    curr_beta,
                                                    curr_phi);

      double denom = 1.0 + eta * psi_vec(0);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi 和 mean_psi_outer
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec psi_vec = get_psi_AD_EY_normal_rcpp(X.row(i).t(), alpha, beta, phi);

    double denom = 1.0 + eta * psi_vec(0);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    mean_psi += psi_vec * inv_denom;
    mean_psi_outer += (psi_vec * psi_vec.t()) * inv_denom;
  }

  mean_psi /= n;
  mean_psi_outer /= n;

  // 4. ADvar
  arma::mat ADvar = mean_psi_outer + sigma * sigma;

  // 5. 數值微分計算 Gradient
  // *** 修改點：使用矩陣格式以保持一致性 ***
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    // 使用 .col(j) 填入，確保結果是 (1 x n_params) 的矩陣
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient,
    Named("ADvar") = ADvar
  );
}

// No shift, EXsubY

arma::vec get_psi_AD_EXsubY_normal_rcpp(const arma::vec & Xrow_i,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const double & sigma,
                                        const arma::mat & phi,     // phi: (n_k, n_p)
                                        const arma::mat & y_pts) { // y_pts: (n_k, 2)

  const arma::uword n_k = phi.n_rows;

  // 1. 計算線性預測與標準化座標
  double SI_i = alpha + arma::dot(Xrow_i, beta);
  arma::mat q_pts = (y_pts - SI_i) / sigma;

  // 2. 計算 CDF 差值 P(a < Y < b)
  arma::mat cdf_pts = arma::normcdf(q_pts);
  arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);

  // 3. 計算距離矩陣 (X - phi)
  // phi: (n_k, n_p)
  // Xrow_i.t(): (1, n_p) -> repmat to (n_k, n_p)
  arma::mat x_phi_dist = arma::repmat(Xrow_i.t(), n_k, 1) - phi;

  // 4. 計算 Psi
  // 公式: (X - phi) * P(a < Y < b)
  // 使用 .each_col() 讓 cdf_dist (n_k x 1) 自動廣播到每一列
  // 最後轉置並拉直成 (n_p * n_k, 1) 向量
  return arma::vectorise((x_phi_dist.each_col() % cdf_dist).t());
}

// [[Rcpp::export]]
List SolveLagrange_AD_EXsubY_normal_rcpp(const arma::mat & X,
                                         const double & alpha,
                                         const arma::vec & beta,
                                         const double & sigma,
                                         const arma::mat & phi,
                                         const arma::mat & y_pts,
                                         const arma::vec & eta_initial,
                                         const size_t & iter_max,
                                         const double & step_rate,
                                         const size_t & step_max,
                                         const double & tol,
                                         const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    // 定義如何取得第 i 個樣本的 Psi
    auto get_psi_lambda = [&](size_t i) -> arma::vec {

      return get_psi_AD_EXsubY_normal_rcpp(X.row(i).t(), alpha, beta, sigma, phi, y_pts);
    };

    // 呼叫 EL Core
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解 (eta 為向量)
  return EL_Core::solve_newton(eta_initial,
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_AD_EXsubY_normal_rcpp(const arma::mat & X,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & sigma,
                                    const arma::mat & phi,
                                    const arma::mat & y_pts,
                                    const arma::vec & eta,
                                    const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k; // Psi 的維度

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi...]
  const arma::uword n_params = 1 + n_p + 1 + n_m;
  arma::vec theta = arma::zeros(n_params);

  // 填入基礎參數
  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta.subvec(n_p + 2, n_params - 1) = arma::vectorise(phi.t());

  // 2. 輔助函數：計算 "加權" 平均 Psi (用於數值微分)
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    // 解包參數
    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_sigma = curr_theta(n_p + 1);
    arma::vec phi_vec = curr_theta.subvec(n_p + 2, n_params - 1);
    arma::mat curr_phi = arma::reshape(phi_vec, n_p, n_k).t();

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {
      arma::vec psi_vec = get_psi_AD_EXsubY_normal_rcpp(X.row(i).t(),
                                                        curr_alpha,
                                                        curr_beta,
                                                        curr_sigma,
                                                        curr_phi,
                                                        y_pts);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A = arma::zeros(n_m);
  arma::mat var_B = arma::zeros(n_m, n_m);

  for (arma::uword i = 0; i < n; ++i) {

    // 取得基準 Psi (使用原始參數)
    arma::vec psi_vec = get_psi_AD_EXsubY_normal_rcpp(X.row(i).t(),
                                                      alpha,
                                                      beta,
                                                      sigma,
                                                      phi,
                                                      y_pts);

    // 計算權重
    double denom = 1.0 + arma::dot(eta, psi_vec);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) 累積 mean_psi 與 mean_psi_outer
    mean_psi += psi_vec * inv_denom;
    mean_psi_outer += (psi_vec * psi_vec.t()) * inv_denom;

    // 準備計算 var_A 與 var_B 的中間變數
    double SI_i = alpha + arma::dot(X.row(i), beta);
    arma::mat q_pts = (y_pts - SI_i) / sigma;
    arma::mat cdf_pts = arma::normcdf(q_pts);
    arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);

    // 計算距離矩陣 (X - phi)
    arma::mat dist = (arma::repmat(X.row(i), n_k, 1) - phi);

    // (B) var_A (diff)
    arma::vec cdf_dist_expanded = arma::vectorise((arma::repmat(cdf_dist, 1, n_p)).t());
    var_A -= cdf_dist_expanded * inv_denom;

    // (C) var_B (Block Diagonal calculation)
    for (size_t k = 0; k < n_k; ++k) {
      arma::vec dist_k = dist.row(k).t();
      var_B.submat(k * n_p, k * n_p, (k + 1) * n_p - 1, (k + 1) * n_p - 1)
        += (dist_k * dist_k.t()) * cdf_dist(k) * inv_denom;
    }
  }

  mean_psi /= n;
  mean_psi_outer /= n;
  var_A /= n;
  var_B /= n;

  // 防呆機制：避免 var_A 元素過小導致除以零
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A(k)) < epsilon) {
      var_A(k) = (var_A(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  // ADvar 計算
  arma::vec var_A_inv = 1.0 / var_A;
  arma::mat ADvar = (var_B.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Cols: alpha, beta, sigma, phi
    Named("ADvar") = ADvar
  );
}

// No shift, EYsubX

arma::vec get_psi_AD_EYsubX_normal_rcpp(const arma::vec & Xrow_i,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const arma::vec & phi,
                                        const arma::vec & inclusion) {

  // 1. 計算線性預測值
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算 Psi
  // 公式: (SI - phi) * inclusion
  return (SI_i - phi) % inclusion;
}

// [[Rcpp::export]]
List SolveLagrange_AD_EYsubX_normal_rcpp(const arma::mat & X,
                                         const double & alpha,
                                         const arma::vec & beta,
                                         const arma::vec & phi,
                                         const arma::umat & inclusion,
                                         const arma::vec & eta_initial,
                                         const size_t & iter_max,
                                         const double & step_rate,
                                         const size_t & step_max,
                                         const double & tol,
                                         const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    // 定義如何取得第 i 個樣本的 Psi
    auto get_psi_lambda = [&](size_t i) -> arma::vec {

      return get_psi_AD_EYsubX_normal_rcpp(
        X.row(i).t(),
        alpha,
        beta,
        phi,
        arma::conv_to<arma::vec>::from(inclusion.row(i).t()));
    };

    // 呼叫 EL Core
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial,
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_AD_EYsubX_normal_rcpp(const arma::mat & X,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & sigma,
                                    const arma::vec & phi,
                                    const arma::umat & inclusion,
                                    const arma::vec & eta,
                                    const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi...]
  const arma::uword n_params = 1 + n_p + 1 + n_m;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta.subvec(n_p + 2, n_params - 1) = phi;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    // 解包參數
    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec curr_phi = curr_theta.subvec(n_p + 2, n_params - 1);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
      arma::vec psi_vec = get_psi_AD_EYsubX_normal_rcpp(X.row(i).t(),
                                                        curr_alpha,
                                                        curr_beta,
                                                        curr_phi,
                                                        inc_vec);

      double denom = 1.0 + arma::dot(eta, psi_vec);

      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A = arma::zeros(n_m);
  arma::mat var_B_a = arma::zeros(n_m, n_m);

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
    arma::vec psi_vec = get_psi_AD_EYsubX_normal_rcpp(X.row(i).t(),
                                                      alpha,
                                                      beta,
                                                      phi,
                                                      inc_vec);

    double denom = 1.0 + arma::dot(eta, psi_vec);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    mean_psi += psi_vec * inv_denom;
    mean_psi_outer += (psi_vec * psi_vec.t()) * inv_denom;

    var_A -= inc_vec * inv_denom;
    var_B_a += (inc_vec * inc_vec.t()) * inv_denom;
  }

  mean_psi /= n;
  mean_psi_outer /= n;
  var_A /= n;
  var_B_a /= n;

  arma::mat var_B = var_B_a * sigma * sigma + mean_psi_outer;

  // 防呆機制
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A(k)) < epsilon) {
      var_A(k) = (var_A(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_inv = 1.0 / var_A;
  arma::mat ADvar = (var_B.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Cols: alpha, beta..., sigma, phi...
    Named("ADvar") = ADvar
  );
}

// Covariate shift, EY

arma::vec get_psi_ADCS_EY_normal_rcpp(const arma::vec & Xrow_i,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & phi,
                                      const arma::vec & CS_beta) {

  // 1. 計算線性預測值
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算 Covariate Shift 權重 (含數值截斷)
  double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
  if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

  // 3. 計算 Psi
  // 公式: w(X) * (SI - phi)
  return arma::vec({eSI_CS_i * (SI_i - phi)});
}

// [[Rcpp::export]]
List SolveLagrange_ADCS_EY_normal_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & phi,
                                       const arma::vec & CS_beta,
                                       const double & eta_initial,
                                       const size_t & iter_max,
                                       const double & step_rate,
                                       const size_t & step_max,
                                       const double & tol,
                                       const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    // 定義如何取得第 i 個樣本的 Psi
    auto get_psi_lambda = [&](size_t i) -> arma::vec {

      return get_psi_ADCS_EY_normal_rcpp(X.row(i).t(), alpha, beta, phi, CS_beta);
    };

    // 呼叫 EL Core (e(0) 是 scalar eta)
    return EL_Core::compute_lagrange(n, arma::vec({e(0)}), get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(arma::vec({eta_initial}),
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADCS_EY_normal_rcpp(const arma::mat & X,
                                  const double & alpha,
                                  const arma::vec & beta,
                                  const double & sigma,
                                  const double & phi,
                                  const arma::vec & CS_beta,
                                  const double & eta,
                                  const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = 1; // Psi 維度

  // 1. 參數打包 theta = [alpha(1), beta(p), sigma(1), phi(1), CS_beta(p)]
  const arma::uword n_params = 1 + n_p + 1 + 1 + n_p;
  arma::vec theta = arma::zeros(n_params);

  // 填入參數
  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta(n_p + 2) = phi;
  theta.subvec(n_p + 3, n_params - 1) = CS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_phi = curr_theta(n_p + 2);
    arma::vec curr_CS_beta = curr_theta.subvec(n_p + 3, n_params - 1);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADCS_EY_normal_rcpp(X.row(i).t(),
                                                      curr_alpha,
                                                      curr_beta,
                                                      curr_phi,
                                                      curr_CS_beta);

      double denom = 1.0 + eta * psi_vec(0);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }

    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  double var_B_Ntor = 0.0;
  double var_B_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {

    arma::vec Xrow_i = X.row(i).t();
    double SI_i = alpha + arma::dot(Xrow_i, beta);
    double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

    double Psi_val = eSI_CS_i * (SI_i - phi);

    // 計算權重
    double denom = 1.0 + eta * Psi_val;
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    mean_psi(0) += Psi_val * inv_denom;
    mean_psi_outer(0,0) += (Psi_val * Psi_val) * inv_denom;

    var_B_Ntor += (SI_i - phi) * (SI_i - phi) * eSI_CS_i * inv_denom;
    var_B_Dtor += eSI_CS_i * inv_denom;
  }

  mean_psi /= n;
  mean_psi_outer /= n;

  // ADvar
  arma::vec ADvar = {0.0};
  if (var_B_Dtor < 1e-10) var_B_Dtor = 1e-10;
  ADvar = (var_B_Ntor / var_B_Dtor) + sigma * sigma;

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi, CS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, sigma, phi, CS_beta
    Named("ADvar") = ADvar
  );
}

// Covariate shift, EXsubY

arma::vec get_psi_ADCS_EXsubY_normal_rcpp(const arma::vec & Xrow_i,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const double & sigma,
                                          const arma::mat & phi,     // phi: (n_k, n_p)
                                          const arma::vec & CS_beta,
                                          const arma::mat & y_pts) {

  const arma::uword n_k = phi.n_rows;

  // 1. 計算線性預測
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算 Covariate Shift 權重 (含截斷)
  double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
  if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

  // 3. 計算機率區間
  arma::mat q_pts = (y_pts - SI_i) / sigma;
  arma::mat cdf_pts = arma::normcdf(q_pts);
  arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);

  // 4. 計算距離矩陣 (X - phi)
  arma::mat x_phi_dist = arma::repmat(Xrow_i.t(), n_k, 1) - phi;

  // 5. 計算 Psi
  // 公式: (X - phi) * cdf_dist * w(X)
  return arma::vectorise((x_phi_dist.each_col() % cdf_dist).t()) * eSI_CS_i;
}

// [[Rcpp::export]]
List SolveLagrange_ADCS_EXsubY_normal_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const double & sigma,
                                           const arma::mat & phi,
                                           const arma::vec & CS_beta,
                                           const arma::mat & y_pts,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol,
                                           const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {

    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADCS_EXsubY_normal_rcpp(X.row(i).t(),
                                             alpha,
                                             beta,
                                             sigma,
                                             phi,
                                             CS_beta, y_pts);
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial,
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADCS_EXsubY_normal_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & sigma,
                                      const arma::mat & phi,
                                      const arma::vec & CS_beta,
                                      const arma::mat & y_pts,
                                      const arma::vec & eta,
                                      const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k; // Psi 的維度

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi..., CS_beta...]
  const arma::uword n_params = 1 + n_p + 1 + n_m + n_p;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta.subvec(n_p + 2, n_p + 1 + n_m) = arma::vectorise(phi.t());
  theta.subvec(n_p + 2 + n_m, n_params - 1) = CS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_sigma = curr_theta(n_p + 1);
    arma::vec phi_vec = curr_theta.subvec(n_p + 2, n_p + 1 + n_m);
    arma::mat curr_phi = arma::reshape(phi_vec, n_p, n_k).t();
    arma::vec curr_CS_beta = curr_theta.subvec(n_p + 2 + n_m, n_params - 1);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADCS_EXsubY_normal_rcpp(X.row(i).t(),
                                                          curr_alpha,
                                                          curr_beta,
                                                          curr_sigma,
                                                          curr_phi,
                                                          curr_CS_beta,
                                                          y_pts);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {

    arma::vec Xrow_i = X.row(i).t();
    double SI_i = alpha + arma::dot(Xrow_i, beta);
    double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

    arma::mat q_pts = (y_pts - SI_i) / sigma;
    arma::mat cdf_pts = arma::normcdf(q_pts);
    arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    arma::mat x_phi_dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);

    // Psi
    arma::vec Psi_i = arma::vectorise((x_phi_dist.each_col() % cdf_dist).t()) * eSI_CS_i;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 自歸一化權重
    double weighted_w = eSI_CS_i * inv_denom;
    var_AB_Dtor += weighted_w;

    // (C) var_A (diff)
    arma::vec cdf_dist_expanded = arma::vectorise((arma::repmat(cdf_dist, 1, n_p)).t());
    var_A_accum -= cdf_dist_expanded * weighted_w;

    // (D) var_B (var)
    for (size_t k = 0; k < n_k; ++k) {
      arma::vec dist_k = x_phi_dist.row(k).t();
      var_B_accum.submat(k * n_p, k * n_p, (k + 1) * n_p - 1, (k + 1) * n_p - 1)
        += (dist_k * dist_k.t()) * cdf_dist(k) * weighted_w;
    }
  }

  mean_psi /= n;
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_AB_Dtor < 1e-10) var_AB_Dtor = 1e-10;
  var_A_accum /= var_AB_Dtor;
  var_B_accum /= var_AB_Dtor;

  // 防呆機制
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B_accum.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi, CS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, sigma, phi, CS_beta
    Named("ADvar") = ADvar
  );
}

// Covariate shift, EYsubX

arma::vec get_psi_ADCS_EYsubX_normal_rcpp(const arma::vec & Xrow_i,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const arma::vec & phi,
                                          const arma::vec & CS_beta,
                                          const arma::vec & inclusion) {

  // 1. 計算線性預測
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算 Covariate Shift 權重 (含截斷)
  double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
  if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

  // 3. 計算 Psi
  // 公式: (SI - phi) * inclusion * w(X)
  return ((SI_i - phi) % inclusion) * eSI_CS_i;
}

// [[Rcpp::export]]
List SolveLagrange_ADCS_EYsubX_normal_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const arma::vec & phi,
                                           const arma::vec & CS_beta,
                                           const arma::umat & inclusion,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol,
                                           const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {

      return get_psi_ADCS_EYsubX_normal_rcpp(
        X.row(i).t(),
        alpha,
        beta,
        phi,
        CS_beta,
        arma::conv_to<arma::vec>::from(inclusion.row(i).t()));
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial,
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADCS_EYsubX_normal_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & sigma,
                                      const arma::vec & phi,
                                      const arma::vec & CS_beta,
                                      const arma::umat & inclusion,
                                      const arma::vec & eta,
                                      const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi..., CS_beta...]
  const arma::uword n_params = 1 + n_p + 1 + n_m + n_p;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta.subvec(n_p + 2, n_p + 1 + n_m) = phi;
  theta.subvec(n_p + 2 + n_m, n_params - 1) = CS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec curr_phi = curr_theta.subvec(n_p + 2, n_p + 1 + n_m);
    arma::vec curr_CS_beta = curr_theta.subvec(n_p + 2 + n_m, n_params - 1);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
      arma::vec psi_vec = get_psi_ADCS_EYsubX_normal_rcpp(X.row(i).t(),
                                                          curr_alpha,
                                                          curr_beta,
                                                          curr_phi,
                                                          curr_CS_beta,
                                                          inc_vec);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_a_accum = arma::zeros(n_m, n_m);
  arma::mat var_B_b_accum = arma::zeros(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {

    arma::vec Xrow_i = X.row(i).t();
    arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());

    double SI_i = alpha + arma::dot(Xrow_i, beta);
    double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

    // Psi
    arma::vec Psi_diff_i = (SI_i - phi) % inc_vec;
    arma::vec Psi_i = Psi_diff_i * eSI_CS_i;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;                 // [新增] 累積加權平均
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 自歸一化權重
    double weighted_w = eSI_CS_i * inv_denom;
    var_AB_Dtor += weighted_w;

    // (C) var_A (diff)
    var_A_accum -= inc_vec * weighted_w;

    // (D) var_B (var)
    var_B_a_accum += (inc_vec * inc_vec.t()) * weighted_w;
    var_B_b_accum += (Psi_diff_i * Psi_diff_i.t()) * weighted_w;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_AB_Dtor < 1e-10) var_AB_Dtor = 1e-10;
  var_A_accum /= var_AB_Dtor;
  var_B_a_accum /= var_AB_Dtor;
  var_B_b_accum /= var_AB_Dtor;

  // 組合 var_B
  arma::mat var_B = var_B_a_accum * sigma * sigma + var_B_b_accum;

  // 防呆機制
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_accum_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B.each_col() % var_A_accum_inv).each_row() % var_A_accum_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi, CS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, sigma, phi, CS_beta
    Named("ADvar") = ADvar
  );
}

// Prior probability shift, EX

arma::vec get_psi_ADPPS_EX_normal_rcpp(const arma::vec & Xrow_i,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & sigma,
                                       const arma::vec & phi,
                                       const double & PPS_beta) {

  // 1. 計算線性預測
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算 PPS 權重
  // 公式: exp( (2*SI + sigma^2*PPS_beta) * PPS_beta / 2 )
  // 展開: exp( SI*PPS_beta + 0.5*sigma^2*PPS_beta^2 )
  double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
  if (e_f > 1e10) e_f = 1e10; // 加上基礎防呆

  // 3. 計算 Psi
  // 公式: (X - phi) * w(X)
  return (Xrow_i - phi) * e_f;
}

// [[Rcpp::export]]
List SolveLagrange_ADPPS_EX_normal_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const double & sigma,
                                        const arma::vec & phi,
                                        const double & PPS_beta,
                                        const arma::vec & eta_initial,
                                        const size_t & iter_max,
                                        const double & step_rate,
                                        const size_t & step_max,
                                        const double & tol,
                                        const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADPPS_EX_normal_rcpp(X.row(i).t(), alpha, beta, sigma, phi, PPS_beta);
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADPPS_EX_normal_rcpp(const arma::mat & X,
                                   const double & alpha,
                                   const arma::vec & beta,
                                   const double & sigma,
                                   const arma::vec & phi,
                                   const double & PPS_beta,
                                   const arma::vec & eta,
                                   const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem; // 應等於 n_p

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi..., PPS_beta]
  const arma::uword n_params = 1 + n_p + 1 + n_m + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta.subvec(n_p + 2, n_p + 1 + n_m) = phi;
  theta(n_p + 2 + n_m) = PPS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_sigma = curr_theta(n_p + 1);
    arma::vec curr_phi = curr_theta.subvec(n_p + 2, n_p + 1 + n_m);
    double curr_PPS_beta = curr_theta(n_p + 2 + n_m);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADPPS_EX_normal_rcpp(X.row(i).t(),
                                                       curr_alpha,
                                                       curr_beta,
                                                       curr_sigma,
                                                       curr_phi,
                                                       curr_PPS_beta);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);
  double var_B_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec Xrow_i = X.row(i).t();

    double SI_i = alpha + arma::dot(Xrow_i, beta);
    double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
    if (e_f > 1e10) e_f = 1e10;

    // Psi
    arma::vec Psi_i = (Xrow_i - phi) * e_f;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;                 // [新增] 累積加權平均
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 自歸一化權重
    double weighted_w = e_f * inv_denom;

    var_B_Dtor += weighted_w;
    arma::vec diff = (Xrow_i - phi);
    var_B_accum += (diff * diff.t()) * weighted_w;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;

  // Normalization
  if (var_B_Dtor < 1e-10) var_B_Dtor = 1e-10;
  var_B_accum /= var_B_Dtor;

  arma::mat ADvar = var_B_accum;

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi, PPS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, sigma, phi, PPS_beta
    Named("ADvar") = ADvar
  );
}

// Prior probability shift, EY

arma::vec get_psi_ADPPS_EY_normal_rcpp(const arma::vec & Xrow_i,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & sigma,
                                       const double & phi,
                                       const double & PPS_beta) {

  // 1. 計算線性預測
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算偏移後的中心 (Shifted Center)
  // 公式: SI + sigma^2 * beta_pps - phi
  double c_f = SI_i + sigma * sigma * PPS_beta - phi;

  // 3. 計算 PPS 權重
  // 公式: exp( (2*SI + sigma^2 * beta_pps) * beta_pps / 2 )
  double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
  if (e_f > 1e10) e_f = 1e10;

  // 4. 計算 Psi = c_f * w
  return arma::vec({c_f * e_f});
}

// [[Rcpp::export]]
List SolveLagrange_ADPPS_EY_normal_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const double & sigma,
                                        const double & phi,
                                        const double & PPS_beta,
                                        const double & eta_initial,
                                        const size_t & iter_max,
                                        const double & step_rate,
                                        const size_t & step_max,
                                        const double & tol,
                                        const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      // eta 是 scalar，但介面需要 vec
      return get_psi_ADPPS_EY_normal_rcpp(X.row(i).t(), alpha, beta, sigma, phi, PPS_beta);
    };
    // 傳入 e(0) 作為 scalar eta，compute_lagrange 會將其轉為 vec
    return EL_Core::compute_lagrange(n, arma::vec({e(0)}), get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(arma::vec({eta_initial}),
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADPPS_EY_normal_rcpp(const arma::mat & X,
                                   const double & alpha,
                                   const arma::vec & beta,
                                   const double & sigma,
                                   const double & phi,
                                   const double & PPS_beta,
                                   const double & eta,
                                   const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = 1; // Psi 維度

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi, PPS_beta]
  const arma::uword n_params = 1 + n_p + 1 + 1 + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta(n_p + 2) = phi;
  theta(n_p + 3) = PPS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi (用於差分計算 Gradient)
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_sigma = curr_theta(n_p + 1);
    double curr_phi = curr_theta(n_p + 2);
    double curr_PPS_beta = curr_theta(n_p + 3);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      // 使用 curr_phi 與 curr_PPS_beta 進行計算
      arma::vec psi_vec = get_psi_ADPPS_EY_normal_rcpp(X.row(i).t(),
                                                       curr_alpha,
                                                       curr_beta,
                                                       curr_sigma,
                                                       curr_phi,
                                                       curr_PPS_beta);

      double denom = 1.0 + eta * psi_vec(0);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }

    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);

  double var_B_Ntor = 0.0;
  double var_B_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    // 準備數據
    arma::vec Xrow_i = X.row(i).t();

    double SI_i = alpha + arma::dot(Xrow_i, beta);
    double c_f = SI_i + sigma * sigma * PPS_beta - phi;

    // PPS 權重
    double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
    if (e_f > 1e10) e_f = 1e10;

    // Psi = c_f * w
    double Psi_val = c_f * e_f;

    // EL 權重計算
    double denom = 1.0 + eta * Psi_val;
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi(0) += Psi_val * inv_denom;            // [新增] 累積加權平均
    mean_psi_outer(0,0) += (Psi_val * Psi_val) * inv_denom;

    // (B) 準備自歸一化變異數計算
    double weighted_w = e_f * inv_denom;

    // 分子: c_f^2 * w * p
    var_B_Ntor += (c_f * c_f) * weighted_w;

    // 分母: w * p
    var_B_Dtor += weighted_w;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;

  // 自歸一化 normalization 並加上 sigma^2
  double var_B_val = 0.0;
  if (var_B_Dtor < 1e-10) var_B_Dtor = 1e-10;

  var_B_val = (var_B_Ntor / var_B_Dtor) + sigma * sigma;

  // 組合 ADvar
  arma::vec ADvar = {var_B_val};

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi, PPS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, sigma, phi, PPS_beta
    Named("ADvar") = ADvar
  );
}

// Prior probability shift, EXsubY

arma::vec get_psi_ADPPS_EXsubY_normal_rcpp(const arma::vec & Xrow_i,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const double & sigma,
                                           const arma::mat & phi,     // n_k * n_p
                                           const double & PPS_beta,
                                           const arma::mat & y_pts) {

  const arma::uword n_k = phi.n_rows;

  // 1. 計算線性預測
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算 PPS 權重
  // 公式: exp( (2*SI + sigma^2 * PPS_beta) * PPS_beta / 2 )
  double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
  if (e_f > 1e10) e_f = 1e10;

  // 3. 計算機率區間 (注意中心點的修正)
  // Shifted mean = SI + sigma^2 * PPS_beta
  double mean_shift = SI_i + sigma * sigma * PPS_beta;
  arma::mat q_pts = (y_pts - mean_shift) / sigma;
  arma::mat cdf_pts = arma::normcdf(q_pts);
  arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);

  // 4. 計算距離矩陣 (X - phi)
  arma::mat x_phi_dist = arma::repmat(Xrow_i.t(), n_k, 1) - phi;

  // 5. 計算 Psi
  // 公式: (X - phi) * cdf_dist * w(X)
  return arma::vectorise((x_phi_dist.each_col() % cdf_dist).t()) * e_f;
}

// [[Rcpp::export]]
List SolveLagrange_ADPPS_EXsubY_normal_rcpp(const arma::mat & X,
                                            const double & alpha,
                                            const arma::vec & beta,
                                            const double & sigma,
                                            const arma::mat & phi,
                                            const double & PPS_beta,
                                            const arma::mat & y_pts,
                                            const arma::vec & eta_initial,
                                            const size_t & iter_max,
                                            const double & step_rate,
                                            const size_t & step_max,
                                            const double & tol,
                                            const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {

    auto get_psi_lambda = [&](size_t i) -> arma::vec {

      return get_psi_ADPPS_EXsubY_normal_rcpp(X.row(i).t(), alpha, beta, sigma, phi, PPS_beta, y_pts);
    };

    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADPPS_EXsubY_normal_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & sigma,
                                       const arma::mat & phi,
                                       const double & PPS_beta,
                                       const arma::mat & y_pts,
                                       const arma::vec & eta,
                                       const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi..., PPS_beta]
  const arma::uword n_params = 1 + n_p + 1 + n_m + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta.subvec(n_p + 2, n_p + 1 + n_m) = arma::vectorise(phi.t());
  theta(n_p + 2 + n_m) = PPS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_sigma = curr_theta(n_p + 1);
    arma::vec phi_vec = curr_theta.subvec(n_p + 2, n_p + 1 + n_m);
    arma::mat curr_phi = arma::reshape(phi_vec, n_p, n_k).t();
    double curr_PPS_beta = curr_theta(n_p + 2 + n_m);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADPPS_EXsubY_normal_rcpp(X.row(i).t(),
                                                           curr_alpha,
                                                           curr_beta,
                                                           curr_sigma,
                                                           curr_phi,
                                                           curr_PPS_beta,
                                                           y_pts);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    // 準備數據
    arma::vec Xrow_i = X.row(i).t();
    double SI_i = alpha + arma::dot(Xrow_i, beta);

    // PPS 權重
    double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
    if (e_f > 1e10) e_f = 1e10;

    // 計算機率區間
    double mean_shift = SI_i + sigma * sigma * PPS_beta;
    arma::mat q_pts = (y_pts - mean_shift) / sigma;
    arma::mat cdf_pts = arma::normcdf(q_pts);
    arma::vec cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    arma::mat x_phi_dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);

    // Psi
    arma::vec Psi_i = arma::vectorise((x_phi_dist.each_col() % cdf_dist).t()) * e_f;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;                 // [新增] 累積加權平均
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 準備自歸一化變異數計算
    double weighted_w = e_f * inv_denom;
    var_AB_Dtor += weighted_w;

    // (C) var_A (diff)
    arma::vec cdf_dist_expanded = arma::vectorise((arma::repmat(cdf_dist, 1, n_p)).t());
    var_A_accum -= cdf_dist_expanded * weighted_w;

    // (D) var_B (var)
    for (size_t k = 0; k < n_k; ++k) {
      arma::vec dist_k = x_phi_dist.row(k).t();
      var_B_accum.submat(k * n_p, k * n_p, (k + 1) * n_p - 1, (k + 1) * n_p - 1)
        += (dist_k * dist_k.t()) * cdf_dist(k) * weighted_w;
    }
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_AB_Dtor < 1e-10) var_AB_Dtor = 1e-10;
  var_A_accum /= var_AB_Dtor;
  var_B_accum /= var_AB_Dtor;

  // 防呆機制
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_accum_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B_accum.each_col() % var_A_accum_inv).each_row() % var_A_accum_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi, PPS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, sigma, phi, PPS_beta
    Named("ADvar") = ADvar
  );
}

// Prior probability shift, EYsubX

arma::vec get_psi_ADPPS_EYsubX_normal_rcpp(const arma::vec & Xrow_i,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const double & sigma,
                                           const arma::vec & phi,
                                           const double & PPS_beta,
                                           const arma::vec & inclusion) {

  // 1. 計算線性預測
  double SI_i = alpha + arma::dot(Xrow_i, beta);

  // 2. 計算偏移後的均值向量 (Shifted Mean Vector)
  // 公式: SI + sigma^2 * PPS_beta - phi
  // 注意: SI_i 是 scalar, phi 是 vector, 結果為 vector
  arma::vec c_f = (SI_i + sigma * sigma * PPS_beta) - phi;

  // 3. 計算 PPS 權重
  // 公式: exp( (2*SI + sigma^2 * PPS_beta) * PPS_beta / 2 )
  double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
  if (e_f > 1e10) e_f = 1e10;

  // 4. 組合並回傳 Psi
  // (Shifted Mean - Phi) * I(Inclusion) * Weight
  return (c_f % inclusion) * e_f;
}

// [[Rcpp::export]]
List SolveLagrange_ADPPS_EYsubX_normal_rcpp(const arma::mat & X,
                                            const double & alpha,
                                            const arma::vec & beta,
                                            const double & sigma,
                                            const arma::vec & phi,
                                            const double & PPS_beta,
                                            const arma::umat & inclusion,
                                            const arma::vec & eta_initial,
                                            const size_t & iter_max,
                                            const double & step_rate,
                                            const size_t & step_max,
                                            const double & tol,
                                            const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADPPS_EYsubX_normal_rcpp(X.row(i).t(),
                                              alpha,
                                              beta,
                                              sigma,
                                              phi,
                                              PPS_beta,
                                              arma::conv_to<arma::vec>::from(inclusion.row(i).t()));
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADPPS_EYsubX_normal_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & sigma,
                                       const arma::vec & phi,
                                       const double & PPS_beta,
                                       const arma::umat & inclusion,
                                       const arma::vec & eta,
                                       const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  // 1. 參數打包 theta = [alpha, beta..., sigma, phi..., PPS_beta]
  const arma::uword n_params = 1 + n_p + 1 + n_m + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = sigma;
  theta.subvec(n_p + 2, n_p + 1 + n_m) = phi;
  theta(n_p + 2 + n_m) = PPS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_sigma = curr_theta(n_p + 1);
    arma::vec curr_phi = curr_theta.subvec(n_p + 2, n_p + 1 + n_m);
    double curr_PPS_beta = curr_theta(n_p + 2 + n_m);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
      arma::vec psi_vec = get_psi_ADPPS_EYsubX_normal_rcpp(X.row(i).t(),
                                                           curr_alpha,
                                                           curr_beta,
                                                           curr_sigma,
                                                           curr_phi,
                                                           curr_PPS_beta,
                                                           inc_vec);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_a_accum = arma::zeros(n_m, n_m);
  arma::mat var_B_b_accum = arma::zeros(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    // 準備數據
    arma::vec Xrow_i = X.row(i).t();
    arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());

    double SI_i = alpha + arma::dot(Xrow_i, beta);

    // PPS 權重
    double e_f = std::exp((SI_i * 2.0 + sigma * sigma * PPS_beta) * PPS_beta / 2.0);
    if (e_f > 1e10) e_f = 1e10;

    // Shifted center: (SI + sigma^2 * PPS_beta) - phi
    arma::vec c_f = (SI_i + sigma * sigma * PPS_beta) - phi;

    // Psi
    arma::vec Psi_diff_i = c_f % inc_vec;
    arma::vec Psi_i = Psi_diff_i * e_f;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;                 // [新增] 累積加權平均
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 準備自歸一化變異數計算
    double weighted_w = e_f * inv_denom;
    var_AB_Dtor += weighted_w;

    // (C) var_A (diff)
    var_A_accum -= inc_vec * weighted_w;

    // (D) var_B (var)
    var_B_a_accum += (inc_vec * inc_vec.t()) * weighted_w;
    var_B_b_accum += (Psi_diff_i * Psi_diff_i.t()) * weighted_w;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_AB_Dtor < 1e-10) var_AB_Dtor = 1e-10;
  var_A_accum /= var_AB_Dtor;
  var_B_a_accum /= var_AB_Dtor;
  var_B_b_accum /= var_AB_Dtor;

  // 組合 var_B
  arma::mat var_B = var_B_a_accum * sigma * sigma + var_B_b_accum;

  // 防呆機制
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_accum_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B.each_col() % var_A_accum_inv).each_row() % var_A_accum_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, sigma, phi, PPS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, sigma, phi, PPS_beta
    Named("ADvar") = ADvar
  );
}

// =========================================================================
// 3. Logistic Regression Models
// =========================================================================

// No shift, EY

arma::vec get_psi_AD_EY_logistic_rcpp(const arma::vec & Xrow_i,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & phi) {

  // 1. 計算指數項 e^(alpha + X*beta)
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  if (eSI_i > 1e10) eSI_i = 1e10;

  // 2. 計算 Psi
  // 公式: ( (1-phi)*eSI - phi ) / (1 + eSI)
  // 這等價於 P(Y=1|X) - phi
  double val = ((1.0 - phi) * eSI_i - phi) / (1.0 + eSI_i);

  return arma::vec({ val });
}

// [[Rcpp::export]]
List SolveLagrange_AD_EY_logistic_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & phi,
                                       const double & eta_initial,
                                       const size_t & iter_max,
                                       const double & step_rate,
                                       const size_t & step_max,
                                       const double & tol,
                                       const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      // eta 是 scalar
      return get_psi_AD_EY_logistic_rcpp(X.row(i).t(), alpha, beta, phi);
    };
    return EL_Core::compute_lagrange(n, arma::vec({e(0)}), get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(arma::vec({eta_initial}),
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_AD_EY_logistic_rcpp(const arma::mat & X,
                                  const double & alpha,
                                  const arma::vec & beta,
                                  const double & phi,
                                  const double & eta,
                                  const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = 1;

  // 1. 參數打包 theta = [alpha, beta..., phi]
  const arma::uword n_params = 1 + n_p + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = phi;

  // 2. 輔助函數：計算 "加權" 平均 Psi (用於差分計算 Gradient)
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_phi = curr_theta(n_p + 1);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_AD_EY_logistic_rcpp(X.row(i).t(),
                                                      curr_alpha,
                                                      curr_beta,
                                                      curr_phi);

      double denom = 1.0 + eta * psi_vec(0);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  double var_B_accum = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec Xrow_i = X.row(i).t();
    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    // Psi
    double Psi_val = ((1.0 - phi) * eSI_i - phi) / (1.0 + eSI_i);

    // 權重計算
    double denom = 1.0 + eta * Psi_val;
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi(0) += Psi_val * inv_denom;            // [新增] 累積加權平均
    mean_psi_outer(0,0) += (Psi_val * Psi_val) * inv_denom;

    // (B) var_B 計算
    double term_numerator = (1.0 - phi) * (1.0 - phi) * eSI_i + phi * phi;
    double term = term_numerator / (1.0 + eSI_i);

    var_B_accum += term * inv_denom;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;
  var_B_accum /= n;

  // 組合 ADvar
  arma::vec ADvar = {var_B_accum};

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi
    Named("ADvar") = ADvar
  );
}

// No shift, EXsubY

arma::vec get_psi_AD_EXsubY_logistic_rcpp(const arma::vec & Xrow_i,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const arma::mat & phi) {

  // phi 為 2 x p 矩陣
  // Row 0 對應 Y=1, Row 1 對應 Y=0

  // 1. 計算機率
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  // 截斷防護
  if (eSI_i > 1e10) eSI_i = 1e10;

  double prob_1 = eSI_i / (1.0 + eSI_i); // P(Y=1)
  double prob_0 = 1.0 / (1.0 + eSI_i);   // P(Y=0)

  // 2. 計算 Psi 分量
  // Part 1: (X - phi_1) * P(Y=1)
  arma::vec psi_1 = (Xrow_i - phi.row(0).t()) * prob_1;

  // Part 2: (X - phi_0) * P(Y=0)
  arma::vec psi_0 = (Xrow_i - phi.row(1).t()) * prob_0;

  // 3. 組合 (長度 2p)
  return arma::join_cols(psi_1, psi_0);
}

// [[Rcpp::export]]
List SolveLagrange_AD_EXsubY_logistic_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const arma::mat & phi,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol,
                                           const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_AD_EXsubY_logistic_rcpp(X.row(i).t(), alpha, beta, phi);
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_AD_EXsubY_logistic_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const arma::mat & phi,
                                      const arma::vec & eta,
                                      const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;

  // 1. 參數打包 theta = [alpha, beta..., phi...]
  const arma::uword n_params = 1 + n_p + n_m;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta.subvec(n_p + 1, n_params - 1) = arma::vectorise(phi.t());

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec phi_vec = curr_theta.subvec(n_p + 1, n_params - 1);
    arma::mat curr_phi = arma::reshape(phi_vec, n_p, n_k).t();
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_AD_EXsubY_logistic_rcpp(X.row(i).t(),
                                                          curr_alpha,
                                                          curr_beta,
                                                          curr_phi);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }

    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec Xrow_i = X.row(i).t();
    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    double prob_1 = eSI_i / (1.0 + eSI_i);
    double prob_0 = 1.0 / (1.0 + eSI_i);

    // Psi (n_m x 1)
    arma::vec psi_1 = (Xrow_i - phi.row(0).t()) * prob_1;
    arma::vec psi_0 = (Xrow_i - phi.row(1).t()) * prob_0;
    arma::vec Psi_i = arma::join_cols(psi_1, psi_0);

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;                 // [新增] 累積加權平均
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) var_A (diff)
    var_A_accum.subvec(0, n_p - 1) -= prob_1 * inv_denom;
    var_A_accum.subvec(n_p, 2 * n_p - 1) -= prob_0 * inv_denom;

    // (C) var_B (var)
    arma::vec diff_1 = (Xrow_i - phi.row(0).t());
    arma::vec diff_0 = (Xrow_i - phi.row(1).t());

    // Block 1 (Top-Left)
    var_B_accum.submat(0, 0, n_p - 1, n_p - 1)
      += (diff_1 * diff_1.t()) * prob_1 * inv_denom;
    // Block 2 (Bottom-Right)
    var_B_accum.submat(n_p, n_p, 2 * n_p - 1, 2 * n_p - 1)
      += (diff_0 * diff_0.t()) * prob_0 * inv_denom;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;
  var_A_accum /= n;
  var_B_accum /= n;

  // 防呆與 ADvar 計算
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B_accum.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi(row-major)
    Named("ADvar") = ADvar
  );
}

// No shift, EYsubX

arma::vec get_psi_AD_EYsubX_logistic_rcpp(const arma::vec & Xrow_i,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const arma::vec & phi,
                                          const arma::vec & inclusion) {

  // 1. 計算 Logistic Probability P(Y=1|X)
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  if (eSI_i > 1e10) eSI_i = 1e10; // 截斷

  // prob_1 = eSI / (1 + eSI)
  double prob_1 = eSI_i / (1.0 + eSI_i);

  // 2. 計算 Psi 向量 (未遮罩)
  // Psi = P(Y=1|X) - phi
  arma::vec psi_raw = prob_1 - phi;

  // 3. 應用遮罩
  // element-wise multiplication
  return psi_raw % inclusion;
}

// [[Rcpp::export]]
List SolveLagrange_AD_EYsubX_logistic_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const arma::vec & phi,
                                           const arma::umat & inclusion,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol,
                                           const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      // 將 inclusion 的 row 轉為 vec
      return get_psi_AD_EYsubX_logistic_rcpp(X.row(i).t(),
                                             alpha,
                                             beta,
                                             phi,
                                             arma::conv_to<arma::vec>::from(inclusion.row(i).t()));
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_AD_EYsubX_logistic_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const arma::vec & phi,
                                      const arma::umat & inclusion,
                                      const arma::vec & eta,
                                      const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  // 1. 參數打包 theta = [alpha, beta..., phi...]
  const arma::uword n_params = 1 + n_p + n_m;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta.subvec(n_p + 1, n_params - 1) = phi;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec curr_phi = curr_theta.subvec(n_p + 1, n_params - 1);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
      arma::vec psi_vec = get_psi_AD_EYsubX_logistic_rcpp(X.row(i).t(),
                                                          curr_alpha,
                                                          curr_beta,
                                                          curr_phi,
                                                          inc_vec);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec Xrow_i = X.row(i).t();
    arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());

    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    double prob_1 = eSI_i / (1.0 + eSI_i);

    // Psi
    arma::vec Psi_i = (prob_1 - phi) % inc_vec;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;                 // [新增] 累積加權平均
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) var_A (diff)
    var_A_accum -= inc_vec * inv_denom;

    // (C) var_B (var)
    double prob_0 = 1.0 / (1.0 + eSI_i);
    arma::vec term1 = (1.0 - phi) % inc_vec;
    arma::vec term2 = phi % inc_vec;

    arma::mat contribution = (term1 * term1.t() * prob_1 + term2 * term2.t() * prob_0) * inv_denom;

    var_B_accum += contribution;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;
  var_A_accum /= n;
  var_B_accum /= n;

  // 防呆與 ADvar 計算
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B_accum.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi
    Named("ADvar") = ADvar
  );
}

// Covariate shift, EY

arma::vec get_psi_ADCS_EY_logistic_rcpp(const arma::vec & Xrow_i,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const double & phi,
                                        const arma::vec & CS_beta) {

  // 1. 計算線性預測與指數項
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  if (eSI_i > 1e10) eSI_i = 1e10; // 截斷

  // 2. 計算 P(Y=1|X)
  double prob_1 = eSI_i / (1.0 + eSI_i);

  // 3. 計算 Covariate Shift 權重
  double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
  if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

  // 4. 計算 Psi
  // 公式: (P(Y=1) - phi) * w(X)
  return arma::vec({(prob_1 - phi) * eSI_CS_i});
}

// [[Rcpp::export]]
List SolveLagrange_ADCS_EY_logistic_rcpp(const arma::mat & X,
                                         const double & alpha,
                                         const arma::vec & beta,
                                         const double & phi,
                                         const arma::vec & CS_beta,
                                         const double & eta_initial,
                                         const size_t & iter_max,
                                         const double & step_rate,
                                         const size_t & step_max,
                                         const double & tol,
                                         const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      // eta 是 scalar
      return get_psi_ADCS_EY_logistic_rcpp(X.row(i).t(), alpha, beta, phi, CS_beta);
    };
    return EL_Core::compute_lagrange(n, arma::vec({e(0)}), get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(arma::vec({eta_initial}),
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADCS_EY_logistic_rcpp(const arma::mat & X,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & phi,
                                    const arma::vec & CS_beta,
                                    const double & eta,
                                    const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = 1; // Psi 維度

  // 1. 參數打包 theta = [alpha, beta..., phi, CS_beta...]
  const arma::uword n_params = 1 + n_p + 1 + n_p;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = phi;
  theta.subvec(n_p + 2, n_params - 1) = CS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_phi = curr_theta(n_p + 1);
    arma::vec curr_CS_beta = curr_theta.subvec(n_p + 2, n_params - 1);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADCS_EY_logistic_rcpp(X.row(i).t(),
                                                        curr_alpha,
                                                        curr_beta,
                                                        curr_phi,
                                                        curr_CS_beta);

      double denom = 1.0 + eta * psi_vec(0);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  double var_B_Ntor = 0.0;
  double var_B_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec Xrow_i = X.row(i).t();

    // Logistic Model
    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    // Covariate Shift
    double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

    // Psi = (P(1) - phi) * w
    double prob_1 = eSI_i / (1.0 + eSI_i);
    double Psi_val = (prob_1 - phi) * eSI_CS_i;

    // 權重計算
    double denom = 1.0 + eta * Psi_val;
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi(0) += Psi_val * inv_denom;            // [新增] 累積加權平均
    mean_psi_outer(0,0) += (Psi_val * Psi_val) * inv_denom;

    // (B) 準備自歸一化變異數計算
    double weighted_w = eSI_CS_i * inv_denom;

    // 分子計算
    double term_numerator = (1.0 - phi) * (1.0 - phi) * eSI_i + phi * phi;
    double term = term_numerator / (1.0 + eSI_i);

    var_B_Ntor += term * weighted_w;
    var_B_Dtor += weighted_w;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_B_Dtor < 1e-10) var_B_Dtor = 1e-10;
  double var_B_val = var_B_Ntor / var_B_Dtor;

  // 組合 ADvar
  arma::vec ADvar = {var_B_val};

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi, CS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi, CS_beta
    Named("ADvar") = ADvar
  );
}

// Covariate shift, EXsubY

arma::vec get_psi_ADCS_EXsubY_logistic_rcpp(const arma::vec & Xrow_i,
                                            const double & alpha,
                                            const arma::vec & beta,
                                            const arma::mat & phi,
                                            const arma::vec & CS_beta) {

  // 1. 計算 Logistic Probability
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  if (eSI_i > 1e10) eSI_i = 1e10; // 截斷

  double prob_1 = eSI_i / (1.0 + eSI_i);
  double prob_0 = 1.0 / (1.0 + eSI_i);

  // 2. 計算 Covariate Shift 權重
  double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
  if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

  // 3. 計算 Psi 分量 (加權)
  // Psi_1 = (X - phi_1) * P(1) * w
  arma::vec psi_1 = (Xrow_i - phi.row(0).t()) * prob_1 * eSI_CS_i;

  // Psi_0 = (X - phi_0) * P(0) * w
  arma::vec psi_0 = (Xrow_i - phi.row(1).t()) * prob_0 * eSI_CS_i;

  // 4. 組合 (長度 2p)
  return arma::join_cols(psi_1, psi_0);
}

// [[Rcpp::export]]
List SolveLagrange_ADCS_EXsubY_logistic_rcpp(const arma::mat & X,
                                             const double & alpha,
                                             const arma::vec & beta,
                                             const arma::mat & phi,
                                             const arma::vec & CS_beta,
                                             const arma::vec & eta_initial,
                                             const size_t & iter_max,
                                             const double & step_rate,
                                             const size_t & step_max,
                                             const double & tol,
                                             const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADCS_EXsubY_logistic_rcpp(X.row(i).t(), alpha, beta, phi, CS_beta);
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADCS_EXsubY_logistic_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const arma::mat & phi,
                                        const arma::vec & CS_beta,
                                        const arma::vec & eta,
                                        const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;

  // 1. 參數打包 theta = [alpha, beta..., phi..., CS_beta...]
  const arma::uword n_params = 1 + n_p + n_m + n_p;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta.subvec(n_p + 1, n_p + n_m) = arma::vectorise(phi.t());
  theta.subvec(n_p + 1 + n_m, n_params - 1) = CS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec phi_vec = curr_theta.subvec(n_p + 1, n_p + n_m);
    arma::mat curr_phi = arma::reshape(phi_vec, n_p, n_k).t();
    arma::vec curr_CS_beta = curr_theta.subvec(n_p + 1 + n_m, n_params - 1);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADCS_EXsubY_logistic_rcpp(X.row(i).t(),
                                                            curr_alpha,
                                                            curr_beta,
                                                            curr_phi,
                                                            curr_CS_beta);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);       // [新增] mean_psi
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec Xrow_i = X.row(i).t();
    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    double prob_1 = eSI_i / (1.0 + eSI_i);
    double prob_0 = 1.0 / (1.0 + eSI_i);

    double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

    // Psi
    arma::vec psi_1 = (Xrow_i - phi.row(0).t()) * prob_1 * eSI_CS_i;
    arma::vec psi_0 = (Xrow_i - phi.row(1).t()) * prob_0 * eSI_CS_i;
    arma::vec Psi_i = arma::join_cols(psi_1, psi_0);

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;                 // [新增] 累積加權平均
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 準備自歸一化變異數計算
    double weighted_w = eSI_CS_i * inv_denom;
    var_AB_Dtor += weighted_w;

    // (C) var_A (diff) - 對角線
    var_A_accum.subvec(0, n_p - 1) -= prob_1 * weighted_w;
    var_A_accum.subvec(n_p, 2 * n_p - 1) -= prob_0 * weighted_w;

    // (D) var_B (var) - Block Diagonal
    arma::vec diff_1 = (Xrow_i - phi.row(0).t());
    arma::vec diff_0 = (Xrow_i - phi.row(1).t());

    // Block 1 (Top-Left)
    var_B_accum.submat(0, 0, n_p - 1, n_p - 1)
      += (diff_1 * diff_1.t()) * prob_1 * weighted_w;

    // Block 2 (Bottom-Right)
    var_B_accum.submat(n_p, n_p, 2 * n_p - 1, 2 * n_p - 1)
      += (diff_0 * diff_0.t()) * prob_0 * weighted_w;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_AB_Dtor < 1e-10) var_AB_Dtor = 1e-10;
  var_A_accum /= var_AB_Dtor;
  var_B_accum /= var_AB_Dtor;

  // 防呆機制
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B_accum.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi, CS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,                  // [新增] 回傳 mean_psi
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi, CS_beta
    Named("ADvar") = ADvar
  );
}

// Covariate shift, EYsubX

arma::vec get_psi_ADCS_EYsubX_logistic_rcpp(const arma::vec & Xrow_i,
                                            const double & alpha,
                                            const arma::vec & beta,
                                            const arma::vec & phi,
                                            const arma::vec & CS_beta,
                                            const arma::vec & inclusion) {

  // 1. 計算 Logistic Probability
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  if (eSI_i > 1e10) eSI_i = 1e10; // 截斷

  double prob_1 = eSI_i / (1.0 + eSI_i);

  // 2. 計算 Covariate Shift 權重
  double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
  if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

  // 3. 計算 Psi (未遮罩)
  arma::vec psi_raw = (prob_1 - phi) * eSI_CS_i;

  // 4. 應用遮罩
  return psi_raw % inclusion;
}

// [[Rcpp::export]]
List SolveLagrange_ADCS_EYsubX_logistic_rcpp(const arma::mat & X,
                                             const double & alpha,
                                             const arma::vec & beta,
                                             const arma::vec & phi,
                                             const arma::vec & CS_beta,
                                             const arma::umat & inclusion,
                                             const arma::vec & eta_initial,
                                             const size_t & iter_max,
                                             const double & step_rate,
                                             const size_t & step_max,
                                             const double & tol,
                                             const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADCS_EYsubX_logistic_rcpp(X.row(i).t(),
                                               alpha,
                                               beta,
                                               phi,
                                               CS_beta,
                                               arma::conv_to<arma::vec>::from(inclusion.row(i).t()));
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADCS_EYsubX_logistic_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const arma::vec & phi,
                                        const arma::vec & CS_beta,
                                        const arma::umat & inclusion,
                                        const arma::vec & eta,
                                        const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  // 1. 參數打包 theta = [alpha, beta..., phi..., CS_beta...]
  const arma::uword n_params = 1 + n_p + n_m + n_p;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta.subvec(n_p + 1, n_p + n_m) = phi;
  theta.subvec(n_p + 1 + n_m, n_params - 1) = CS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec curr_phi = curr_theta.subvec(n_p + 1, n_p + n_m);
    arma::vec curr_CS_beta = curr_theta.subvec(n_p + 1 + n_m, n_params - 1);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
      arma::vec psi_vec = get_psi_ADCS_EYsubX_logistic_rcpp(X.row(i).t(),
                                                            curr_alpha,
                                                            curr_beta,
                                                            curr_phi,
                                                            curr_CS_beta,
                                                            inc_vec);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec Xrow_i = X.row(i).t();
    arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());

    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    double eSI_CS_i = std::exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) eSI_CS_i = 1e10;

    double prob_1 = eSI_i / (1.0 + eSI_i);

    // Psi
    arma::vec Psi_i = (prob_1 - phi) % inc_vec * eSI_CS_i;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 準備自歸一化變異數計算
    double weighted_w = eSI_CS_i * inv_denom;
    var_AB_Dtor += weighted_w;

    // (C) var_A (diff)
    var_A_accum -= inc_vec * weighted_w;

    // (D) var_B (var)
    double prob_0 = 1.0 / (1.0 + eSI_i);
    arma::vec term1 = (1.0 - phi) % inc_vec;
    arma::vec term2 = phi % inc_vec;

    arma::mat contribution = (term1 * term1.t() * prob_1 + term2 * term2.t() * prob_0) * weighted_w;

    var_B_accum += contribution;
  }

  mean_psi /= n;
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_AB_Dtor < 1e-10) var_AB_Dtor = 1e-10;
  var_A_accum /= var_AB_Dtor;
  var_B_accum /= var_AB_Dtor;

  // 防呆機制
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B_accum.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi, CS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi, CS_beta
    Named("ADvar") = ADvar
  );
}

// Prior probability shift, EX

arma::vec get_psi_ADPPS_EX_logistic_rcpp(const arma::vec & Xrow_i,
                                         const double & alpha,
                                         const arma::vec & beta,
                                         const arma::vec & phi,
                                         const double & PPS_beta) {

  // 1. 計算 eSI
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  // 截斷防護
  if (eSI_i > 1e10) eSI_i = 1e10;

  // 2. 計算 PPS 權重
  double ePPS_beta = std::exp(PPS_beta);
  double weight = (1.0 + eSI_i * ePPS_beta) / (1.0 + eSI_i);

  // 3. 計算 Psi
  // (X - phi) * w
  return (Xrow_i - phi) * weight;
}

// [[Rcpp::export]]
List SolveLagrange_ADPPS_EX_logistic_rcpp(const arma::mat & X,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const arma::vec & phi,
                                          const double & PPS_beta,
                                          const arma::vec & eta_initial,
                                          const size_t & iter_max,
                                          const double & step_rate,
                                          const size_t & step_max,
                                          const double & tol,
                                          const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADPPS_EX_logistic_rcpp(X.row(i).t(), alpha, beta, phi, PPS_beta);
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADPPS_EX_logistic_rcpp(const arma::mat & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const arma::vec & phi,
                                     const double & PPS_beta,
                                     const arma::vec & eta,
                                     const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem; // 應等於 n_p

  // 1. 參數打包 theta = [alpha, beta..., phi..., PPS_beta]
  const arma::uword n_params = 1 + n_p + n_m + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta.subvec(n_p + 1, n_p + n_m) = phi; // 加入 phi
  theta(n_p + 1 + n_m) = PPS_beta;        // PPS_beta 移至最後

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec curr_phi = curr_theta.subvec(n_p + 1, n_p + n_m);
    double curr_PPS_beta = curr_theta(n_p + 1 + n_m);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADPPS_EX_logistic_rcpp(X.row(i).t(),
                                                         curr_alpha,
                                                         curr_beta,
                                                         curr_phi,
                                                         curr_PPS_beta);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }

    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_B (使用原始參數)
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);
  double var_B_Dtor = 0.0;

  double ePPS_beta = std::exp(PPS_beta); // 預先計算常數

  for (arma::uword i = 0; i < n; ++i) {
    // 準備數據
    arma::vec Xrow_i = X.row(i).t();

    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    // PPS 權重
    double weight = (1.0 + eSI_i * ePPS_beta) / (1.0 + eSI_i);

    // Psi
    arma::vec Psi_i = (Xrow_i - phi) * weight;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) 準備自歸一化變異數計算
    double weighted_w = weight * inv_denom;

    var_B_Dtor += weighted_w;

    // 累積 var_B 分子
    arma::vec diff = (Xrow_i - phi);
    var_B_accum += (diff * diff.t()) * weighted_w;
  }

  mean_psi /= n;
  mean_psi_outer /= n;

  // 自歸一化 normalization
  if (var_B_Dtor < 1e-10) var_B_Dtor = 1e-10;
  var_B_accum /= var_B_Dtor;

  // 組合 ADvar
  arma::mat ADvar = var_B_accum;

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi, PPS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi, PPS_beta
    Named("ADvar") = ADvar
  );
}

// Prior probability shift, EY

arma::vec get_psi_ADPPS_EY_logistic_rcpp(const arma::vec & Xrow_i,
                                         const double & alpha,
                                         const arma::vec & beta,
                                         const double & phi,
                                         const double & PPS_beta) {

  // 1. 預先計算常數 ePPS
  double ePPS_beta = std::exp(PPS_beta);

  // 2. 計算 Logistic Probability components
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  if (eSI_i > 1e10) eSI_i = 1e10; // 截斷

  double prob_1 = eSI_i / (1.0 + eSI_i);
  double prob_0 = 1.0 - prob_1;

  // 3. 計算 Psi
  double psi_val = (1.0 - phi) * ePPS_beta * prob_1 - phi * prob_0;

  return arma::vec({psi_val});
}

// [[Rcpp::export]]
List SolveLagrange_ADPPS_EY_logistic_rcpp(const arma::mat & X,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const double & phi,
                                          const double & PPS_beta,
                                          const double & eta_initial,
                                          const size_t & iter_max,
                                          const double & step_rate,
                                          const size_t & step_max,
                                          const double & tol,
                                          const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADPPS_EY_logistic_rcpp(X.row(i).t(), alpha, beta, phi, PPS_beta);
    };
    return EL_Core::compute_lagrange(n, arma::vec({e(0)}), get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(arma::vec({eta_initial}),
                               obj_func,
                               iter_max,
                               step_rate,
                               step_max,
                               tol,
                               eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADPPS_EY_logistic_rcpp(const arma::mat & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const double & phi,
                                     const double & PPS_beta,
                                     const double & eta,
                                     const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = 1; // Psi 維度

  // 1. 參數打包 theta = [alpha, beta..., phi, PPS_beta]
  const arma::uword n_params = 1 + n_p + 1 + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta(n_p + 1) = phi;
  theta(n_p + 2) = PPS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    double curr_phi = curr_theta(n_p + 1);
    double curr_PPS_beta = curr_theta(n_p + 2);
    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec psi_vec = get_psi_ADPPS_EY_logistic_rcpp(X.row(i).t(),
                                                         curr_alpha,
                                                         curr_beta,
                                                         curr_phi,
                                                         curr_PPS_beta);

      double denom = 1.0 + eta * psi_vec(0);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }

    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_B
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  double var_B_accum = 0.0;

  // 預先計算 PPS 常數
  double ePPS_beta = std::exp(PPS_beta);

  for (arma::uword i = 0; i < n; ++i) {
    // 準備數據
    arma::vec Xrow_i = X.row(i).t();
    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    double prob_1 = eSI_i / (1.0 + eSI_i);
    double prob_0 = 1.0 - prob_1;
    double Psi_val = (1.0 - phi) * ePPS_beta * prob_1 - phi * prob_0;

    // 權重計算
    double denom = 1.0 + eta * Psi_val;
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi_outer & mean_psi
    mean_psi(0) += Psi_val * inv_denom;
    mean_psi_outer(0,0) += (Psi_val * Psi_val) * inv_denom;

    // (B) var_B 計算
    double term = (1.0 - phi) * (1.0 - phi) * ePPS_beta * prob_1 + phi * phi * prob_0;

    var_B_accum += term * inv_denom;
  }

  mean_psi /= n;        // [新增] 取平均
  mean_psi_outer /= n;
  var_B_accum /= n;

  // 組合 ADvar
  arma::vec ADvar = {var_B_accum};

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi, PPS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi, PPS_beta
    Named("ADvar") = ADvar
  );
}

// Prior probability shift, EYsubX

arma::vec get_psi_ADPPS_EYsubX_logistic_rcpp(const arma::vec & Xrow_i,
                                             const double & alpha,
                                             const arma::vec & beta,
                                             const arma::vec & phi,
                                             const double & PPS_beta,
                                             const arma::vec & inclusion) {

  // 1. 計算 PPS 常數
  double ePPS_beta = std::exp(PPS_beta);

  // 2. 計算 Logistic Probabilities
  double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
  if (eSI_i > 1e10) eSI_i = 1e10; // 截斷

  double prob_1 = eSI_i / (1.0 + eSI_i);
  double prob_0 = 1.0 - prob_1;

  // 3. 計算 Psi (未遮罩)
  arma::vec term1 = (1.0 - phi) * (ePPS_beta * prob_1);
  arma::vec term2 = phi * prob_0;
  arma::vec psi_raw = term1 - term2;

  // 4. 應用遮罩
  return psi_raw % inclusion;
}

// [[Rcpp::export]]
List SolveLagrange_ADPPS_EYsubX_logistic_rcpp(const arma::mat & X,
                                              const double & alpha,
                                              const arma::vec & beta,
                                              const arma::vec & phi,
                                              const double & PPS_beta,
                                              const arma::umat & inclusion,
                                              const arma::vec & eta_initial,
                                              const size_t & iter_max,
                                              const double & step_rate,
                                              const size_t & step_max,
                                              const double & tol,
                                              const double & eps_inv) {

  const arma::uword n = X.n_rows;

  // 定義目標函數呼叫器
  auto obj_func = [&](const arma::vec & e) {
    auto get_psi_lambda = [&](size_t i) -> arma::vec {
      return get_psi_ADPPS_EYsubX_logistic_rcpp(X.row(i).t(),
                                                alpha,
                                                beta,
                                                phi,
                                                PPS_beta,
                                                arma::conv_to<arma::vec>::from(inclusion.row(i).t()));
    };
    return EL_Core::compute_lagrange(n, e, get_psi_lambda);
  };

  // 執行牛頓法求解
  return EL_Core::solve_newton(eta_initial, obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}

// [[Rcpp::export]]
List grad_psi_ADPPS_EYsubX_logistic_rcpp(const arma::mat & X,
                                         const double & alpha,
                                         const arma::vec & beta,
                                         const arma::vec & phi,
                                         const double & PPS_beta,
                                         const arma::umat & inclusion,
                                         const arma::vec & eta,
                                         const double & diff_size) {

  const arma::uword n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  // 1. 參數打包 theta = [alpha, beta..., phi..., PPS_beta]
  const arma::uword n_params = 1 + n_p + n_m + 1;
  arma::vec theta = arma::zeros(n_params);

  theta(0) = alpha;
  theta.subvec(1, n_p) = beta;
  theta.subvec(n_p + 1, n_p + n_m) = phi;
  theta(n_p + 1 + n_m) = PPS_beta;

  // 2. 輔助函數：計算 "加權" 平均 Psi
  auto compute_weighted_avg_psi = [&](const arma::vec & curr_theta) -> arma::vec {

    double curr_alpha = curr_theta(0);
    arma::vec curr_beta = curr_theta.subvec(1, n_p);
    arma::vec curr_phi = curr_theta.subvec(n_p + 1, n_p + n_m);
    double curr_PPS_beta = curr_theta(n_p + 1 + n_m);

    arma::vec sum_val = arma::zeros(n_m);

    for (arma::uword i = 0; i < n; ++i) {

      arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
      arma::vec psi_vec = get_psi_ADPPS_EYsubX_logistic_rcpp(X.row(i).t(),
                                                             curr_alpha,
                                                             curr_beta,
                                                             curr_phi,
                                                             curr_PPS_beta,
                                                             inc_vec);

      double denom = 1.0 + arma::dot(eta, psi_vec);
      if (denom < 1e-10) denom = 1e-10;

      sum_val += psi_vec / denom;
    }
    return sum_val / n;
  };

  // 3. 主迴圈：計算 mean_psi, mean_psi_outer, var_A, var_B
  arma::vec mean_psi = arma::zeros(n_m);
  arma::mat mean_psi_outer = arma::zeros(n_m, n_m);
  arma::vec var_A_accum = arma::zeros(n_m);
  arma::mat var_B_accum = arma::zeros(n_m, n_m);

  // 預計算 PPS 參數
  double ePPS_beta = std::exp(PPS_beta);

  for (arma::uword i = 0; i < n; ++i) {
    // 準備數據
    arma::vec Xrow_i = X.row(i).t();
    arma::vec inc_vec = arma::conv_to<arma::vec>::from(inclusion.row(i).t());

    double eSI_i = std::exp(alpha + arma::dot(Xrow_i, beta));
    if (eSI_i > 1e10) eSI_i = 1e10;

    double prob_1 = eSI_i / (1.0 + eSI_i);
    double prob_0 = 1.0 - prob_1;

    // Psi
    arma::vec term1 = (1.0 - phi) * (ePPS_beta * prob_1);
    arma::vec term2 = phi * prob_0;
    arma::vec Psi_i = (term1 - term2) % inc_vec;

    // 權重計算
    double denom = 1.0 + arma::dot(eta, Psi_i);
    if (denom < 1e-10) denom = 1e-10;
    double inv_denom = 1.0 / denom;

    // (A) mean_psi & mean_psi_outer
    mean_psi += Psi_i * inv_denom;
    mean_psi_outer += (Psi_i * Psi_i.t()) * inv_denom;

    // (B) var_A (diff) - 對角線
    var_A_accum -= inc_vec * inv_denom;

    // (C) var_B (var)
    arma::vec vec_1_phi = (1.0 - phi) % inc_vec;
    arma::vec vec_phi   = phi % inc_vec;

    arma::mat contrib = (vec_1_phi * vec_1_phi.t()) * (prob_1 * ePPS_beta) +
      (vec_phi * vec_phi.t()) * prob_0;

    var_B_accum += contrib * inv_denom;
  }

  mean_psi /= n;
  mean_psi_outer /= n;
  var_A_accum /= n;
  var_B_accum /= n;

  // 防呆機制與 ADvar 計算
  double epsilon = 1e-10;
  for(arma::uword k = 0; k < n_m; ++k) {
    if (std::abs(var_A_accum(k)) < epsilon) {
      var_A_accum(k) = (var_A_accum(k) < 0.0) ? -epsilon : epsilon;
    }
  }

  arma::vec var_A_inv = 1.0 / var_A_accum;
  arma::mat ADvar = (var_B_accum.each_col() % var_A_inv).each_row() % var_A_inv.t();

  // 4. 數值微分計算 Gradient (對 alpha, beta, phi, PPS_beta)
  arma::mat gradient = arma::zeros(n_m, n_params);

  for (arma::uword j = 0; j < n_params; ++j) {
    arma::vec theta_plus = theta;
    arma::vec theta_minus = theta;

    theta_plus(j) += diff_size;
    theta_minus(j) -= diff_size;

    arma::vec f_plus = compute_weighted_avg_psi(theta_plus);
    arma::vec f_minus = compute_weighted_avg_psi(theta_minus);

    // Central Difference
    gradient.col(j) = (f_plus - f_minus) / (2.0 * diff_size);
  }

  return List::create(
    Named("mean_psi") = mean_psi,
    Named("mean_psi_outer") = mean_psi_outer,
    Named("mean_psi_gradient") = gradient, // Order: alpha, beta, phi, PPS_beta
    Named("ADvar") = ADvar
  );
}

// =========================================================================
// 4. Gamma Regression Models
// =========================================================================

// No shift, EY

// No shift, EXsubY

// No shift, EYsubX

// Covariate shift, EY

// Covariate shift, EXsubY

// Covariate shift, EYsubX

// Prior probability shift, EX

// Prior probability shift, EY

// Prior probability shift, EXsubY

// Prior probability shift, EYsubX

















//===========================================================================


// [[Rcpp::export]]
List AD_EY_Gamma_Lagrange_rcpp(const arma::mat & X,
                               double alpha,
                               const arma::vec & beta,
                               double nu,
                               double phi,
                               double eta) {
  const arma::uword n = X.n_rows;
  auto get_psi = [&](size_t i) -> arma::vec {
    double eSI_i = std::exp(alpha + arma::dot(X.row(i), beta));
    return arma::vec({ eSI_i - phi });
  };
  return EL_Core::compute_lagrange(n, arma::vec({eta}), get_psi);
}

// [[Rcpp::export]]
List AD_EY_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                    double alpha,
                                    const arma::vec & beta,
                                    double nu,
                                    double phi,
                                    double eta_initial,
                                    size_t iter_max,
                                    double step_rate,
                                    size_t step_max,
                                    double tol,
                                    double eps_inv) {
  auto obj_func = [&](const arma::vec& e) {
    return AD_EY_Gamma_Lagrange_rcpp(X, alpha, beta, nu, phi, e(0));
  };
  return EL_Core::solve_newton(arma::vec({eta_initial}), obj_func, iter_max, step_rate, step_max, tol, eps_inv);
}


// [[Rcpp::export]]
List AD_EXsubY_Gamma_Lagrange_rcpp(const arma::mat & X,
                                   const double & alpha,
                                   const arma::vec & beta,
                                   const double & nu,
                                   const arma::mat & phi,
                                   const arma::mat & y_pts,
                                   const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double lambda_i_inv = 1 / lambda_i;
    arma::mat Psi_i_mat(n_k, n_p);

    for (size_t k = 0; k < n_k; ++k) {

      Psi_i_mat.row(k) = (Xrow_i.t() - phi.row(k)) *
        (R::pgamma(y_pts(k, 1), nu, lambda_i_inv, TRUE, FALSE) -
        R::pgamma(y_pts(k, 0), nu, lambda_i_inv, TRUE, FALSE));
    }

    const arma::vec Psi_i = arma::reshape(
      Psi_i_mat.t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !gradient.is_finite() ||
      !hessian.is_finite()) {

    value = -10000.0 * n_n;
    gradient = arma::vec(n_m);
    hessian = arma::mat(n_m, n_m);
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List AD_EXsubY_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const double & nu,
                                        const arma::mat & phi,
                                        const arma::mat & y_pts,
                                        const arma::vec & eta_initial,
                                        const size_t & iter_max,
                                        const double & step_rate,
                                        const size_t & step_max,
                                        const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = AD_EXsubY_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, y_pts, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const arma::vec gradient = step["gradient"];
    const arma::mat hessian = step["hessian"];
    const double ind_NT_GD = arma::rcond(hessian);

    arma::vec direction_step(n_m);

    if (ind_NT_GD > tol) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = AD_EXsubY_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, y_pts, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EXsubY_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, y_pts, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// [[Rcpp::export]]
List AD_EXsubY_Gamma_rcpp(const arma::mat & X,
                          const double & alpha,
                          const arma::vec & beta,
                          const double & nu,
                          const arma::mat & phi,
                          const arma::mat & y_pts) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 2);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double lambda_i_inv = 1 / lambda_i;
    arma::mat Psi_i_mat(n_k, n_p);

    for (size_t k = 0; k < n_k; ++k) {

      Psi_i_mat.row(k) = (Xrow_i.t() - phi.row(k)) *
        (R::pgamma(y_pts(k, 1), nu, lambda_i_inv, TRUE, FALSE) -
        R::pgamma(y_pts(k, 0), nu, lambda_i_inv, TRUE, FALSE));
    }

    const arma::vec Psi_i = arma::reshape(
      Psi_i_mat.t(), n_m, 1);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
  }










  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List AD_EYsubX_Gamma_Lagrange_rcpp(const arma::mat & X,
                                   const double & alpha,
                                   const arma::vec & beta,
                                   const double & nu,
                                   const arma::vec & phi,
                                   const arma::umat & inclusion,
                                   const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += eSI_i - phi(m);
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !gradient.is_finite() ||
      !hessian.is_finite()) {

    value = -10000.0 * n_n;
    gradient = arma::vec(n_m);
    hessian = arma::mat(n_m, n_m);
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List AD_EYsubX_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const double & nu,
                                        const arma::vec & phi,
                                        const arma::umat & inclusion,
                                        const arma::vec & eta_initial,
                                        const size_t & iter_max,
                                        const double & step_rate,
                                        const size_t & step_max,
                                        const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = AD_EYsubX_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, inclusion, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const arma::vec gradient = step["gradient"];
    const arma::mat hessian = step["hessian"];
    const double ind_NT_GD = arma::rcond(hessian);

    arma::vec direction_step(n_m);

    if (ind_NT_GD > tol) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = AD_EYsubX_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EYsubX_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, inclusion, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// Prior probability shift, Gamma regression model

// [[Rcpp::export]]
List ADPPS_EX_Gamma_Lagrange_rcpp(const arma::mat & X,
                                  const double & alpha,
                                  const arma::vec & beta,
                                  const double & nu,
                                  const arma::vec & phi,
                                  const double & PPS_beta,
                                  const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  arma::vec gradient(n_p);
  arma::mat hessian(n_p, n_p);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const arma::vec Psi_i = (Xrow_i - phi) *
      pow(lambda_i / (lambda_i - PPS_beta), nu);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !gradient.is_finite() ||
      !hessian.is_finite()) {

    value = -10000.0 * n_n;
    gradient = arma::vec(n_p);
    hessian = arma::mat(n_p, n_p);
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List ADPPS_EX_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & nu,
                                       const arma::vec & phi,
                                       const double & PPS_beta,
                                       const arma::vec & eta_initial,
                                       const size_t & iter_max,
                                       const double & step_rate,
                                       const size_t & step_max,
                                       const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADPPS_EX_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, PPS_beta, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const arma::vec gradient = step["gradient"];
    const arma::mat hessian = step["hessian"];
    const double ind_NT_GD = arma::rcond(hessian);

    arma::vec direction_step(n_m);

    if (ind_NT_GD > tol) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = ADPPS_EX_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, PPS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EX_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, PPS_beta, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// [[Rcpp::export]]
List ADPPS_EY_Gamma_Lagrange_rcpp(const arma::mat & X,
                                  const double & alpha,
                                  const arma::vec & beta,
                                  const double & nu,
                                  const double & phi,
                                  const double & PPS_beta,
                                  const double & eta) {

  const arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double Psi_i = (nu / (lambda_i - PPS_beta) - phi) *
      pow(lambda_i / (lambda_i - PPS_beta), nu);

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !std::isfinite(gradient) ||
      !std::isfinite(hessian)) {

    value = -10000.0 * n_n;
    gradient = 0.0;
    hessian = 0.0;
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List ADPPS_EY_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & nu,
                                       const double & phi,
                                       const double & PPS_beta,
                                       const double & eta_initial,
                                       const size_t & iter_max,
                                       const double & step_rate,
                                       const size_t & step_max,
                                       const double & tol) {

  double eta = eta_initial;
  List step = ADPPS_EY_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, PPS_beta, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const double gradient = step["gradient"];
    const double hessian = step["hessian"];

    double direction_step;

    if (hessian != 0) {

      direction_step = gradient / hessian;

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    double eta_new = eta - direction_step;
    List step_new = ADPPS_EY_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, PPS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EY_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, PPS_beta, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// [[Rcpp::export]]
List ADPPS_EXsubY_Gamma_Lagrange_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & nu,
                                      const arma::mat & phi,
                                      const double & PPS_beta,
                                      const arma::mat & y_pts,
                                      const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double lambda_i_PPS_beta = lambda_i - PPS_beta;
    const double lambda_i_PPS_beta_inv = 1 / lambda_i_PPS_beta;
    const double ratio_i = pow(lambda_i / lambda_i_PPS_beta, nu);
    arma::mat Psi_i_mat(n_k, n_p);

    for (size_t k = 0; k < n_k; ++k) {

      Psi_i_mat.row(k) = (Xrow_i.t() - phi.row(k)) * ratio_i *
        (R::pgamma(y_pts(k, 1), nu, lambda_i_PPS_beta_inv, TRUE, FALSE) -
        R::pgamma(y_pts(k, 0), nu, lambda_i_PPS_beta_inv, TRUE, FALSE));
    }

    const arma::vec Psi_i = arma::reshape(
      Psi_i_mat.t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !gradient.is_finite() ||
      !hessian.is_finite()) {

    value = -10000.0 * n_n;
    gradient = arma::vec(n_m);
    hessian = arma::mat(n_m, n_m);
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List ADPPS_EXsubY_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const double & nu,
                                           const arma::mat & phi,
                                           const double & PPS_beta,
                                           const arma::mat & y_pts,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADPPS_EXsubY_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, PPS_beta, y_pts, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const arma::vec gradient = step["gradient"];
    const arma::mat hessian = step["hessian"];
    const double ind_NT_GD = arma::rcond(hessian);

    arma::vec direction_step(n_m);

    if (ind_NT_GD > tol) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = ADPPS_EXsubY_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, PPS_beta, y_pts, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EXsubY_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, PPS_beta, y_pts, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// [[Rcpp::export]]
List ADPPS_EYsubX_Gamma_Lagrange_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & nu,
                                      const arma::vec & phi,
                                      const double & PPS_beta,
                                      const arma::umat & inclusion,
                                      const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += (nu / (lambda_i - PPS_beta) - phi(m)) *
          pow(lambda_i / (lambda_i - PPS_beta), nu);
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !gradient.is_finite() ||
      !hessian.is_finite()) {

    value = -10000.0 * n_n;
    gradient = arma::vec(n_m);
    hessian = arma::mat(n_m, n_m);
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List ADPPS_EYsubX_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const double & nu,
                                           const arma::vec & phi,
                                           const double & PPS_beta,
                                           const arma::umat & inclusion,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADPPS_EYsubX_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, PPS_beta, inclusion, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const arma::vec gradient = step["gradient"];
    const arma::mat hessian = step["hessian"];
    const double ind_NT_GD = arma::rcond(hessian);

    arma::vec direction_step(n_m);

    if (ind_NT_GD > tol) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = ADPPS_EYsubX_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, PPS_beta, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EYsubX_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, PPS_beta, inclusion, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// Covariate shift, Gamma regression model

// [[Rcpp::export]]
List ADCS_EY_Gamma_Lagrange_rcpp(const arma::mat & X,
                                 const double & alpha,
                                 const arma::vec & beta,
                                 const double & nu,
                                 const double & phi,
                                 const arma::vec & CS_beta,
                                 const double & eta) {

  const arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    const double Psi_i = eSI_CS_i * (eSI_i - phi);

    const double denominator = 1.0 + eta * Psi_i;
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -pow(Psi_i, 2) / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !std::isfinite(gradient) ||
      !std::isfinite(hessian)) {

    value = -10000.0 * n_n;
    gradient = 0.0;
    hessian = 0.0;
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List ADCS_EY_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & nu,
                                      const double & phi,
                                      const arma::vec & CS_beta,
                                      const double & eta_initial,
                                      const size_t & iter_max,
                                      const double & step_rate,
                                      const size_t & step_max,
                                      const double & tol) {

  double eta = eta_initial;
  List step = ADCS_EY_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, CS_beta, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const double gradient = step["gradient"];
    const double hessian = step["hessian"];

    double direction_step;

    if (hessian != 0) {

      direction_step = gradient / hessian;

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    double eta_new = eta - direction_step;
    List step_new = ADCS_EY_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, CS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EY_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, CS_beta, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// [[Rcpp::export]]
List ADCS_EXsubY_Gamma_Lagrange_rcpp(const arma::mat & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const double & nu,
                                     const arma::mat & phi,
                                     const arma::vec & CS_beta,
                                     const arma::mat & y_pts,
                                     const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double lambda_i = nu / eSI_i;
    const double lambda_i_inv = 1 / lambda_i;
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    arma::mat Psi_i_mat(n_k, n_p);

    for (size_t k = 0; k < n_k; ++k) {

      Psi_i_mat.row(k) = (Xrow_i.t() - phi.row(k)) * eSI_CS_i *
        (R::pgamma(y_pts(k, 1), nu, lambda_i_inv, TRUE, FALSE) -
        R::pgamma(y_pts(k, 0), nu, lambda_i_inv, TRUE, FALSE));
    }

    const arma::vec Psi_i = arma::reshape(
      Psi_i_mat.t(), n_m, 1);

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !gradient.is_finite() ||
      !hessian.is_finite()) {

    value = -10000.0 * n_n;
    gradient = arma::vec(n_m);
    hessian = arma::mat(n_m, n_m);
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List ADCS_EXsubY_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const double & nu,
                                          const arma::mat & phi,
                                          const arma::vec & CS_beta,
                                          const arma::mat & y_pts,
                                          const arma::vec & eta_initial,
                                          const size_t & iter_max,
                                          const double & step_rate,
                                          const size_t & step_max,
                                          const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADCS_EXsubY_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, CS_beta, y_pts, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const arma::vec gradient = step["gradient"];
    const arma::mat hessian = step["hessian"];
    const double ind_NT_GD = arma::rcond(hessian);

    arma::vec direction_step(n_m);

    if (ind_NT_GD > tol) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = ADCS_EXsubY_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, CS_beta, y_pts, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EXsubY_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, CS_beta, y_pts, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}

// [[Rcpp::export]]
List ADCS_EYsubX_Gamma_Lagrange_rcpp(const arma::mat  & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const double & nu,
                                     const arma::vec & phi,
                                     const arma::vec & CS_beta,
                                     const arma::umat & inclusion,
                                     const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += eSI_CS_i * (eSI_i - phi(m));
      }
    }

    const double denominator = 1.0 + arma::dot(eta, Psi_i);
    if (denominator > 0.0) {

      value += log(denominator);
      gradient += Psi_i / denominator;
      hessian += -Psi_i * Psi_i.t() / pow(denominator, 2);

    } else {

      value += -10000.0;
    }
  }

  if (!std::isfinite(value) ||
      !gradient.is_finite() ||
      !hessian.is_finite()) {

    value = -10000.0 * n_n;
    gradient = arma::vec(n_m);
    hessian = arma::mat(n_m, n_m);
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
List ADCS_EYsubX_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const double & nu,
                                          const arma::vec & phi,
                                          const arma::vec & CS_beta,
                                          const arma::umat & inclusion,
                                          const arma::vec & eta_initial,
                                          const size_t & iter_max,
                                          const double & step_rate,
                                          const size_t & step_max,
                                          const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADCS_EYsubX_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, CS_beta, inclusion, eta);

  for (size_t k = 0; k < iter_max; k++) {

    const double step_value = step["value"];
    const arma::vec gradient = step["gradient"];
    const arma::mat hessian = step["hessian"];
    const double ind_NT_GD = arma::rcond(hessian);

    arma::vec direction_step(n_m);

    if (ind_NT_GD > tol) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    double step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = ADCS_EYsubX_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, CS_beta, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EYsubX_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, CS_beta, inclusion, eta_new);

      } else {

        break;
      }
    }

    const double step_new_value = step_new["value"];

    if (step_new_value > step_value + tol) {

      eta = eta_new;
      step = step_new;

    } else {

      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("value") = step["value"],
    Rcpp::Named("gradient") = step["gradient"],
    Rcpp::Named("hessian") = step["hessian"]
  );
}


