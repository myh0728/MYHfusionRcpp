#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include "head_algebra.h"

// [[Rcpp::export]]
double lL_normal_rcpp(const arma::mat & X,
                      const arma::vec & Y,
                      const double & alpha,
                      const arma::vec & beta,
                      const double & sigma,
                      Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n_n = X.n_rows;

  // 1. 向量化計算殘差 (Residuals)
  // R: res <- Y - alpha - X %*% beta
  arma::vec res = Y - (alpha + X * beta);

  // 2. 準備常數項
  double sigma2 = pow(sigma, 2);
  double log_const = log(2.0 * arma::datum::pi * sigma2);
  double lL = 0.0;

  // 3. 根據是否有權重進行分支計算
  if (w.isNotNull()) {

    // --- 有權重 (Weighted) ---
    arma::vec weights = as<arma::vec>(w);

    // 防呆：檢查維度
    if (weights.n_elem != n_n) stop("Weights length differs from data rows.");

    // Weighted SSE = sum(w * res^2)
    // 使用 dot product 加速計算
    double w_sse = dot(weights, square(res));
    double sum_w = sum(weights);

    // Log-Likelihood
    lL = -0.5 * (w_sse / sigma2 + sum_w * log_const);

  } else {

    // --- 無權重 (Unweighted) ---

    // SSE = sum(res^2)
    double sse = sum(square(res));

    // Log-Likelihood
    lL = -0.5 * (sse / sigma2 + (double)n_n * log_const);
  }

  // 4. 數值穩定性檢查
  if (!std::isfinite(lL)) {
    return -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_normal_rcpp(const arma::mat & X,
                         const arma::vec & Y,
                         const double & alpha,
                         const arma::vec & beta,
                         const double & sigma,
                         Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

  // 1. 準備常數與殘差
  // res = Y - alpha - X * beta
  arma::vec res = Y - (alpha + X * beta);

  // --- 優化：改用直接乘法 ---
  double sigma2 = sigma * sigma;
  double sigma3 = sigma2 * sigma;
  double sigma4 = sigma2 * sigma2; // sigma^4 = (sigma^2)^2

  // 準備輸出容器
  // Index: 0 -> alpha, 1..p -> beta, p+1 -> sigma
  arma::vec gradient(p + 2, arma::fill::zeros);
  arma::mat hessian(p + 2, p + 2, arma::fill::zeros);

  // 2. 計算統計量 (根據有無權重)
  double sum_w = 0.0;     // sum(w) or n
  double sum_wr2 = 0.0;   // sum(w * res^2)
  arma::vec w_res;        // w * res
  arma::vec sum_wX;       // sum(w * X)
  arma::mat Xt_WX;        // X' * W * X

  if (w.isNotNull()) {

    // --- 有權重 (Weighted) ---
    arma::vec weights = as<arma::vec>(w);

    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    sum_w = arma::sum(weights);
    w_res = weights % res;            // element-wise multiply
    sum_wr2 = arma::dot(w_res, res);  // dot product

    // 計算 Hessian 需要的矩陣項
    // H_ab: sum(w * X) => (w^T * X)^T
    sum_wX = (weights.t() * X).t();

    // H_bb: X' * diag(w) * X
    // 優化技巧：避免建立巨大的對角矩陣
    arma::mat X_weighted = X;
    X_weighted.each_col() %= weights; // column-wise multiplication
    Xt_WX = X.t() * X_weighted;

  } else {

    // --- 無權重 (Unweighted) ---
    sum_w = (double)n;
    w_res = res;
    sum_wr2 = arma::dot(res, res);

    sum_wX = arma::sum(X, 0).t(); // sum columns, transform to col vector
    Xt_WX = X.t() * X;
  }

  // 3. 組合 Gradient
  // dL/d_alpha = sum(w * res) / sigma^2
  // dL/d_beta  = X' * (w * res) / sigma^2
  // dL/d_sigma = -sum(w)/sigma + sum(w * res^2)/sigma^3

  gradient(0) = arma::sum(w_res) / sigma2;                // alpha
  gradient.subvec(1, p) = (X.t() * w_res) / sigma2;       // beta
  gradient(p + 1) = -sum_w / sigma + sum_wr2 / sigma3;    // sigma

  // 4. 組合 Hessian
  // 結構:
  // [aa  ab  as]
  // [ba  bb  bs]
  // [sa  sb  ss]

  // (1) Alpha-Alpha block: -sum(w) / sigma^2
  hessian(0, 0) = -sum_w / sigma2;

  // (2) Beta-Beta block: -X'WX / sigma^2
  hessian.submat(1, 1, p, p) = -Xt_WX / sigma2;

  // (3) Alpha-Beta blocks: -sum(wX) / sigma^2
  // sum_wX 已經是 column vector (p x 1)
  arma::vec H_ab = -sum_wX / sigma2;
  hessian.submat(0, 1, 0, p) = H_ab.t();
  hessian.submat(1, 0, p, 0) = H_ab;

  // (4) Interaction with Sigma (Alpha-Sigma, Beta-Sigma)
  // H_as = -2 * sum(w * res) / sigma^3
  // H_bs = -2 * X'(w * res) / sigma^3
  double factor_sigma = -2.0 / sigma3;

  double H_as = factor_sigma * arma::sum(w_res);
  arma::vec H_bs = factor_sigma * (X.t() * w_res);

  hessian(0, p + 1) = H_as;
  hessian(p + 1, 0) = H_as;
  hessian.submat(1, p + 1, p, p + 1) = H_bs;
  hessian.submat(p + 1, 1, p + 1, p) = H_bs.t();

  // (5) Sigma-Sigma block
  // H_ss = sum(w)/sigma^2 - 3 * sum(w * res^2) / sigma^4
  hessian(p + 1, p + 1) = (sum_w / sigma2) - (3.0 * sum_wr2 / sigma4);

  // 5. 正規化 (除以 n) 與 檢查
  gradient /= (double)n;
  hessian /= (double)n;

  if (!gradient.is_finite()) gradient.zeros();
  if (!hessian.is_finite()) hessian.zeros();

  return List::create(
    Named("gradient") = gradient,
    Named("hessian")  = hessian
  );
}

// [[Rcpp::export]]
double lL_logistic_rcpp(const arma::mat & X,
                        const arma::vec & Y,
                        const double & alpha,
                        const arma::vec & beta,
                        Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;
  double lL = 0.0;

  // 1. 向量化計算 Linear Predictor (這是最大的效能提升點！)
  // 直接一次算出所有 eta，利用 BLAS 優化
  arma::vec eta = alpha + X * beta;

  // 2. 根據是否有權重，進入不同的迴圈 (Loop Unswitching)
  if (w.isNotNull()) {

    // --- 有權重迴圈 ---
    arma::vec weights = as<arma::vec>(w);
    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    for (arma::uword i = 0; i < n; ++i) {

      double term = 0.0;
      // 存取預先算好的 eta(i)
      if (Y(i) == 1) {

        term = -log1pexp_neg(eta(i));

      } else {

        term = -log1pexp(eta(i));
      }
      lL += term * weights(i); // 乘上權重
    }

  } else {

    // --- 無權重迴圈 ---
    for (arma::uword i = 0; i < n; ++i) {

      double term = 0.0;
      if (Y(i) == 1) {

        term = -log1pexp_neg(eta(i));

      } else {

        term = -log1pexp(eta(i));
      }
      lL += term; // 直接加總
    }
  }

  // 檢查數值有限性
  if (!std::isfinite(lL)) {
    lL = -10000.0 * (double)n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_logistic_rcpp(const arma::mat & X,
                           const arma::vec & Y,
                           const double & alpha,
                           const arma::vec & beta,
                           Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

  // 1. 計算機率向量 p
  // eta = alpha + X * beta
  arma::vec eta = alpha + X * beta;
  // p = 1 / (1 + exp(-eta))
  arma::vec p_vec = 1.0 / (1.0 + arma::exp(-eta));

  // 2. 準備權重與殘差
  // Gradient weight factor: w (or 1)
  // Hessian weight factor: w * p * (1-p) (or p * (1-p))

  arma::vec res = Y - p_vec; // Unweighted residuals initially
  arma::vec hess_w;          // Hessian diagonal weights

  if (w.isNotNull()) {

    // --- 有權重 ---
    arma::vec weights = as<arma::vec>(w);
    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    // Gradient 部分: res 變為 w * (Y - p)
    res %= weights;

    // Hessian 部分: w * p * (1-p)
    hess_w = weights % p_vec % (1.0 - p_vec);

  } else {

    // --- 無權重 ---
    // Hessian 部分: p * (1-p)
    hess_w = p_vec % (1.0 - p_vec);
  }

  // 3. 輸出容器初始化
  arma::vec gradient(p + 1, arma::fill::zeros);
  arma::mat hessian(p + 1, p + 1, arma::fill::zeros);

  // 4. 計算 Gradient
  // dL/d_alpha = sum(w * (Y - p))
  // dL/d_beta  = X' * (w * (Y - p))

  gradient(0) = arma::sum(res);                   // alpha
  gradient.subvec(1, p) = X.t() * res;            // beta

  // 5. 計算 Hessian
  // Hessian = - X_aug' * W * X_aug
  // W = diag(hess_w)

  // (1) Alpha-Alpha block: -sum(hess_w)
  hessian(0, 0) = -arma::sum(hess_w);

  // (2) Alpha-Beta blocks: -sum(hess_w * x_ij)
  // hess_w' * X  ->  1 x p vector
  arma::rowvec H_ab = - (hess_w.t() * X);
  hessian.submat(0, 1, 0, p) = H_ab;
  hessian.submat(1, 0, p, 0) = H_ab.t();

  // (3) Beta-Beta block: -X' * diag(hess_w) * X
  // 優化技巧：避免建立對角矩陣
  // 計算 X_weighted = X 乘上 column-wise weights
  arma::mat X_weighted = X;
  X_weighted.each_col() %= hess_w;

  hessian.submat(1, 1, p, p) = -X.t() * X_weighted;

  // 6. 正規化 (除以 n) 與 檢查
  gradient /= (double)n;
  hessian /= (double)n;

  if (!gradient.is_finite()) gradient.zeros();
  if (!hessian.is_finite()) hessian.zeros();

  return List::create(
    Named("gradient") = gradient,
    Named("hessian")  = hessian
  );
}

// [[Rcpp::export]]
double lL_Gamma_rcpp(const arma::mat & X,
                     const arma::vec & Y,
                     const double & alpha,
                     const arma::vec & beta,
                     const double & nu,
                     Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;

  // 1. 向量化計算 Linear Predictor
  // eta = alpha + X * beta
  arma::vec eta = alpha + X * beta;

  // 2. 計算中間變數 (向量化)
  // lambda = nu * exp(-eta)
  arma::vec lambda = nu * arma::exp(-eta);

  // log(lambda) = log(nu * exp(-eta)) = log(nu) - eta
  // 優化技巧：直接用減法，避免對整個向量做 log運算，速度大幅提升
  arma::vec log_lambda = std::log(nu) - eta;

  // log(Y)
  arma::vec log_Y = arma::log(Y);

  // 3. 準備統計量變數
  double term_log_lambda = 0.0; // sum(w * log(lambda))
  double term_log_Y      = 0.0; // sum(w * log(Y))
  double term_lambda_Y   = 0.0; // sum(w * lambda * Y)
  double sum_w           = 0.0; // sum(w) for lgamma term

  // 4. 根據權重計算各項總和 (使用 dot product 加速)
  if (w.isNotNull()) {

    // --- 有權重 ---
    arma::vec weights = as<arma::vec>(w);
    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    sum_w = arma::sum(weights);

    // 利用內積計算加權和，取代迴圈
    term_log_lambda = arma::dot(weights, log_lambda);
    term_log_Y      = arma::dot(weights, log_Y);

    // lambda * Y * w
    // 先算 (lambda % Y) 得到 element-wise 乘積，再跟 weights 做 dot
    term_lambda_Y   = arma::dot(weights, lambda % Y);

  } else {

    // --- 無權重 ---
    sum_w = (double)n;

    term_log_lambda = arma::sum(log_lambda);
    term_log_Y      = arma::sum(log_Y);
    term_lambda_Y   = arma::dot(lambda, Y); // dot(a, b) = sum(a * b)
  }

  // 5. 組合 Log-Likelihood
  // Formula: nu * sum(log_lambda) + (nu - 1) * sum(log_Y) - sum(lambda * Y) - sum(w)*lgamma(nu)
  double lL = nu * term_log_lambda +
    (nu - 1.0) * term_log_Y -
    term_lambda_Y -
    sum_w * std::lgamma(nu);

  // 6. 數值檢查
  if (!std::isfinite(lL)) {
    lL = -10000.0 * (double)n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_Gamma_rcpp(const arma::mat & X,
                        const arma::vec & Y,
                        const double & alpha,
                        const arma::vec & beta,
                        const double & nu,
                        Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

  // 1. 向量化計算中間變數
  // eta = alpha + X * beta
  arma::vec eta = alpha + X * beta;

  // eSI = exp(eta)
  arma::vec eSI = arma::exp(eta);

  // lambda = nu / eSI
  arma::vec lambda = nu / eSI;

  // 2. 準備權重與統計量
  double sum_w = 0.0;     // sum(w)
  arma::vec w_vec;        // 用於儲存 weights (若有)
  arma::vec res_beta;     // 用於 beta gradient 的殘差項
  arma::vec term_nu;      // 用於 nu gradient 的項
  arma::vec hess_w_beta;  // 用於 beta hessian 的權重

  // 基礎項計算
  // res_beta_base = Y * lambda - nu
  arma::vec res_beta_base = Y % lambda - nu;

  // term_nu_base = log(lambda * Y) - Y / eSI
  // 優化：log(lambda * Y) = log(nu/eSI * Y) = log(nu) - eta + log(Y)
  arma::vec term_nu_base = std::log(nu) - eta + arma::log(Y) - (Y / eSI);

  if (w.isNotNull()) {

    // --- 有權重 ---
    w_vec = as<arma::vec>(w);
    if (w_vec.n_elem != n) stop("Weights length differs from X rows.");

    sum_w = arma::sum(w_vec);

    // Gradient parts applied with weights
    res_beta = w_vec % res_beta_base;
    term_nu  = w_vec % term_nu_base;

    // Hessian weight for beta: w * Y * lambda
    hess_w_beta = w_vec % Y % lambda;

  } else {

    // --- 無權重 ---
    sum_w = (double)n;

    res_beta = res_beta_base;
    term_nu  = term_nu_base;

    // Hessian weight for beta: Y * lambda
    hess_w_beta = Y % lambda;
  }

  // 3. 輸出容器初始化
  arma::vec gradient(p + 2, arma::fill::zeros);
  arma::mat hessian(p + 2, p + 2, arma::fill::zeros);

  // 4. 計算 Gradient
  // ---------------------------------------------------
  // (1) Alpha (Intercept): sum(w * (Y*lambda - nu))
  gradient(0) = arma::sum(res_beta);

  // (2) Beta: X' * (w * (Y*lambda - nu))
  gradient.subvec(1, p) = X.t() * res_beta;

  // (3) Nu: sum(w * term_nu) + sum(w) * (1 - digamma(nu))
  // 原代碼是在最後加 (1 - digamma)，這裡對應 sum_w
  gradient(p + 1) = arma::sum(term_nu) + sum_w * (1.0 - R::digamma(nu));


  // 5. 計算 Hessian
  // ---------------------------------------------------
  // H_beta_block = - X_aug' * W * X_aug
  // Weight for Hessian: hess_w_beta = w * Y * lambda

  // (1) Alpha-Alpha: -sum(hess_w_beta)
  hessian(0, 0) = -arma::sum(hess_w_beta);

  // (2) Alpha-Beta: -sum(hess_w_beta * X)
  // hess_w_beta' * X -> 1 x p vector
  arma::rowvec H_ab = - (hess_w_beta.t() * X);
  hessian.submat(0, 1, 0, p) = H_ab;
  hessian.submat(1, 0, p, 0) = H_ab.t();

  // (3) Beta-Beta: -X' * diag(hess_w_beta) * X
  arma::mat X_weighted = X;
  X_weighted.each_col() %= hess_w_beta; // column-wise multiplication
  hessian.submat(1, 1, p, p) = -X.t() * X_weighted;

  // (4) Nu-Nu: sum(w) * (1/nu - trigamma(nu))
  hessian(p + 1, p + 1) = sum_w * (1.0 / nu - R::trigamma(nu));


  // 6. 正規化 (除以 n) 與 檢查
  gradient /= (double)n;
  hessian /= (double)n;

  if (!gradient.is_finite()) gradient.zeros();
  if (!hessian.is_finite()) hessian.zeros();

  return List::create(
    Named("gradient") = gradient,
    Named("hessian")  = hessian
  );
}

// [[Rcpp::export]]
double lL_binomial_rcpp(const arma::mat & X,
                        const arma::vec & Y,
                        const double & alpha,
                        const arma::vec & beta,
                        const double & N, // 建議用 double 避免與 R 互動時的整數溢位問題
                        Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;
  double lL = 0.0;

  // 1. 向量化計算 Linear Predictor
  // eta = alpha + X * beta
  arma::vec eta = alpha + X * beta;

  // 2. Loop Unswitching (將權重判斷移出迴圈)
  if (w.isNotNull()) {

    // --- 有權重 ---
    arma::vec weights = as<arma::vec>(w);
    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    for (arma::uword i = 0; i < n; ++i) {
      // 公式：w * [ Y * eta - N * log(1 + exp(eta)) ]
      double term = Y(i) * eta(i) - N * log1pexp(eta(i));
      lL += term * weights(i);
    }

  } else {

    // --- 無權重 ---
    for (arma::uword i = 0; i < n; ++i) {
      // 公式：Y * eta - N * log(1 + exp(eta))
      double term = Y(i) * eta(i) - N * log1pexp(eta(i));
      lL += term;
    }
  }

  // 3. 數值檢查
  if (!std::isfinite(lL)) {
    lL = -10000.0 * (double)n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_binomial_rcpp(const arma::mat & X,
                           const arma::vec & Y,
                           const double & alpha,
                           const arma::vec & beta,
                           const double & N, // 改用 double 避免計算時的轉型問題
                           Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

  // 1. 向量化計算機率
  // eta = alpha + X * beta
  arma::vec eta = alpha + X * beta;

  // p = 1 / (1 + exp(-eta))
  arma::vec p_vec = 1.0 / (1.0 + arma::exp(-eta));

  // 2. 準備權重與中間變數
  // res: 用於 gradient 的殘差向量 (Y - N*p)
  // hess_w: 用於 hessian 的對角線權重 N * p * (1-p)

  arma::vec res;
  arma::vec hess_w;

  // 基礎殘差與變異數項
  // Residual base = Y - N * p
  arma::vec res_base = Y - N * p_vec;
  // Variance base = N * p * (1-p)
  arma::vec var_base = N * p_vec % (1.0 - p_vec);

  if (w.isNotNull()) {

    // --- 有權重 ---
    arma::vec weights = as<arma::vec>(w);
    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    // Gradient 殘差: w * (Y - N*p)
    res = weights % res_base;

    // Hessian 權重: w * N * p * (1-p)
    hess_w = weights % var_base;

  } else {

    // --- 無權重 ---
    res = res_base;
    hess_w = var_base;
  }

  // 3. 輸出容器初始化
  arma::vec gradient(p + 1, arma::fill::zeros);
  arma::mat hessian(p + 1, p + 1, arma::fill::zeros);

  // 4. 計算 Gradient (利用矩陣運算)
  // dL/d_alpha = sum(res)
  // dL/d_beta  = X' * res

  gradient(0) = arma::sum(res);                   // alpha
  gradient.subvec(1, p) = X.t() * res;            // beta

  // 5. 計算 Hessian (利用矩陣運算)
  // Hessian = - X_aug' * W * X_aug
  // W = diag(hess_w)

  // (1) Alpha-Alpha block: -sum(hess_w)
  hessian(0, 0) = -arma::sum(hess_w);

  // (2) Alpha-Beta blocks: -sum(hess_w * X)
  // hess_w' * X -> 1 x p row vector
  arma::rowvec H_ab = - (hess_w.t() * X);
  hessian.submat(0, 1, 0, p) = H_ab;
  hessian.submat(1, 0, p, 0) = H_ab.t();

  // (3) Beta-Beta block: -X' * diag(hess_w) * X
  // 優化技巧：計算 X_weighted = X 乘上 column-wise weights
  arma::mat X_weighted = X;
  X_weighted.each_col() %= hess_w;

  hessian.submat(1, 1, p, p) = -X.t() * X_weighted;

  // 6. 正規化 (除以 n) 與 檢查
  gradient /= (double)n;
  hessian /= (double)n;

  if (!gradient.is_finite()) gradient.zeros();
  if (!hessian.is_finite()) hessian.zeros();

  return List::create(
    Named("gradient") = gradient,
    Named("hessian")  = hessian
  );
}

// [[Rcpp::export]]
double lL_Poisson_rcpp(const arma::mat & X,
                       const arma::vec & Y,
                       const double & alpha,
                       const arma::vec & beta,
                       Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;

  // 1. 向量化計算 Linear Predictor
  // SI (eta) = alpha + X * beta
  // 這裡利用矩陣乘法，比迴圈快
  arma::vec SI = alpha + X * beta;

  // 2. 向量化計算 Lambda
  // lambda = exp(SI)
  arma::vec lambda = arma::exp(SI);

  // 3. 計算 Log-Likelihood
  double lL = 0.0;

  if (w.isNotNull()) {

    // --- 有權重 ---
    arma::vec weights = as<arma::vec>(w);
    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    // 公式: sum( w * (Y * SI - lambda) )
    // 這可以用內積 (dot product) 一步完成
    // Y % SI 是 element-wise 乘法
    arma::vec term = Y % SI - lambda;
    lL = arma::dot(weights, term);

  } else {

    // --- 無權重 ---
    // 公式: sum( Y * SI - lambda )
    // 拆解為 dot(Y, SI) - sum(lambda) 會比先相減再相加稍快一些
    lL = arma::dot(Y, SI) - arma::sum(lambda);
  }

  // 4. 數值檢查
  if (!std::isfinite(lL)) {
    lL = -10000.0 * (double)n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_Poisson_rcpp(const arma::mat & X,
                          const arma::vec & Y,
                          const double & alpha,
                          const arma::vec & beta,
                          Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

  // 1. 向量化計算 Lambda
  // eta = alpha + X * beta
  arma::vec eta = alpha + X * beta;

  // lambda = exp(eta)
  arma::vec lambda = arma::exp(eta);

  // 2. 準備權重與中間變數
  // res: 用於 Gradient 的殘差項 (Y - lambda)
  // hess_w: 用於 Hessian 的權重項 lambda

  arma::vec res;
  arma::vec hess_w;

  // 基礎殘差
  arma::vec res_base = Y - lambda;

  if (w.isNotNull()) {

    // --- 有權重 ---
    arma::vec weights = as<arma::vec>(w);
    if (weights.n_elem != n) stop("Weights length differs from X rows.");

    // Gradient 殘差: w * (Y - lambda)
    res = weights % res_base;

    // Hessian 權重: w * lambda
    hess_w = weights % lambda;

  } else {

    // --- 無權重 ---
    res = res_base;
    hess_w = lambda;
  }

  // 3. 輸出容器初始化
  arma::vec gradient(p + 1, arma::fill::zeros);
  arma::mat hessian(p + 1, p + 1, arma::fill::zeros);

  // 4. 計算 Gradient (矩陣運算)
  // dL/d_alpha = sum(res)
  // dL/d_beta  = X' * res

  gradient(0) = arma::sum(res);                   // alpha
  gradient.subvec(1, p) = X.t() * res;            // beta

  // 5. 計算 Hessian (矩陣運算)
  // Hessian = - X_aug' * W * X_aug
  // W = diag(hess_w)

  // (1) Alpha-Alpha block: -sum(hess_w)
  hessian(0, 0) = -arma::sum(hess_w);

  // (2) Alpha-Beta blocks: -sum(hess_w * X)
  // hess_w' * X -> 1 x p row vector
  arma::rowvec H_ab = - (hess_w.t() * X);
  hessian.submat(0, 1, 0, p) = H_ab;
  hessian.submat(1, 0, p, 0) = H_ab.t();

  // (3) Beta-Beta block: -X' * diag(hess_w) * X
  // 優化技巧：避免建立巨大對角矩陣
  // 計算 X_weighted = X 乘上 column-wise weights
  arma::mat X_weighted = X;
  X_weighted.each_col() %= hess_w;

  hessian.submat(1, 1, p, p) = -X.t() * X_weighted;

  // 6. 正規化 (除以 n) 與 檢查
  gradient /= (double)n;
  hessian /= (double)n;

  if (!gradient.is_finite()) gradient.zeros();
  if (!hessian.is_finite()) hessian.zeros();

  return List::create(
    Named("gradient") = gradient,
    Named("hessian")  = hessian
  );
}
















//============================================================================

// [[Rcpp::export]]
double lL_Poisson_rcpp_v1(const arma::mat & X,
                          const arma::vec & Y,
                          const double & alpha,
                          const arma::vec & beta) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  double SI_i = 0.0;
  double lambda_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    SI_i = alpha + arma::dot(X.row(i).t(), beta);
    lambda_i = exp(SI_i);

    lL += Y(i) * SI_i - lambda_i;
  }

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
double lL_Poisson_w_rcpp_v1(const arma::mat & X,
                            const arma::vec & Y,
                            const double & alpha,
                            const arma::vec & beta,
                            const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  double SI_i = 0.0;
  double lambda_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    SI_i = alpha + arma::dot(X.row(i).t(), beta);
    lambda_i = exp(SI_i);

    lL += (Y(i) * SI_i - lambda_i) * w(i);
  }

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_Poisson_rcpp_v1(const arma::mat & X,
                             const arma::vec & Y,
                             const double & alpha,
                             const arma::vec & beta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec gradient(n_p + 1);
  arma::mat hessian(n_p + 1, n_p + 1);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double lambda_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    lambda_i = exp(alpha + arma::dot(X.row(i).t(), beta));

    gradient += extXrow_i * (Y(i) - lambda_i);
    hessian += -extXrow_i * extXrow_i.t() * lambda_i;
  }

  gradient /= n_n;
  hessian /= n_n;

  if (!gradient.is_finite()) {

    gradient = arma::vec(n_p + 1);
  }

  if (!hessian.is_finite()) {

    hessian = arma::mat(n_p + 1, n_p + 1);
  }

  return Rcpp::List::create(
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}
