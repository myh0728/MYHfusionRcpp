#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List cumuSIR_rcpp(const arma::mat & X,
                  const arma::vec & Y,
                  const double & eps_inv) {

  const arma::uword n_n = X.n_rows;

  // 1. 標準化 X (Standardization)
  // 計算 Cov(X) 並進行特徵值分解
  arma::vec varX_eigen_value;
  arma::mat varX_eigen_vector;
  arma::mat covX = arma::cov(X);

  // 處理數值穩定性：如果 Cov 矩陣接近奇異矩陣
  eig_sym(varX_eigen_value, varX_eigen_vector, covX);

  // 處理過小的特徵值 (Regularization)
  varX_eigen_value.elem(arma::find(varX_eigen_value < eps_inv)).fill(eps_inv);

  // 計算 Sigma^{-1/2}
  arma::mat normalizing = varX_eigen_vector *
    arma::diagmat(pow(varX_eigen_value, -0.5)) * varX_eigen_vector.t();

  arma::rowvec meanX = arma::mean(X, 0);
  arma::mat X_cs = (X.each_row() - meanX) * normalizing;

  // 2. 排序 (Sorting) - 這是加速的關鍵
  // 找出 Y 從小到大的索引
  arma::uvec sorted_idx = arma::sort_index(Y);

  // 根據 Y 的順序重新排列 X_cs
  arma::mat X_sorted = X_cs.rows(sorted_idx);

  // 3. 計算候選矩陣 (Candidate Matrix)
  arma::mat Km(X.n_cols, X.n_cols, arma::fill::zeros);
  arma::vec cumsum_x(X.n_cols, arma::fill::zeros);

  double inv_n = 1.0 / n_n; // 預先計算倒數

  for (arma::uword j = 0; j < n_n; ++j) {
    // 累積 X 的值
    cumsum_x += X_sorted.row(j).t();
    arma::vec m_y = cumsum_x * inv_n;

    // 累積到 Km
    Km += m_y * m_y.t() * inv_n; // 這裡多乘一次 1/n_n 是對外層迴圈的平均
  }

  // 4. 對 Km 進行特徵值分解
  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, Km);

  // 排序特徵值由大到小
  eigval = arma::reverse(eigval);
  eigvec = arma::reverse(eigvec, 1);

  // 轉換回原始空間
  arma::mat Bhat = normalizing * eigvec;

  return Rcpp::List::create(
    Rcpp::Named("basis") = Bhat,
    Rcpp::Named("values") = eigval
  );
}

