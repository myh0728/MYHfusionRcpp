#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

extern List eigen_rcpp(const arma::mat & M_S);

// [[Rcpp::export]]
List cumuSIR_rcpp(const arma::mat & X,
                  const arma::vec & Y,
                  const double & eps_inv) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;

  arma::vec m_y(n_p);
  arma::mat Km(n_p, n_p);

  List varX_eigen = eigen_rcpp(arma::cov(X));
  arma::vec varX_eigen_value = varX_eigen["value"];
  arma::mat varX_eigen_vector = varX_eigen["vector"];
  varX_eigen_value.elem(arma::find(varX_eigen_value < eps_inv)).fill(eps_inv);
  arma::mat normalizing = varX_eigen_vector *
    arma::diagmat(pow(varX_eigen_value, -0.5)) * varX_eigen_vector.t();
  arma::rowvec meanX = arma::mean(X, 0);
  arma::mat X_cs = (X.each_row() - meanX) * normalizing.t();

  for (size_t j = 0; j < n_n; ++j) {

    m_y.zeros();

    for (size_t i = 0; i < n_n; ++i) {

      m_y += X_cs.row(i).t() * (Y(i) <= Y(j)) / n_n;
    }

    Km += m_y * m_y.t() / n_n;
  }

  List RR = eigen_rcpp(Km);
  arma::mat Bhat = normalizing * as<arma::mat>(RR["vector"]);

  return Rcpp::List::create(
    Rcpp::Named("basis") = Bhat,
    Rcpp::Named("values") = RR["value"]
  );
}








