#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double lL_normal_rcpp(const arma::mat & X,
                      const arma::vec & Y,
                      const double & alpha,
                      const arma::vec & beta,
                      const double & sigma) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    lL += pow(Y(i) - alpha - arma::dot(X.row(i).t(), beta), 2);
  }

  lL = -lL / (2 * pow(sigma, 2)) -
    n_n * log(2 * arma::datum::pi * pow(sigma, 2)) / 2;

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
double lL_normal_w_rcpp(const arma::mat & X,
                        const arma::vec & Y,
                        const double & alpha,
                        const arma::vec & beta,
                        const double & sigma,
                        const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    lL += pow(Y(i) - alpha - arma::dot(X.row(i).t(), beta), 2) * w(i);
  }

  lL = -lL / (2 * pow(sigma, 2)) -
    log(2 * arma::datum::pi * pow(sigma, 2)) * sum(w) / 2;

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_normal_rcpp(const arma::mat & X,
                         const arma::vec & Y,
                         const double & alpha,
                         const arma::vec & beta,
                         const double & sigma) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec gradient(n_p + 2);
  arma::mat hessian(n_p + 2, n_p + 2);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double res_i = 0.0;
  double res_i_square = 0.0;
  arma::vec cross(n_p + 1);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    res_i = Y(i) - alpha - arma::dot(Xrow_i, beta);
    res_i_square = pow(res_i, 2);
    cross = extXrow_i * res_i;

    gradient.subvec(0, n_p) += cross;
    gradient(n_p + 1) += res_i_square;
    hessian.submat(0, 0, n_p, n_p) += -extXrow_i * extXrow_i.t();
    hessian(n_p + 1, n_p + 1) += -res_i_square;
    hessian.submat(0, n_p + 1, n_p, n_p + 1) += -cross;
    hessian.submat(n_p + 1, 0, n_p + 1, n_p) += -cross.t();
  }

  gradient /= n_n;
  gradient.subvec(0, n_p) /= pow(sigma, 2);
  gradient(n_p + 1) /= pow(sigma, 3);
  gradient(n_p + 1) -= pow(sigma, -1);
  hessian /= n_n;
  hessian.submat(0, 0, n_p, n_p) /= pow(sigma, 2);
  hessian(n_p + 1, n_p + 1) *= 3 / pow(sigma, 4);
  hessian(n_p + 1, n_p + 1) += pow(sigma, -2);
  hessian.submat(0, n_p + 1, n_p, n_p + 1) *= 2 / pow(sigma, 3);
  hessian.submat(n_p + 1, 0, n_p + 1, n_p) *= 2 / pow(sigma, 3);

  if (!gradient.is_finite()) {

    gradient = arma::vec(n_p + 2);
  }

  if (!hessian.is_finite()) {

    hessian = arma::mat(n_p + 2, n_p + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
double lL_logistic_rcpp(const arma::mat & X,
                        const arma::vec & Y,
                        const double & alpha,
                        const arma::vec & beta) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    if (Y(i) == 1) {

      lL += -log(1 + exp(-(alpha + arma::dot(X.row(i).t(), beta))));

    }else{

      lL += -log(1 + exp(alpha + arma::dot(X.row(i).t(), beta)));
    }
  }

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
double lL_logistic_w_rcpp(const arma::mat & X,
                          const arma::vec & Y,
                          const double & alpha,
                          const arma::vec & beta,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    if (Y(i) == 1) {

      lL += -log(1 + exp(-(alpha + arma::dot(X.row(i).t(), beta)))) * w(i);

    }else{

      lL += -log(1 + exp(alpha + arma::dot(X.row(i).t(), beta))) * w(i);
    }
  }

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_logistic_rcpp(const arma::mat & X,
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
  double eSI_i = 0.0;
  double p1_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    p1_i = eSI_i / (1 + eSI_i);

    gradient += extXrow_i * (Y(i) - p1_i);
    hessian += -extXrow_i * extXrow_i.t() * p1_i / (1 + eSI_i);
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

// [[Rcpp::export]]
double lL_Gamma_rcpp(const arma::mat & X,
                     const arma::vec & Y,
                     const double & alpha,
                     const arma::vec & beta,
                     const double & nu) {

  const arma::uword n_n = X.n_rows;
  double lL_lambda = 0.0;
  double lL_Y = 0.0;
  double lL_lambdaY = 0.0;
  double lL = 0.0;

  double lambda_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    lambda_i = nu * exp(-(alpha + arma::dot(X.row(i).t(), beta)));
    lL_lambda += log(lambda_i);
    lL_Y += log(Y(i));
    lL_lambdaY += lambda_i * Y(i);
  }

  lL = nu * lL_lambda + (nu - 1) * lL_Y - lL_lambdaY - n_n * lgamma(nu);

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
double lL_Gamma_w_rcpp(const arma::mat & X,
                       const arma::vec & Y,
                       const double & alpha,
                       const arma::vec & beta,
                       const double & nu,
                       const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  double lL_lambda = 0.0;
  double lL_Y = 0.0;
  double lL_lambdaY = 0.0;
  double lL = 0.0;

  double lambda_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    lambda_i = nu * exp(-(alpha + arma::dot(X.row(i).t(), beta)));
    lL_lambda += log(lambda_i) * w(i);
    lL_Y += log(Y(i)) * w(i);
    lL_lambdaY += lambda_i * Y(i) * w(i);
  }

  lL = nu * lL_lambda + (nu - 1) * lL_Y - lL_lambdaY - sum(w) * lgamma(nu);

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_Gamma_rcpp(const arma::mat & X,
                        const arma::vec & Y,
                        const double & alpha,
                        const arma::vec & beta,
                        const double & nu) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec gradient(n_p + 2);
  arma::mat hessian(n_p + 2, n_p + 2);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double lambda_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(X.row(i).t(), beta));
    lambda_i = nu / eSI_i;

    gradient.subvec(0, n_p) += extXrow_i * (Y(i) * lambda_i - nu);
    gradient(n_p + 1) += log(lambda_i * Y(i)) - Y(i) / eSI_i;
    hessian.submat(0, 0, n_p, n_p) += -Y(i) * lambda_i * extXrow_i * extXrow_i.t();
  }

  gradient /= n_n;
  gradient(n_p + 1) += 1 - R::digamma(nu);
  hessian /= n_n;
  hessian(n_p + 1, n_p + 1) = 1 / nu - R::trigamma(nu);

  if (!gradient.is_finite()) {

    gradient = arma::vec(n_p + 2);
  }

  if (!hessian.is_finite()) {

    hessian = arma::mat(n_p + 2, n_p + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}

// [[Rcpp::export]]
double lL_binomial_rcpp(const arma::mat & X,
                        const arma::vec & Y,
                        const double & alpha,
                        const arma::vec & beta,
                        const size_t & N) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  double p_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    p_i = 1 / (1 + exp(-(alpha + arma::dot(X.row(i).t(), beta))));
    lL += Y(i) * log(p_i) + (N - Y(i)) * log(1 - p_i);
  }

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
double lL_binomial_w_rcpp(const arma::mat & X,
                          const arma::vec & Y,
                          const double & alpha,
                          const arma::vec & beta,
                          const size_t & N,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  double lL = 0.0;

  double p_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    p_i = 1 / (1 + exp(-(alpha + arma::dot(X.row(i).t(), beta))));
    lL += (Y(i) * log(p_i) + (N - Y(i)) * log(1 - p_i)) * w(i);
  }

  if (!std::isfinite(lL)) {

    lL = -10000.0 * n_n;
  }

  return lL;
}

// [[Rcpp::export]]
List diff_lL_binomial_rcpp(const arma::mat & X,
                           const arma::vec & Y,
                           const double & alpha,
                           const arma::vec & beta,
                           const size_t & N) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec gradient(n_p + 1);
  arma::mat hessian(n_p + 1, n_p + 1);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double p_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    p_i = 1 / (1 + exp(-(alpha + arma::dot(X.row(i).t(), beta))));

    gradient += extXrow_i * (Y(i) - N * p_i);
    hessian += -extXrow_i * extXrow_i.t() * p_i * (1 - p_i);
  }

  gradient /= n_n;
  hessian *= N;
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

// [[Rcpp::export]]
double lL_Poisson_rcpp(const arma::mat & X,
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
double lL_Poisson_w_rcpp(const arma::mat & X,
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
List diff_lL_Poisson_rcpp(const arma::mat & X,
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



// For comparison

// [[Rcpp::export]]
double lL_normal_rcpp_fun(const arma::mat & X,
                          const arma::vec & Y,
                          const double & alpha,
                          const arma::vec & beta,
                          const double & sigma) {

  const arma::vec residual = Y - X * beta - alpha;
  const arma::vec log_likelihood = arma::log_normpdf(residual, 0, sigma);
  double lL = sum(log_likelihood);

  return lL;
}

// [[Rcpp::export]]
double lL_normal_w_rcpp_fun(const arma::mat & X,
                            const arma::vec & Y,
                            const double & alpha,
                            const arma::vec & beta,
                            const double & sigma,
                            const arma::vec & w) {

  const arma::vec residual = Y - X * beta - alpha;
  const arma::vec log_likelihood = arma::log_normpdf(residual, 0, sigma);
  double lL = arma::dot(log_likelihood, w);

  return lL;
}

// [[Rcpp::export]]
List diff_lL_normal_rcpp_v1(const arma::mat & X,
                            const arma::vec & Y,
                            const double & alpha,
                            const arma::vec & beta,
                            const double & sigma) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec gradient(n_p + 2);
  arma::mat hessian(n_p + 2, n_p + 2);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double res_i = Y(i) - alpha - arma::dot(Xrow_i, beta);
    const double res_i_square = pow(res_i, 2);
    const arma::vec cross = extXrow_i * res_i;
    gradient.subvec(0, n_p) += cross;
    gradient(n_p + 1) += res_i_square;
    hessian.submat(0, 0, n_p, n_p) += -extXrow_i * extXrow_i.t();
    hessian(n_p + 1, n_p + 1) += -res_i_square;
    hessian.submat(0, n_p + 1, n_p, n_p + 1) += -cross;
    hessian.submat(n_p + 1, 0, n_p + 1, n_p) += -cross.t();
  }

  gradient /= n_n;
  gradient.subvec(0, n_p) /= pow(sigma, 2);
  gradient(n_p + 1) /= pow(sigma, 3);
  gradient(n_p + 1) -= pow(sigma, -1);
  hessian /= n_n;
  hessian.submat(0, 0, n_p, n_p) /= pow(sigma, 2);
  hessian(n_p + 1, n_p + 1) *= 3 / pow(sigma, 4);
  hessian(n_p + 1, n_p + 1) += pow(sigma, -2);
  hessian.submat(0, n_p + 1, n_p, n_p + 1) *= 2 / pow(sigma, 3);
  hessian.submat(n_p + 1, 0, n_p + 1, n_p) *= 2 / pow(sigma, 3);

  if (!gradient.is_finite()) {

    gradient = arma::vec(n_p + 2);
  }

  if (!hessian.is_finite()) {

    hessian = arma::mat(n_p + 2, n_p + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = hessian
  );
}





