#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// No shift, normal regression model

// [[Rcpp::export]]
List AD_EY_normal_Lagrange_rcpp(const arma::mat & X,
                                const double & alpha,
                                const arma::vec & beta,
                                const double & sigma,
                                const double & phi,
                                const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double Psi_i = 0.0;
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    Psi_i = SI_i - phi;

    denominator = 1.0 + eta * Psi_i;
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
List AD_EY_normal_SolveLagrange_rcpp(const arma::mat & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const double & sigma,
                                     const double & phi,
                                     const double & eta_initial,
                                     const size_t & iter_max,
                                     const double & step_rate,
                                     const size_t & step_max,
                                     const double & tol,
                                     const double & eps_inv) {

  double eta = eta_initial;
  List step = AD_EY_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, eta);

  double step_value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  double direction_step = 0.0;
  double step_size = 0.0;
  double eta_new = 0.0;
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = step["gradient"];
    hessian = step["hessian"];

    if (abs(hessian) > eps_inv) {

      direction_step = gradient / hessian;

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = AD_EY_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EY_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List AD_EY_normal_rcpp(const arma::mat & X,
                       const double & alpha,
                       const arma::vec & beta,
                       const double & sigma,
                       const double & phi) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p + 3);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double Psi_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    Psi_i = SI_i - phi;

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i;
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t();
    Psi_gradient(0, n_p + 2) += -1;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(1);
    Psi_square = arma::mat(1, 1);
    Psi_gradient = arma::mat(1, n_p + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADvar_EY_normal_rcpp(const arma::mat & X,
                          const double & alpha,
                          const arma::vec & beta,
                          const double & sigma,
                          const double & phi,
                          const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double Psi_i = 0.0;
  double var_A = 0.0;
  double var_B = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    Psi_i = SI_i - phi;
    var_B += pow(Psi_i, 2) / n_n / (1.0 + Psi_i * eta);
  }

  var_A = -1.0;
  var_B += pow(sigma, 2);

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List AD_EXsubY_normal_Lagrange_rcpp(const arma::mat & X,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & sigma,
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

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::mat dist_e_f(n_k, n_p);
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    q_pts = (y_pts - SI_i) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    dist_e_f = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    Psi_i = arma::reshape(
      (dist_e_f % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List AD_EXsubY_normal_SolveLagrange_rcpp(const arma::mat & X,
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

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = AD_EXsubY_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, y_pts, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = AD_EXsubY_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, y_pts, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EXsubY_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, y_pts, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List AD_EXsubY_normal_rcpp(const arma::mat & X,
                           const double & alpha,
                           const arma::vec & beta,
                           const double & sigma,
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

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::mat pdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::vec pdf_dist(n_k);
  arma::mat dist(n_k, n_p);
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    q_pts = (y_pts - SI_i) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    pdf_pts = arma::normpdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    Psi_i = arma::reshape((dist % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
    Psi_gradient.cols(0, n_p) += arma::reshape((dist %
      arma::repmat(-pdf_dist / sigma, 1, n_p)).t(), n_m, 1) *
      extXrow_i.t();
    Psi_gradient.col(n_p + 1) += arma::reshape((dist %
      arma::repmat(-(q_pts.col(1) / sigma) % pdf_pts.col(1) +
      (q_pts.col(0) / sigma) % pdf_pts.col(0), 1, n_p)).t(), n_m, 1);
    Psi_gradient.cols(n_p + 2, n_p + n_m + 1).diag() +=
      -arma::reshape(arma::repmat(cdf_dist, 1, n_p).t(), n_m, 1);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p + n_m + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADvar_EXsubY_normal_rcpp(const arma::mat & X,
                              const double & alpha,
                              const arma::vec & beta,
                              const double & sigma,
                              const arma::mat & phi,
                              const arma::mat & y_pts,
                              const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::mat pdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::vec pdf_dist(n_k);
  arma::mat dist(n_k, n_p);
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    q_pts = (y_pts - SI_i) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    pdf_pts = arma::normpdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    Psi_i = arma::reshape((dist % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    var_A.diag() += -arma::reshape(arma::repmat(
      cdf_dist, 1, n_p).t(), n_m, 1) * p_i;
    for (size_t k = 0; k < n_k; ++k) {

      var_B.submat(k * n_p, k * n_p,
                   (k + 1) * n_p - 1,
                   (k + 1) * n_p - 1) += dist.row(k).t() *
                     dist.row(k) * cdf_dist(k) * p_i;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List AD_EYsubX_normal_Lagrange_rcpp(const arma::mat & X,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & sigma,
                                    const arma::vec & phi,
                                    const arma::umat & inclusion,
                                    const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    Psi_i = (SI_i - phi) % inclusion.row(i).t();

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List AD_EYsubX_normal_SolveLagrange_rcpp(const arma::mat & X,
                                         const double & alpha,
                                         const arma::vec & beta,
                                         const double & sigma,
                                         const arma::vec & phi,
                                         const arma::umat & inclusion,
                                         const arma::vec & eta_initial,
                                         const size_t & iter_max,
                                         const double & step_rate,
                                         const size_t & step_max,
                                         const double & tol,
                                         const double & eps_inv) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = AD_EYsubX_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, inclusion, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = AD_EYsubX_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EYsubX_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, inclusion, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List AD_EYsubX_normal_rcpp(const arma::mat & X,
                           const double & alpha,
                           const arma::vec & beta,
                           const double & sigma,
                           const arma::vec & phi,
                           const arma::umat & inclusion) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 2);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    Psi_i.zeros();

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) = SI_i - phi(m);
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t();
        Psi_gradient(m, n_p + 2 + m) += -1;
      }
    }

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p + n_m + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADvar_EYsubX_normal_rcpp(const arma::mat & X,
                              const double & alpha,
                              const arma::vec & beta,
                              const double & sigma,
                              const arma::vec & phi,
                              const arma::umat & inclusion,
                              const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::vec inclusion_i(n_m);
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);
  arma::mat var_B_a(n_m, n_m);
  arma::mat var_B_b(n_m, n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    Psi_i = (SI_i - phi) % inclusion.row(i).t();
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        var_A(m, m) += -p_i;
      }
    }
    inclusion_i = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
    var_B_a += inclusion_i * inclusion_i.t() * p_i;
    var_B_b += Psi_i * Psi_i.t() * p_i;
  }
  var_B = var_B_a * pow(sigma, 2) + var_B_b;

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// Covariate shift, normal regression model

// [[Rcpp::export]]
List ADCS_EY_normal_Lagrange_rcpp(const arma::mat & X,
                                  const double & alpha,
                                  const arma::vec & beta,
                                  const double & sigma,
                                  const double & phi,
                                  const arma::vec & CS_beta,
                                  const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  double Psi_i = 0.0;
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i = eSI_CS_i * (SI_i - phi);

    denominator = 1.0 + eta * Psi_i;
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
List ADCS_EY_normal_SolveLagrange_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & sigma,
                                       const double & phi,
                                       const arma::vec & CS_beta,
                                       const double & eta_initial,
                                       const size_t & iter_max,
                                       const double & step_rate,
                                       const size_t & step_max,
                                       const double & tol,
                                       const double & eps_inv) {

  double eta = eta_initial;
  List step = ADCS_EY_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, CS_beta, eta);
  double step_value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  double direction_step = 0.0;
  double step_size = 0.0;
  double eta_new = 0.0;
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = step["gradient"];
    hessian = step["hessian"];

    if (abs(hessian) > eps_inv) {

      direction_step = gradient / hessian;

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADCS_EY_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, CS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EY_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, CS_beta, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADCS_EY_normal_rcpp(const arma::mat & X,
                         const double & alpha,
                         const arma::vec & beta,
                         const double & sigma,
                         const double & phi,
                         const arma::vec & CS_beta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p * 2 + 3);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }

    Psi += eSI_CS_i * (SI_i - phi);
    Psi_square += Psi * Psi.t();
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t() * eSI_CS_i;
    Psi_gradient(0, n_p + 2) += -eSI_CS_i;
    Psi_gradient.submat(0, n_p + 3, 0, n_p * 2 + 2) += Xrow_i.t() * eSI_CS_i * (SI_i - phi);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

      Psi = arma::vec(1);
    Psi_square = arma::mat(1, 1);
    Psi_gradient = arma::mat(1, n_p * 2 + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADCSvar_EY_normal_rcpp(const arma::mat & X,
                            const double & alpha,
                            const arma::vec & beta,
                            const double & sigma,
                            const double & phi,
                            const arma::vec & CS_beta,
                            const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  double Psi_i = 0.0;
  double p_i = 0.0;
  double var_A = 0.0;
  double var_B = 0.0;
  double var_B_Ntor = 0.0;
  double var_B_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i = eSI_CS_i * (SI_i - phi);
    p_i = 1 / (1 + Psi_i * eta) / n_n;

    var_B_Ntor += pow(SI_i - phi, 2) * eSI_CS_i * p_i;
    var_B_Dtor += eSI_CS_i * p_i;
  }

  var_A = -1.0;
  if (var_B_Dtor > 0) {

    var_B = var_B_Ntor / var_B_Dtor + pow(sigma, 2);

  } else {

    var_B = var_B_Ntor + pow(sigma, 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List ADCS_EXsubY_normal_Lagrange_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & sigma,
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

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::mat pdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::vec pdf_dist(n_k);
  arma::mat x_phi_dist(n_k, n_p);
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    q_pts = (y_pts - SI_i) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    pdf_pts = arma::normpdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    x_phi_dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    Psi_i = arma::reshape(
      (x_phi_dist % arma::repmat(cdf_dist, 1, n_p) * eSI_CS_i).t(), n_m, 1);

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List ADCS_EXsubY_normal_SolveLagrange_rcpp(const arma::mat & X,
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

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADCS_EXsubY_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, CS_beta, y_pts, eta);
  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADCS_EXsubY_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, CS_beta, y_pts, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EXsubY_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, CS_beta, y_pts, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADCS_EXsubY_normal_rcpp(const arma::mat & X,
                             const double & alpha,
                             const arma::vec & beta,
                             const double & sigma,
                             const arma::mat & phi,
                             const arma::vec & CS_beta,
                             const arma::mat & y_pts) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p * 2 + n_m + 2);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::mat pdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::vec pdf_dist(n_k);
  arma::mat x_phi_dist(n_k, n_p);
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    q_pts = (y_pts - SI_i) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    pdf_pts = arma::normpdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    x_phi_dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    Psi_i = arma::reshape((x_phi_dist % arma::repmat(cdf_dist, 1, n_p) *
      eSI_CS_i).t(), n_m, 1);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
    Psi_gradient.cols(0, n_p) += -arma::reshape((x_phi_dist %
      arma::repmat(pdf_dist / sigma, 1, n_p)).t(), n_m, 1) *
      extXrow_i.t() * eSI_CS_i;
    Psi_gradient.col(n_p + 1) += -arma::reshape((x_phi_dist %
      arma::repmat((q_pts.col(1) / sigma) % pdf_pts.col(1) +
      (q_pts.col(0) / sigma) % pdf_pts.col(0), 1, n_p)).t(), n_m, 1) * eSI_CS_i;
    Psi_gradient.cols(n_p + 2, n_p + n_m + 1).diag() +=
      arma::reshape(arma::repmat(cdf_dist, 1, n_p).t(), n_m, 1) * (-eSI_CS_i);
    Psi_gradient.cols(n_p + n_m + 2, n_p * 2 + n_m + 1) += Psi_i * Xrow_i.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

      Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p * 2 + n_m + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADCSvar_EXsubY_normal_rcpp(const arma::mat & X,
                                const double & alpha,
                                const arma::vec & beta,
                                const double & sigma,
                                const arma::mat & phi,
                                const arma::vec & CS_beta,
                                const arma::mat & y_pts,
                                const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::mat pdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::vec pdf_dist(n_k);
  arma::mat x_phi_dist(n_k, n_p);
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    q_pts = (y_pts - SI_i) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    pdf_pts = arma::normpdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    x_phi_dist = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    Psi_i = arma::reshape((x_phi_dist % arma::repmat(cdf_dist, 1, n_p) *
      eSI_CS_i).t(), n_m, 1);
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    var_A.diag() += -arma::reshape(arma::repmat(
      cdf_dist, 1, n_p).t(), n_m, 1) * eSI_CS_i * p_i;
    for (size_t k = 0; k < n_k; ++k) {

      var_B.submat(k * n_p, k * n_p,
                   (k + 1) * n_p - 1,
                   (k + 1) * n_p - 1) += x_phi_dist.row(k).t() *
                     x_phi_dist.row(k) * cdf_dist(k) * eSI_CS_i * p_i;
    }
    var_AB_Dtor += eSI_CS_i * p_i;
  }

  if (var_AB_Dtor > 0) {

    var_A /= var_AB_Dtor;
    var_B /= var_AB_Dtor;
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List ADCS_EYsubX_normal_Lagrange_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const double & sigma,
                                      const arma::vec & phi,
                                      const arma::vec & CS_beta,
                                      const arma::umat & inclusion,
                                      const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i = (SI_i - phi) % inclusion.row(i).t() * eSI_CS_i;

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List ADCS_EYsubX_normal_SolveLagrange_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const double & sigma,
                                           const arma::vec & phi,
                                           const arma::vec & CS_beta,
                                           const arma::umat & inclusion,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol,
                                           const double & eps_inv) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADCS_EYsubX_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, CS_beta, inclusion, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADCS_EYsubX_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, CS_beta, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EYsubX_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, CS_beta, inclusion, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADCS_EYsubX_normal_rcpp(const arma::mat & X,
                             const double & alpha,
                             const arma::vec & beta,
                             const double & sigma,
                             const arma::vec & phi,
                             const arma::vec & CS_beta,
                             const arma::umat & inclusion) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p * 2 + n_m + 2);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i.zeros();

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) = (SI_i - phi(m)) * eSI_CS_i;
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t() * eSI_CS_i;
        Psi_gradient(m, n_p + 2 + m) += -eSI_CS_i;
        Psi_gradient.submat(m, n_p + n_m + 2, m, n_p * 2 + n_m + 1) += Xrow_i.t() *
          eSI_CS_i * (SI_i - phi(m));
      }
    }

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p * 2 + n_m + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADCSvar_EYsubX_normal_rcpp(const arma::mat & X,
                                const double & alpha,
                                const arma::vec & beta,
                                const double & sigma,
                                const arma::vec & phi,
                                const arma::vec & CS_beta,
                                const arma::umat & inclusion,
                                const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::vec Psi_diff_i(n_m);
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::vec inclusion_i(n_m);
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);
  arma::mat var_B_a(n_m, n_m);
  arma::mat var_B_b(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_diff_i = (SI_i - phi) % inclusion.row(i).t();
    Psi_i = Psi_diff_i * eSI_CS_i;
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        var_A(m, m) += -eSI_CS_i * p_i;
      }
    }
    inclusion_i = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
    var_B_a += inclusion_i * inclusion_i.t() * eSI_CS_i * p_i;
    var_B_b += Psi_diff_i * Psi_diff_i.t() * eSI_CS_i * p_i;
    var_AB_Dtor += eSI_CS_i * p_i;
  }

  var_B = var_B_a * pow(sigma, 2) + var_B_b;

  if (var_AB_Dtor > 0) {

    var_A /= var_AB_Dtor;
    var_B /= var_AB_Dtor;

  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// Prior probability shift, normal regression model

// [[Rcpp::export]]
List ADPPS_EX_normal_Lagrange_rcpp(const arma::mat & X,
                                   const double & alpha,
                                   const arma::vec & beta,
                                   const double & sigma,
                                   const arma::vec & phi,
                                   const double & PPS_beta,
                                   const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  arma::vec gradient(n_p);
  arma::mat hessian(n_p, n_p);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double e_f = 0.0;
  arma::vec Psi_i(n_p);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_i = (Xrow_i - phi) * e_f;

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List ADPPS_EX_normal_SolveLagrange_rcpp(const arma::mat & X,
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

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADPPS_EX_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, PPS_beta, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADPPS_EX_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, PPS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EX_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, PPS_beta, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADPPS_EX_normal_rcpp(const arma::mat & X,
                          const double & alpha,
                          const arma::vec & beta,
                          const double & sigma,
                          const arma::vec & phi,
                          const double & PPS_beta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Psi(n_p);
  arma::mat Psi_square(n_p, n_p);
  arma::mat Psi_gradient(n_p, n_p * 2 + 3);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double e_f = 0.0;
  arma::vec Psi_i(n_p);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_i = (Xrow_i - phi) * e_f;

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
    Psi_gradient.cols(0, n_p) += Psi_i * extXrow_i.t() * PPS_beta;
    Psi_gradient.col(n_p + 1) += Psi_i * pow(PPS_beta, 2) * sigma;
    Psi_gradient.cols(n_p + 2, n_p * 2 + 1).diag() += -e_f;
    Psi_gradient.col(n_p * 2 + 2) += Psi_i * pow(sigma, 2) * PPS_beta;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_p);
    Psi_square = arma::mat(n_p, n_p);
    Psi_gradient = arma::mat(n_p, n_p * 2 + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADPPSvar_EX_normal_rcpp(const arma::mat & X,
                             const double & alpha,
                             const arma::vec & beta,
                             const double & sigma,
                             const arma::vec & phi,
                             const double & PPS_beta,
                             const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double e_f = 0.0;
  arma::vec Psi_i(n_p);
  double p_i = 0.0;
  arma::mat var_A(n_p, n_p);
  arma::mat var_B(n_p, n_p);
  double var_B_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_i = (Xrow_i - phi) * e_f;
    p_i = 1 / (1 + arma::dot(Psi_i, eta)) / n_n;

    var_B += (Xrow_i - phi) * (Xrow_i - phi).t() * e_f * p_i;
    var_B_Dtor += e_f * p_i;
  }

  var_A.diag().fill(-1);
  if (var_B_Dtor < 0) {

    var_B /= var_B_Dtor;
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List ADPPS_EY_normal_Lagrange_rcpp(const arma::mat & X,
                                   const double & alpha,
                                   const arma::vec & beta,
                                   const double & sigma,
                                   const double & phi,
                                   const double & PPS_beta,
                                   const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double c_f = 0.0;
  double e_f = 0.0;
  double Psi_i = 0.0;
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    c_f = SI_i + pow(sigma, 2) * PPS_beta - phi;
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_i = c_f * e_f;

    denominator = 1.0 + eta * Psi_i;
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
List ADPPS_EY_normal_SolveLagrange_rcpp(const arma::mat & X,
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

  double eta = eta_initial;
  List step = ADPPS_EY_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, PPS_beta, eta);

  double step_value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  double direction_step = 0.0;
  double step_size = 0.0;
  double eta_new = 0.0;
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = step["gradient"];
    hessian = step["hessian"];

    if (abs(hessian) > eps_inv) {

      direction_step = gradient / hessian;

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADPPS_EY_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, PPS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EY_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, PPS_beta, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADPPS_EY_normal_rcpp(const arma::mat & X,
                          const double & alpha,
                          const arma::vec & beta,
                          const double & sigma,
                          const double & phi,
                          const double & PPS_beta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p + 4);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double c_f = 0.0;
  double e_f = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    c_f = SI_i + pow(sigma, 2) * PPS_beta - phi;
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);

    Psi += c_f * e_f;
    Psi_square += Psi * Psi.t();
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t() * (c_f * PPS_beta + 1) * e_f;
    Psi_gradient(0, n_p + 1) += (c_f * PPS_beta + 2) * sigma * PPS_beta * e_f;
    Psi_gradient(0, n_p + 2) += -e_f;
    Psi_gradient(0, n_p + 3) += (c_f * PPS_beta + 1) * pow(sigma, 2) * e_f;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(1);
    Psi_square = arma::mat(1, 1);
    Psi_gradient = arma::mat(1, n_p + 4);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADPPSvar_EY_normal_rcpp(const arma::mat & X,
                             const double & alpha,
                             const arma::vec & beta,
                             const double & sigma,
                             const double & phi,
                             const double & PPS_beta,
                             const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double c_f = 0.0;
  double e_f = 0.0;
  double Psi_i = 0.0;
  double p_i = 0.0;
  double var_A = 0.0;
  double var_B = 0.0;
  double var_B_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    c_f = SI_i + pow(sigma, 2) * PPS_beta - phi;
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_i = c_f * e_f;
    p_i = 1 / (1 + Psi_i * eta) / n_n;
    var_B += pow(c_f, 2) * e_f * p_i;
    var_B_Dtor += e_f * p_i;
  }

  var_A = -1.0;
  if (var_B_Dtor > 0) {

    var_B /=var_B_Dtor;
  }
  var_B += pow(sigma, 2);

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List ADPPS_EXsubY_normal_Lagrange_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & sigma,
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

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double e_f = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::mat dist_e_f(n_k, n_p);
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    q_pts = (y_pts - SI_i - pow(sigma, 2) * PPS_beta) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    dist_e_f = (arma::repmat(Xrow_i.t(), n_k, 1) - phi) * e_f;
    Psi_i = arma::reshape(
      (dist_e_f % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List ADPPS_EXsubY_normal_SolveLagrange_rcpp(const arma::mat & X,
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

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADPPS_EXsubY_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, PPS_beta, y_pts, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADPPS_EXsubY_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, PPS_beta, y_pts, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EXsubY_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, PPS_beta, y_pts, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADPPS_EXsubY_normal_rcpp(const arma::mat & X,
                              const double & alpha,
                              const arma::vec & beta,
                              const double & sigma,
                              const arma::mat & phi,
                              const double & PPS_beta,
                              const arma::mat & y_pts) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 3);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double e_f = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::mat pdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::vec pdf_dist(n_k);
  arma::mat dist_e_f(n_k, n_p);
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    q_pts = (y_pts - SI_i - pow(sigma, 2) * PPS_beta) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    pdf_pts = arma::normpdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    dist_e_f = (arma::repmat(Xrow_i.t(), n_k, 1) - phi) * e_f;
    Psi_i = arma::reshape((dist_e_f % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
    Psi_gradient.cols(0, n_p) += arma::reshape((dist_e_f %
      arma::repmat(cdf_dist * PPS_beta - pdf_dist / sigma, 1, n_p)).t(), n_m, 1) *
      extXrow_i.t();
    Psi_gradient.col(n_p + 1) += arma::reshape((dist_e_f %
      arma::repmat(cdf_dist * sigma * pow(PPS_beta, 2) -
      (q_pts.col(1) / sigma + PPS_beta * 2) % pdf_pts.col(1) +
      (q_pts.col(0) / sigma + PPS_beta * 2) % pdf_pts.col(0), 1, n_p)).t(), n_m, 1);
    Psi_gradient.cols(n_p + 2, n_p + n_m + 1).diag() +=
      arma::reshape(arma::repmat(cdf_dist, 1, n_p).t(), n_m, 1) * (-e_f);
    Psi_gradient.col(n_p + n_m + 2) +=
      arma::reshape((dist_e_f % arma::repmat(cdf_dist * PPS_beta * pow(sigma, 2) -
      pdf_dist * sigma, 1, n_p)).t(), n_m, 1);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p + n_m + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADPPSvar_EXsubY_normal_rcpp(const arma::mat & X,
                                 const double & alpha,
                                 const arma::vec & beta,
                                 const double & sigma,
                                 const arma::mat & phi,
                                 const double & PPS_beta,
                                 const arma::mat & y_pts,
                                 const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = phi.n_rows;
  const arma::uword n_m = n_p * n_k;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  double e_f = 0.0;
  arma::mat q_pts(n_k, 2);
  arma::mat cdf_pts(n_k, 2);
  arma::mat pdf_pts(n_k, 2);
  arma::vec cdf_dist(n_k);
  arma::vec pdf_dist(n_k);
  arma::mat dist_x_phi(n_k, n_p);
  arma::mat dist_e_f(n_k, n_p);
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    q_pts = (y_pts - SI_i - pow(sigma, 2) * PPS_beta) / sigma;
    cdf_pts = arma::normcdf(q_pts);
    pdf_pts = arma::normpdf(q_pts);
    cdf_dist = cdf_pts.col(1) - cdf_pts.col(0);
    pdf_dist = pdf_pts.col(1) - pdf_pts.col(0);
    dist_x_phi = (arma::repmat(Xrow_i.t(), n_k, 1) - phi);
    dist_e_f = dist_x_phi * e_f;
    Psi_i = arma::reshape(
      (dist_e_f % arma::repmat(cdf_dist, 1, n_p)).t(), n_m, 1);
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    var_A.diag() += -arma::reshape(arma::repmat(
      cdf_dist, 1, n_p).t(), n_m, 1) * e_f * p_i;
    for (size_t k = 0; k < n_k; ++k) {

      var_B.submat(k * n_p, k * n_p,
                   (k + 1) * n_p - 1,
                   (k + 1) * n_p - 1) += dist_x_phi.row(k).t() *
                     dist_x_phi.row(k) * cdf_dist(k) * e_f * p_i;
    }
    var_AB_Dtor += e_f * p_i;
  }

  if (var_AB_Dtor > 0) {

    var_A /= var_AB_Dtor;
    var_B /= var_AB_Dtor;
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List ADPPS_EYsubX_normal_Lagrange_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & sigma,
                                       const arma::vec & phi,
                                       const double & PPS_beta,
                                       const arma::umat & inclusion,
                                       const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::vec c_f(n_m);
  double e_f = 0.0;
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    c_f = SI_i + pow(sigma, 2) * PPS_beta - phi;
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_i = c_f % inclusion.row(i).t() * e_f;

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List ADPPS_EYsubX_normal_SolveLagrange_rcpp(const arma::mat & X,
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

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADPPS_EYsubX_normal_Lagrange_rcpp(
    X, alpha, beta, sigma, phi, PPS_beta, inclusion, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADPPS_EYsubX_normal_Lagrange_rcpp(
      X, alpha, beta, sigma, phi, PPS_beta, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EYsubX_normal_Lagrange_rcpp(
          X, alpha, beta, sigma, phi, PPS_beta, inclusion, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADPPS_EYsubX_normal_rcpp(const arma::mat & X,
                              const double & alpha,
                              const arma::vec & beta,
                              const double & sigma,
                              const arma::vec & phi,
                              const double & PPS_beta,
                              const arma::umat & inclusion) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 3);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::vec c_f(n_m);
  double e_f = 0.0;
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    c_f = SI_i + pow(sigma, 2) * PPS_beta - phi;
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_i.zeros();

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) = c_f(m) * e_f;
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t() *
          (c_f(m) * PPS_beta + 1) * e_f;
        Psi_gradient(m, n_p + 1) += (c_f(m) * PPS_beta + 2) *
          sigma * PPS_beta * e_f;
        Psi_gradient(m, n_p + 2 + m) += -e_f;
        Psi_gradient(m, n_p + n_m + 2) += (c_f(m) * PPS_beta + 1) *
          pow(sigma, 2) * e_f;
      }
    }

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p + n_m + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADPPSvar_EYsubX_normal_rcpp(const arma::mat & X,
                                 const double & alpha,
                                 const arma::vec & beta,
                                 const double & sigma,
                                 const arma::vec & phi,
                                 const double & PPS_beta,
                                 const arma::umat & inclusion,
                                 const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double SI_i = 0.0;
  arma::vec c_f(n_m);
  double e_f = 0.0;
  arma::vec Psi_diff_i(n_m);
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::vec inclusion_i(n_m);
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);
  arma::mat var_B_a(n_m, n_m);
  arma::mat var_B_b(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    SI_i = alpha + arma::dot(Xrow_i, beta);
    c_f = SI_i + pow(sigma, 2) * PPS_beta - phi;
    e_f = exp((SI_i * 2 + pow(sigma, 2) * PPS_beta) * PPS_beta / 2.0);
    Psi_diff_i = c_f % inclusion.row(i).t();
    Psi_i = Psi_diff_i * e_f;
    p_i = 1 / (1 + arma::dot(Psi_i, eta)) / n_n;

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        var_A(m, m) += -e_f * p_i;
      }
    }
    inclusion_i = arma::conv_to<arma::vec>::from(inclusion.row(i).t());
    var_B_a += inclusion_i * inclusion_i.t() * e_f * p_i;
    var_B_b += Psi_diff_i * Psi_diff_i.t() * e_f * p_i;
    var_AB_Dtor += e_f * p_i;
  }

  var_B = var_B_a * pow(sigma, 2) + var_B_b;

  if (var_AB_Dtor > 0) {

    var_A /= var_AB_Dtor;
    var_B /= var_AB_Dtor;
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// No shift, logistic regression model

// [[Rcpp::export]]
List AD_EY_logistic_Lagrange_rcpp(const arma::mat & X,
                                  const double & alpha,
                                  const arma::vec & beta,
                                  const double & phi,
                                  const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double Psi_i = 0.0;
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i = ((1 - phi) * eSI_i - phi) / (1 + eSI_i);

    denominator = 1.0 + eta * Psi_i;
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
List AD_EY_logistic_SolveLagrange_rcpp(const arma::mat & X,
                                       const double & alpha,
                                       const arma::vec & beta,
                                       const double & phi,
                                       const double & eta_initial,
                                       const size_t & iter_max,
                                       const double & step_rate,
                                       const size_t & step_max,
                                       const double & tol,
                                       const double & eps_inv) {

  double eta = eta_initial;
  List step = AD_EY_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, eta);

  double step_value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  double direction_step = 0.0;
  double step_size = 0.0;
  double eta_new = 0.0;
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = step["gradient"];
    hessian = step["hessian"];

    if (abs(hessian) > eps_inv) {

      direction_step = gradient / hessian;

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = AD_EY_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EY_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List AD_EY_logistic_rcpp(const arma::mat & X,
                         const double & alpha,
                         const arma::vec & beta,
                         const double & phi) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p + 2);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double Psi_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i = ((1 - phi) * eSI_i - phi) / (1 + eSI_i);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i;
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t() * eSI_i / pow(1 + eSI_i, 2);
    Psi_gradient(0, n_p + 1) += -1;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(1);
    Psi_square = arma::mat(1, 1);
    Psi_gradient = arma::mat(1, n_p + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADvar_EY_logistic_rcpp(const arma::mat & X,
                            const double & alpha,
                            const arma::vec & beta,
                            const double & phi,
                            const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double Psi_i = 0.0;
  double var_A = 0.0;
  double var_B = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i = ((1 - phi) * eSI_i - phi) / (1 + eSI_i);
    var_B += (pow(1 - phi, 2) * eSI_i + pow(phi, 2)) / (1 + eSI_i)
      / n_n / (1.0 + Psi_i * eta);
  }

  var_A = -1.0;

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List AD_EXsubY_logistic_Lagrange_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const arma::mat & phi,
                                      const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = n_p * 2;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  arma::mat Psi_i_mat(2, n_p);
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i_mat.row(0) = (Xrow_i.t() - phi.row(0)) * eSI_i / (1 + eSI_i);
    Psi_i_mat.row(1) = (Xrow_i.t() - phi.row(1)) / (1 + eSI_i);
    Psi_i = arma::reshape(Psi_i_mat.t(), n_m, 1);

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List AD_EXsubY_logistic_SolveLagrange_rcpp(const arma::mat & X,
                                           const double & alpha,
                                           const arma::vec & beta,
                                           const arma::mat & phi,
                                           const arma::vec & eta_initial,
                                           const size_t & iter_max,
                                           const double & step_rate,
                                           const size_t & step_max,
                                           const double & tol,
                                           const double & eps_inv) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = AD_EXsubY_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = AD_EXsubY_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EXsubY_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List AD_EXsubY_logistic_rcpp(const arma::mat & X,
                             const double & alpha,
                             const arma::vec & beta,
                             const arma::mat & phi) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = n_p * 2;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 1);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  arma::mat Psi_i_mat(2, n_p);
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i_mat.row(0) = (Xrow_i.t() - phi.row(0)) * eSI_i / (1 + eSI_i);
    Psi_i_mat.row(1) = (Xrow_i.t() - phi.row(1)) / (1 + eSI_i);
    Psi_i = arma::reshape(Psi_i_mat.t(), n_m, 1);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
    Psi_i_mat.row(0) /= (1 + eSI_i);
    Psi_i_mat.row(1) *= -eSI_i / (1 + eSI_i);
    Psi_gradient.cols(0, n_p) += arma::reshape(Psi_i_mat.t(), n_m, 1) *
      extXrow_i.t();
    Psi_gradient.submat(0, n_p + 1, n_p - 1, n_p * 2).diag() += -eSI_i / (1 + eSI_i);
    Psi_gradient.submat(n_p, n_p * 2 + 1, n_p * 2 - 1, n_p * 3).diag() += -1 / (1 + eSI_i);
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p + n_m + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADvar_EXsubY_logistic_rcpp(const arma::mat & X,
                                const double & alpha,
                                const arma::vec & beta,
                                const arma::mat & phi,
                                const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = n_p * 2;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  arma::mat Psi_i_mat(2, n_p);
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i_mat.row(0) = (Xrow_i.t() - phi.row(0)) * eSI_i / (1 + eSI_i);
    Psi_i_mat.row(1) = (Xrow_i.t() - phi.row(1)) / (1 + eSI_i);
    Psi_i = arma::reshape(Psi_i_mat.t(), n_m, 1);
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    var_A.submat(0, 0, n_p - 1, n_p - 1).diag() += -eSI_i / (1 + eSI_i) * p_i;
    var_A.submat(n_p, n_p, n_p * 2 - 1, n_p * 2 - 1).diag() += -1 / (1 + eSI_i) * p_i;
    var_B.submat(0, 0, n_p - 1, n_p - 1) += (Xrow_i.t() - phi.row(0)).t() *
      (Xrow_i.t() - phi.row(0)) * eSI_i / (1 + eSI_i) * p_i;
    var_B.submat(n_p, n_p, n_p * 2 - 1, n_p * 2 - 1) += (Xrow_i.t() - phi.row(1)).t() *
      (Xrow_i.t() - phi.row(1)) / (1 + eSI_i) * p_i;
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List AD_EYsubX_logistic_Lagrange_rcpp(const arma::mat & X,
                                      const double & alpha,
                                      const arma::vec & beta,
                                      const arma::vec & phi,
                                      const arma::umat & inclusion,
                                      const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i.zeros();

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += ((1 - phi(m)) * eSI_i - phi(m)) / (1 + eSI_i);
      }
    }

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List AD_EYsubX_logistic_SolveLagrange_rcpp(const arma::mat & X,
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

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = AD_EYsubX_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, inclusion, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = AD_EYsubX_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EYsubX_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, inclusion, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List AD_EYsubX_logistic_rcpp(const arma::mat & X,
                             const double & alpha,
                             const arma::vec & beta,
                             const arma::vec & phi,
                             const arma::umat & inclusion) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p + n_m + 1);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i.zeros();

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += ((1 - phi(m)) * eSI_i - phi(m)) / (1 + eSI_i);
        Psi_gradient.submat(m, 0, m, n_p) += extXrow_i.t() * eSI_i / pow(1 + eSI_i, 2);
        Psi_gradient(m, n_p + 1 + m) += -1;
      }
    }

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p + n_m + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADvar_EYsubX_logistic_rcpp(const arma::mat & X,
                                const double & alpha,
                                const arma::vec & beta,
                                const arma::vec & phi,
                                const arma::umat & inclusion,
                                const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = phi.n_elem;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);
  double var_B_off_diag = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    Psi_i.zeros();
    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += ((1 - phi(m)) * eSI_i - phi(m)) / (1 + eSI_i);
      }
    }
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        var_A(m, m) += -p_i;
        var_B(m, m) += (pow(1 - phi(m), 2) * eSI_i +
          pow(phi(m), 2)) / (1 + eSI_i) * p_i;

        for (size_t l = 0; l < m; ++l) {

          if (inclusion(i, l) == 1) {

            var_B_off_diag = ((1 - phi(l)) * (1 - phi(m)) * eSI_i +
              phi(l) * phi(m)) / (1 + eSI_i) * p_i;
            var_B(l, m) += var_B_off_diag;
            var_B(m, l) += var_B_off_diag;
          }
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// Covariate shift, logistic regression model

// [[Rcpp::export]]
List ADCS_EY_logistic_Lagrange_rcpp(const arma::mat & X,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & phi,
                                    const arma::vec & CS_beta,
                                    const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double eSI_CS_i = 0.0;
  double Psi_i = 0.0;
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i = ((1 - phi) * eSI_i - phi) * eSI_CS_i / (1 + eSI_i);

    denominator = 1.0 + eta * Psi_i;
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
List ADCS_EY_logistic_SolveLagrange_rcpp(const arma::mat & X,
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

  double eta = eta_initial;
  List step = ADCS_EY_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, CS_beta, eta);

  double step_value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;
  double direction_step = 0.0;
  double step_size = 0.0;
  double eta_new = 0.0;
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = step["gradient"];
    hessian = step["hessian"];

    if (abs(hessian) > eps_inv) {

      direction_step = gradient / hessian;

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    eta_new = eta - direction_step;
    List step_new = ADCS_EY_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, CS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EY_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, CS_beta, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADCS_EY_logistic_rcpp(const arma::mat & X,
                           const double & alpha,
                           const arma::vec & beta,
                           const double & phi,
                           const arma::vec & CS_beta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p * 2 + 2);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double eSI_CS_i = 0.0;
  double Psi_i = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i = ((1 - phi) * eSI_i - phi) * eSI_CS_i / (1 + eSI_i);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i;
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t() * eSI_i /
      pow(1 + eSI_i, 2) * eSI_CS_i;
    Psi_gradient(0, n_p + 1) += -eSI_CS_i;
    Psi_gradient.submat(0, n_p + 2, 0, n_p * 2 + 1) += Xrow_i.t() * Psi_i;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(1);
    Psi_square = arma::mat(1, 1);
    Psi_gradient = arma::mat(1, n_p + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADCSvar_EY_logistic_rcpp(const arma::mat & X,
                              const double & alpha,
                              const arma::vec & beta,
                              const double & phi,
                              const arma::vec & CS_beta,
                              const double & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double eSI_CS_i = 0.0;
  double Psi_i = 0.0;
  double p_i = 0.0;
  double var_A = 0.0;
  double var_B = 0.0;
  double var_B_Ntor = 0.0;
  double var_B_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i = ((1 - phi) * eSI_i - phi) * eSI_CS_i / (1 + eSI_i);
    p_i = 1 / (1.0 + Psi_i * eta) / n_n;
    var_B_Ntor += (pow(1 - phi, 2) * eSI_i + pow(phi, 2)) / (1 + eSI_i) *
      eSI_CS_i * p_i;
    var_B_Dtor += eSI_CS_i * p_i;
  }

  var_A = -1.0;
  if (var_B_Dtor > 0) {

    var_B = var_B_Ntor / var_B_Dtor;

  } else {

    var_B = var_B_Ntor;
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

// [[Rcpp::export]]
List ADCS_EXsubY_logistic_Lagrange_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const arma::mat & phi,
                                        const arma::vec & CS_beta,
                                        const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = n_p * 2;
  double value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::mat Psi_i_mat(2, n_p);
  arma::vec Psi_i(n_m);
  double denominator = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i_mat.row(0) = (Xrow_i.t() - phi.row(0)) * eSI_i * eSI_CS_i / (1 + eSI_i);
    Psi_i_mat.row(1) = (Xrow_i.t() - phi.row(1)) * eSI_CS_i / (1 + eSI_i);
    Psi_i = arma::reshape(Psi_i_mat.t(), n_m, 1);

    denominator = 1.0 + arma::dot(eta, Psi_i);
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
List ADCS_EXsubY_logistic_SolveLagrange_rcpp(const arma::mat & X,
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

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADCS_EXsubY_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, CS_beta, eta);

  double step_value = 0.0;
  arma::vec gradient(n_m);
  arma::mat hessian(n_m, n_m);
  arma::vec direction_step(n_m);
  double step_size = 0.0;
  arma::vec eta_new(n_m);
  double step_new_value = 0.0;

  for (size_t k = 0; k < iter_max; k++) {

    step_value = step["value"];
    gradient = as<arma::vec>(step["gradient"]);
    hessian = as<arma::mat>(step["hessian"]);

    if (arma::rcond(hessian) > eps_inv) {

      direction_step = arma::solve(hessian, gradient);

    } else {

      direction_step = -gradient;
    }

    step_size = 1.0;
    arma::vec eta_new = eta - direction_step;
    List step_new = ADCS_EXsubY_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, CS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EXsubY_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, CS_beta, eta_new);

      } else {

        break;
      }
    }

    step_new_value = step_new["value"];

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
List ADCS_EXsubY_logistic_rcpp(const arma::mat & X,
                               const double & alpha,
                               const arma::vec & beta,
                               const arma::mat & phi,
                               const arma::vec & CS_beta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = n_p * 2;
  arma::vec Psi(n_m);
  arma::mat Psi_square(n_m, n_m);
  arma::mat Psi_gradient(n_m, n_p * 2 + n_m + 1);

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::mat Psi_i_mat(2, n_p);
  arma::vec Psi_i(n_m);

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i_mat.row(0) = (Xrow_i.t() - phi.row(0)) * eSI_i / (1 + eSI_i) * eSI_CS_i;
    Psi_i_mat.row(1) = (Xrow_i.t() - phi.row(1)) / (1 + eSI_i) * eSI_CS_i;
    Psi_i = arma::reshape(Psi_i_mat.t(), n_m, 1);

    Psi += Psi_i;
    Psi_square += Psi_i * Psi_i.t();
    Psi_i_mat.row(0) /= (1 + eSI_i);
    Psi_i_mat.row(1) *= -eSI_i / (1 + eSI_i);
    Psi_gradient.cols(0, n_p) += arma::reshape(Psi_i_mat.t(), n_m, 1) *
      extXrow_i.t() * eSI_CS_i;
    Psi_gradient.submat(0, n_p + 1, n_p - 1, n_p * 2).diag() += -eSI_i / (1 + eSI_i) * eSI_CS_i;
    Psi_gradient.submat(n_p, n_p * 2 + 1, n_p * 2 - 1, n_p * 3).diag() += -1 / (1 + eSI_i) * eSI_CS_i;
    Psi_gradient.cols(n_p * 3 + 1, n_p * 4) += Psi_i * Xrow_i.t();
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(n_m);
    Psi_square = arma::mat(n_m, n_m);
    Psi_gradient = arma::mat(n_m, n_p + n_m + 2);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
}

// [[Rcpp::export]]
List ADCSvar_EXsubY_logistic_rcpp(const arma::mat & X,
                                  const double & alpha,
                                  const arma::vec & beta,
                                  const arma::mat & phi,
                                  const arma::vec & CS_beta,
                                  const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_m = n_p * 2;

  const arma::vec X_intercept = arma::vec({1});

  arma::vec Xrow_i(n_p);
  arma::vec extXrow_i(n_p + 1);
  double eSI_i = 0.0;
  double eSI_CS_i = 0.0;
  arma::mat Psi_i_mat(2, n_p);
  arma::vec Psi_i(n_m);
  double p_i = 0.0;
  arma::mat var_A(n_m, n_m);
  arma::mat var_B(n_m, n_m);
  double var_AB_Dtor = 0.0;

  for (size_t i = 0; i < n_n; ++i) {

    Xrow_i = X.row(i).t();
    extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    if (eSI_CS_i > 1e10) {

      eSI_CS_i = 1e10;
    }
    Psi_i_mat.row(0) = (Xrow_i.t() - phi.row(0)) * eSI_i / (1 + eSI_i) * eSI_CS_i;
    Psi_i_mat.row(1) = (Xrow_i.t() - phi.row(1)) / (1 + eSI_i) * eSI_CS_i;
    Psi_i = arma::reshape(Psi_i_mat.t(), n_m, 1);
    p_i = 1.0 / (1.0 + arma::dot(Psi_i, eta)) / n_n;

    var_A.submat(0, 0, n_p - 1, n_p - 1).diag() += -eSI_i / (1 + eSI_i) * eSI_CS_i * p_i;
    var_A.submat(n_p, n_p, n_p * 2 - 1, n_p * 2 - 1).diag() += -1 / (1 + eSI_i) * eSI_CS_i * p_i;
    var_B.submat(0, 0, n_p - 1, n_p - 1) += (Xrow_i.t() - phi.row(0)).t() *
      (Xrow_i.t() - phi.row(0)) * eSI_i / (1 + eSI_i) * eSI_CS_i * p_i;
    var_B.submat(n_p, n_p, n_p * 2 - 1, n_p * 2 - 1) += (Xrow_i.t() - phi.row(1)).t() *
      (Xrow_i.t() - phi.row(1)) / (1 + eSI_i) * eSI_CS_i * p_i;
    var_AB_Dtor += eSI_CS_i * p_i;
  }

  if (var_AB_Dtor > 0) {

    var_A /= var_AB_Dtor;
    var_B /= var_AB_Dtor;
  }

  return Rcpp::List::create(
    Rcpp::Named("diff") = var_A,
    Rcpp::Named("var") = var_B
  );
}

































// [[Rcpp::export]]
List ADCS_EYsubX_logistic_Lagrange_rcpp(const arma::mat & X,
                                        const double & alpha,
                                        const arma::vec & beta,
                                        const arma::vec & phi,
                                        const arma::vec & CS_beta,
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
    const double eSI_CS_i = exp(arma::dot(Xrow_i, CS_beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += ((1 - phi(m)) * eSI_i * eSI_CS_i - phi(m)) / (1 + eSI_i);
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
List ADCS_EYsubX_logistic_SolveLagrange_rcpp(const arma::mat & X,
                                             const double & alpha,
                                             const arma::vec & beta,
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
  List step = ADCS_EYsubX_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, CS_beta, inclusion, eta);

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
    List step_new = ADCS_EYsubX_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, CS_beta, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADCS_EYsubX_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, CS_beta, inclusion, eta_new);

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

// Prior probability shift, logistic regression model

// [[Rcpp::export]]
List ADPPS_EX_logistic_Lagrange_rcpp(const arma::mat & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const arma::vec & phi,
                                     const double & PPS_beta,
                                     const arma::vec & eta) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  double value = 0.0;
  arma::vec gradient(n_p);
  arma::mat hessian(n_p, n_p);

  const arma::vec X_intercept = arma::vec({1});
  const double ePPS_beta = exp(PPS_beta);

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const arma::vec Psi_i = (Xrow_i - phi) *
      (1 + eSI_i * ePPS_beta) / (1 + eSI_i);

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
List ADPPS_EX_logistic_SolveLagrange_rcpp(const arma::mat & X,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const arma::vec & phi,
                                          const double & PPS_beta,
                                          const arma::vec & eta_initial,
                                          const size_t & iter_max,
                                          const double & step_rate,
                                          const size_t & step_max,
                                          const double & tol) {

  const arma::uword n_m = eta_initial.n_elem;
  arma::vec eta = eta_initial;
  List step = ADPPS_EX_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, PPS_beta, eta);

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
    List step_new = ADPPS_EX_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, PPS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EX_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, PPS_beta, eta_new);

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
List ADPPS_EY_logistic_Lagrange_rcpp(const arma::mat & X,
                                     const double & alpha,
                                     const arma::vec & beta,
                                     const double & phi,
                                     const double & PPS_beta,
                                     const double & eta) {

  const arma::uword n_n = X.n_rows;
  double value = 0.0;
  double gradient = 0.0;
  double hessian = 0.0;

  const arma::vec X_intercept = arma::vec({1});
  const double ePPS_beta = exp(PPS_beta);

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double Psi_i = ((1 - phi) * eSI_i * ePPS_beta - phi) / (1 + eSI_i);

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
List ADPPS_EY_logistic_SolveLagrange_rcpp(const arma::mat & X,
                                          const double & alpha,
                                          const arma::vec & beta,
                                          const double & phi,
                                          const double & PPS_beta,
                                          const double & eta_initial,
                                          const size_t & iter_max,
                                          const double & step_rate,
                                          const size_t & step_max,
                                          const double & tol) {

  double eta = eta_initial;
  List step = ADPPS_EY_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, PPS_beta, eta);

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
    List step_new = ADPPS_EY_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, PPS_beta, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EY_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, PPS_beta, eta_new);

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
List ADPPS_EYsubX_logistic_Lagrange_rcpp(const arma::mat & X,
                                         const double & alpha,
                                         const arma::vec & beta,
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
  const double ePPS_beta = exp(PPS_beta);

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    arma::vec Psi_i(n_m);

    for (size_t m = 0; m < n_m; ++m) {

      if (inclusion(i, m) == 1) {

        Psi_i(m) += ((1 - phi(m)) * eSI_i * ePPS_beta - phi(m)) / (1 + eSI_i);
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
List ADPPS_EYsubX_logistic_SolveLagrange_rcpp(const arma::mat & X,
                                              const double & alpha,
                                              const arma::vec & beta,
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
  List step = ADPPS_EYsubX_logistic_Lagrange_rcpp(
    X, alpha, beta, phi, PPS_beta, inclusion, eta);

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
    List step_new = ADPPS_EYsubX_logistic_Lagrange_rcpp(
      X, alpha, beta, phi, PPS_beta, inclusion, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = ADPPS_EYsubX_logistic_Lagrange_rcpp(
          X, alpha, beta, phi, PPS_beta, inclusion, eta_new);

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










// No shift, Gamma regression model

// [[Rcpp::export]]
List AD_EY_Gamma_Lagrange_rcpp(const arma::mat & X,
                               const double & alpha,
                               const arma::vec & beta,
                               const double & nu,
                               const double & phi,
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
    const double Psi_i = eSI_i - phi;

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
List AD_EY_Gamma_SolveLagrange_rcpp(const arma::mat & X,
                                    const double & alpha,
                                    const arma::vec & beta,
                                    const double & nu,
                                    const double & phi,
                                    const double & eta_initial,
                                    const size_t & iter_max,
                                    const double & step_rate,
                                    const size_t & step_max,
                                    const double & tol) {

  double eta = eta_initial;
  List step = AD_EY_Gamma_Lagrange_rcpp(
    X, alpha, beta, nu, phi, eta);

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
    List step_new = AD_EY_Gamma_Lagrange_rcpp(
      X, alpha, beta, nu, phi, eta_new);

    for (size_t iter_step = 0; iter_step < step_max; iter_step++) {

      const double step_new_value = step_new["value"];

      if (step_new_value <= step_value + tol) {

        step_size /= step_rate;
        eta_new = eta - direction_step * step_size;
        step_new = AD_EY_Gamma_Lagrange_rcpp(
          X, alpha, beta, nu, phi, eta_new);

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
List AD_EY_Gamma_rcpp(const arma::mat & X,
                      const double & alpha,
                      const arma::vec & beta,
                      const double & nu,
                      const double & phi) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Psi(1);
  arma::mat Psi_square(1, 1);
  arma::mat Psi_gradient(1, n_p + 3);

  const arma::vec X_intercept = arma::vec({1});

  for (size_t i = 0; i < n_n; ++i) {

    const arma::vec Xrow_i = X.row(i).t();
    const arma::vec extXrow_i = arma::join_vert(X_intercept, Xrow_i);
    const double eSI_i = exp(alpha + arma::dot(Xrow_i, beta));
    const double c_f = eSI_i - phi;

    Psi += c_f;
    Psi_square += Psi * Psi.t();
    Psi_gradient.submat(0, 0, 0, n_p) += extXrow_i.t() * eSI_i;
    Psi_gradient(0, n_p + 2) += -1;
  }

  Psi /= n_n;
  Psi_square /= n_n;
  Psi_gradient /= n_n;

  if (!Psi.is_finite() ||
      !Psi_square.is_finite() ||
      !Psi_gradient.is_finite()) {

    Psi = arma::vec(1);
    Psi_square = arma::mat(1, 1);
    Psi_gradient = arma::mat(1, n_p + 3);
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = Psi,
    Rcpp::Named("score_square") = Psi_square,
    Rcpp::Named("score_gradient") = Psi_gradient
  );
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





































