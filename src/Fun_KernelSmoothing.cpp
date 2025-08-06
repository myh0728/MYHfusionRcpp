#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// kernel functions

// [[Rcpp::export]]
double K2_Ep_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = 0.75 * (1.0 - pow(u, 2));
  }

  return K_val;
}

// [[Rcpp::export]]
double K2_Ep_diff1_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = -1.5 * u;
  }

  return K_val;
}

// [[Rcpp::export]]
double K2_Ep_diff2_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = -1.5;
  }

  return K_val;
}

// [[Rcpp::export]]
double K2_Bw_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = 0.9375 * pow((1.0 - pow(u, 2)), 2);
  }

  return K_val;
}

// [[Rcpp::export]]
double K2_Bw_diff1_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = -3.75 * u * (1.0 - pow(u, 2));
  }

  return K_val;
}

// [[Rcpp::export]]
double K2_Bw_diff2_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = -3.75 * (1.0 - 3.0 * pow(u, 2));
  }

  return K_val;
}

// [[Rcpp::export]]
double K4_Bw_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = 1.640625 * (1.0 - 3.0 * pow(u, 2)) * pow((1.0 - pow(u, 2)), 2);
  }

  return K_val;
}

// [[Rcpp::export]]
double K4_Bw_diff1_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = 3.28125 * u * (1.0 - pow(u, 2)) * (9.0 * pow(u, 2) - 5.0);
  }

  return K_val;
}

// [[Rcpp::export]]
double K4_Bw_diff2_rcpp(const double & u) {

  double K_val = 0.0;

  if (abs(u) < 1.0) {

    K_val = 3.28125 * (-5.0 + 42.0 * pow(u, 2) - 45.0 * pow(u, 4));
  }

  return K_val;
}

// kernel density estimation

// [[Rcpp::export]]
arma::vec KDE_K2Ep_rcpp(const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K2Ep_w_rcpp(const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K2Bw_rcpp(const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K2Bw_w_rcpp(const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K4Bw_rcpp(const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KDE_K4Bw_w_rcpp(const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
    }
  }

  Dhat /= n_n;
  return Dhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Ep_rcpp(const arma::vec & Y,
                        const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
      Nhat(k) += Kik_h * Y(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Ep_w_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Ep_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
      Nhat(k) += Kik_h * Y(i) * w(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Bw_rcpp(const arma::vec & Y,
                        const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
      Nhat(k) += Kik_h * Y(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K2Bw_w_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K2_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
      Nhat(k) += Kik_h * Y(i) * w(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K4Bw_rcpp(const arma::vec & Y,
                        const arma::mat & X,
                        const arma::mat & x,
                        const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h;
      Nhat(k) += Kik_h * Y(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
arma::vec KNW_K4Bw_w_rcpp(const arma::vec & Y,
                          const arma::mat & X,
                          const arma::mat & x,
                          const arma::vec & h,
                          const arma::vec & w) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::vec Dhat(n_k);
  arma::vec Nhat(n_k);
  arma::vec mhat(n_k);

  arma::vec xrow_k(n_p);
  arma::vec Xrow_i(n_p);
  arma::vec Dik_h(n_p);
  double Kik_h = 1.0;

  for (size_t k = 0; k < n_k; ++k) {

    xrow_k = x.row(k).t();

    for (size_t i = 0; i < n_n; ++i) {

      Xrow_i = X.row(i).t();
      Dik_h = (Xrow_i - xrow_k) / h;
      Kik_h = 1.0;

      for (size_t l = 0; l < n_p; ++l) {

        Kik_h *= K4_Bw_rcpp(Dik_h(l)) / h(l);
      }

      Dhat(k) += Kik_h * w(i);
      Nhat(k) += Kik_h * Y(i) * w(i);
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    if (Dhat(k) != 0) {

      mhat(k) = Nhat(k) / Dhat(k);
    }
  }

  return mhat;
}

// [[Rcpp::export]]
double CVKNW_K2Ep_rcpp(const arma::vec & Y,
                       const arma::mat & X,
                       const arma::vec & h) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  arma::vec Dhat(n_n);
  arma::vec Nhat(n_n);
  arma::vec mhat(n_n);
  double cv_value = 0.0;

  arma::vec Xrow_i(n_p);
  arma::vec Xrow_j(n_p);
  arma::vec Dij_h(n_p);
  double Kij_h = 1.0;

  for (size_t j = 0; j < n_n; ++j) {

    Xrow_j = X.row(j).t();

    for (size_t i = 0; i < n_n; ++i) {

      if (i != j) {

        Xrow_i = X.row(i).t();
        Dij_h = (Xrow_i - Xrow_j) / h;
        Kij_h = 1.0;

        for (size_t l = 0; l < n_p; ++l) {

          Kij_h *= K2_Ep_rcpp(Dij_h(l)) / h(l);
        }

        Dhat(j) += Kij_h;
        Nhat(j) += Kij_h * Y(i);
      }
    }

    if (Dhat(j) != 0) {

      mhat(j) = Nhat(j) / Dhat(j);
    }

    cv_value += pow(Y(j) - mhat(j), 2);
  }

  cv_value /= n_n;

  return cv_value;
}









// For comparison

// [[Rcpp::export]]
double K2_Ep_rcpp_v1(const double & u) {

  double K_val = 0.0;

  K_val = 0.75 * (1.0 - pow(u, 2)) * (abs(u) < 1.0);

  return K_val;
}

