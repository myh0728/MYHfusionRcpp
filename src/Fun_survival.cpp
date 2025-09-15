#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

extern double K2_Ep_rcpp(const double & u);
extern double K2_Bw_rcpp(const double & u);
extern double K4_Bw_rcpp(const double & u);

// [[Rcpp::export]]
arma::vec KME_rcpp(const arma::vec & time_last,
                   const arma::uvec & is_event,
                   const arma::vec & time_event) {

  const arma::uword n_n = time_last.n_elem;
  const arma::uword n_event = time_event.n_elem;
  arma::vec Nhat(n_event);
  arma::vec Dhat(n_event);

  for (size_t i = 0;  i < n_n; ++i) {

    if (is_event(i) == 1) {

      for (size_t k = 0; k < n_event; ++k) {

        if (time_event(k) <= time_last(i)) {

          Dhat(k) += 1;

          if (time_event(k) == time_last(i)){

            Nhat(k) += 1;
          }

        } else {

          break;
        }
      }
    }

    if (is_event(i) == 0) {

      for (size_t k = 0; k < n_event; ++k) {

        if (time_event(k) <= time_last(i)) {

          Dhat(k) += 1;

        } else {

          break;
        }
      }
    }
  }

  const arma::vec dLhat = Nhat / Dhat;
  return(dLhat);
}

// [[Rcpp::export]]
arma::mat SKME_K2Ep_rcpp(const arma::vec & time_last,
                         const arma::uvec & is_event,
                         const arma::vec & time_event,
                         const arma::mat & X,
                         const arma::mat & x,
                         const arma::vec & h) {

  const arma::uword n_n = time_last.n_elem;
  const arma::uword n_event = time_event.n_elem;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::mat Dhat(n_k, n_event);
  arma::mat Nhat(n_k, n_event);
  arma::mat mhat(n_k, n_event);

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

      if (is_event(i) == 1) {

        for (size_t l = 0; l < n_event; ++l) {

          if (time_event(l) <= time_last(i)) {

            Dhat(k, l) += Kik_h;

            if (time_event(l) == time_last(i)){

              Nhat(k, l) += Kik_h;
            }

          } else {

            break;
          }
        }
      }

      if (is_event(i) == 0) {

        for (size_t l = 0; l < n_event; ++l) {

          if (time_event(l) <= time_last(i)) {

            Dhat(k, l) += Kik_h;

          } else {

            break;
          }
        }
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    for (size_t l = 0; l < n_event; ++l) {

      if (Dhat(k, l) != 0) {

        mhat(k, l) = Nhat(k, l) / Dhat(k, l);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
arma::mat SKME_K2Bw_rcpp(const arma::vec & time_last,
                         const arma::uvec & is_event,
                         const arma::vec & time_event,
                         const arma::mat & X,
                         const arma::mat & x,
                         const arma::vec & h) {

  const arma::uword n_n = time_last.n_elem;
  const arma::uword n_event = time_event.n_elem;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::mat Dhat(n_k, n_event);
  arma::mat Nhat(n_k, n_event);
  arma::mat mhat(n_k, n_event);

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

      if (is_event(i) == 1) {

        for (size_t l = 0; l < n_event; ++l) {

          if (time_event(l) <= time_last(i)) {

            Dhat(k, l) += Kik_h;

            if (time_event(l) == time_last(i)){

              Nhat(k, l) += Kik_h;
            }

          } else {

            break;
          }
        }
      }

      if (is_event(i) == 0) {

        for (size_t l = 0; l < n_event; ++l) {

          if (time_event(l) <= time_last(i)) {

            Dhat(k, l) += Kik_h;

          } else {

            break;
          }
        }
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    for (size_t l = 0; l < n_event; ++l) {

      if (Dhat(k, l) != 0) {

        mhat(k, l) = Nhat(k, l) / Dhat(k, l);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}

// [[Rcpp::export]]
arma::mat SKME_K4Bw_rcpp(const arma::vec & time_last,
                         const arma::uvec & is_event,
                         const arma::vec & time_event,
                         const arma::mat & X,
                         const arma::mat & x,
                         const arma::vec & h) {

  const arma::uword n_n = time_last.n_elem;
  const arma::uword n_event = time_event.n_elem;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  arma::mat Dhat(n_k, n_event);
  arma::mat Nhat(n_k, n_event);
  arma::mat mhat(n_k, n_event);

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

      if (is_event(i) == 1) {

        for (size_t l = 0; l < n_event; ++l) {

          if (time_event(l) <= time_last(i)) {

            Dhat(k, l) += Kik_h;

            if (time_event(l) == time_last(i)){

              Nhat(k, l) += Kik_h;
            }

          } else {

            break;
          }
        }
      }

      if (is_event(i) == 0) {

        for (size_t l = 0; l < n_event; ++l) {

          if (time_event(l) <= time_last(i)) {

            Dhat(k, l) += Kik_h;

          } else {

            break;
          }
        }
      }
    }
  }

  for (size_t k = 0; k < n_k; ++k) {

    for (size_t l = 0; l < n_event; ++l) {

      if (Dhat(k, l) != 0) {

        mhat(k, l) = Nhat(k, l) / Dhat(k, l);
      }
    }
  }

  mhat.elem(arma::find(mhat < 0)).fill(0);
  mhat.elem(arma::find(mhat > 1)).fill(1);

  return mhat;
}




