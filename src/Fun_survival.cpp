#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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




