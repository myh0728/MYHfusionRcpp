#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//============================================================================
// Armadillo functions transformer (to R functions)

// [[Rcpp::export]]
arma::mat pinv_rcpp(const arma::mat & M_A) {

  arma::mat M_X = arma::pinv(M_A);

  return M_X;
}

// [[Rcpp::export]]
arma::mat solve_rcpp(const arma::mat & M_A,
                     const arma::mat & M_B) {

  arma::mat M_X = arma::solve(M_A, M_B);

  return M_X;
}

// [[Rcpp::export]]
arma::mat inv_sympd_rcpp(const arma::mat & M_S) {

  arma::mat M_S_inv = arma::inv_sympd(M_S);

  return M_S_inv;
}

// [[Rcpp::export]]
List eigen_rcpp(const arma::mat & M_S) {

  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, M_S);
  eigval = arma::reverse(eigval);
  eigvec = arma::reverse(eigvec, 1);

  return Rcpp::List::create(
    Rcpp::Named("value") = eigval,
    Rcpp::Named("vector") = eigvec
  );
}

// [[Rcpp::export]]
arma::mat chol_rcpp(const arma::mat & M_S) {

  arma::mat A = arma::chol(M_S);

  return A;
}

// [[Rcpp::export]]
double rcond_rcpp(const arma::mat & A) {

  double value = arma::rcond(A);

  return value;
}

//============================================================================
// My own functions

// [[Rcpp::export]]
arma::mat GroupSum_rcpp(const arma::mat & MM,
                        const arma::uvec & id) {

  const arma::uvec id_unique = arma::unique(id);
  const arma::uword n_group = id_unique.n_elem;
  const arma::uword n_obs = MM.n_rows;
  const arma::uword n_p = MM.n_cols;
  arma::mat result(n_group, n_p);

  for (size_t i = 0; i < n_obs; ++i) {

    for (size_t j = 0; j < n_group; ++j) {

      if (id(i) == id_unique(j)) {

        result.row(j) += MM.row(i);
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::mat eXsq_rcpp(const arma::mat & data_X) {

  const arma::uword n_p = data_X.n_cols;
  arma::mat eXsq(n_p, n_p);

  for (size_t i = 0; i < n_p; ++i) {

    for (size_t j = 0; j <= i; ++j) {

      eXsq(i, j) = eXsq(j, i) = mean(data_X.col(i) % data_X.col(j));
    }
  }

  return eXsq;
}

// [[Rcpp::export]]
arma::mat eXsq_w_rcpp(const arma::mat & data_X,
                      const arma::vec & weight) {

  const arma::uword n_p = data_X.n_cols;
  arma::mat eXsq(n_p, n_p);

  for (size_t i = 0; i < n_p; ++i) {

    for (size_t j = 0; j <= i; ++j) {

      eXsq(i, j) = eXsq(j, i) = mean(data_X.col(i) % data_X.col(j) % weight);
    }
  }

  return eXsq;
}

// [[Rcpp::export]]
arma::mat Xsq_lowtri_rcpp(const arma::mat & data_X) {

  const arma::uword n_n = data_X.n_rows;
  const arma::uword n_p = data_X.n_cols;
  arma::mat Xsq_lowtri(n_n, (n_p + 1) * n_p / 2);
  arma::uword k = 0;

  for (size_t i = 0; i < n_p; ++i) {

    for (size_t j = i; j < n_p; ++j) {

      Xsq_lowtri.col(k) = data_X.col(i) % data_X.col(j);
      k += 1;
    }
  }

  return Xsq_lowtri;
}

// [[Rcpp::export]]
arma::mat twoXYsym_lowtri_rcpp(const arma::mat & data_X,
                               const arma::mat & data_Y) {

  const arma::uword n_n = data_X.n_rows;
  const arma::uword n_p = data_X.n_cols;
  arma::mat twoXYsym_lowtri(n_n, (n_p + 1) * n_p / 2);
  arma::uword k = 0;

  for (size_t i = 0; i < n_p; ++i) {

    for (size_t j = i; j < n_p; ++j) {

      if (i == j) {

        twoXYsym_lowtri.col(k) = data_X.col(i) % data_Y.col(j) * 2;

      } else {

        twoXYsym_lowtri.col(k) = data_X.col(i) % data_Y.col(j) +
          data_X.col(j) % data_Y.col(i);
      }

      k += 1;
    }
  }

  return twoXYsym_lowtri;
}

// [[Rcpp::export]]
arma::mat ctingP_rcpp(const arma::mat & Y,
                      const arma::mat & y) {

  const arma::uword n_n = Y.n_rows;
  const arma::uword n_k = y.n_rows;
  arma::mat CP(n_n, n_k);

  for (size_t i = 0; i < n_n; ++i) {

    for (size_t k = 0; k < n_k; ++k) {

      CP(i, k) = arma::prod(Y.row(i) <= y.row(k));
    }
  }

  return CP;
}

// [[Rcpp::export]]
arma::mat ctingP_uni_rcpp(const arma::vec & Y,
                          const arma::vec & y) {

  const arma::uword n_n = Y.n_elem;
  const arma::uword n_k = y.n_elem;
  arma::mat CP(n_n, n_k);

  for (size_t i = 0; i < n_n; ++i) {

    for (size_t k = 0; k < n_k; ++k) {

      CP(i, k) = (Y(i) <= y(k));
    }
  }

  return CP;
}

// [[Rcpp::export]]
arma::vec countAinB_rcpp(const arma::vec & A,
                         const arma::vec & B) {

  const arma::uword l_A = A.n_elem;
  const arma::uword l_B = B.n_elem;
  arma::vec result(l_A);

  for (size_t i = 0; i < l_A; ++i) {

    for (size_t j = 0; j < l_B; ++j) {

      if (A(i) == B(j)) {

        result(i) += 1;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::vec countAinB_W_rcpp(const arma::vec & A,
                           const arma::vec & B,
                           const arma::vec & W_B) {

  const arma::uword l_A = A.n_elem;
  const arma::uword l_B = B.n_elem;
  arma::vec result(l_A);

  for (size_t i = 0; i < l_A; ++i) {

    for (size_t j = 0; j < l_B; ++j) {

      if (A(i) == B(j)) {

        result(i) += W_B(j);
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::vec rankAinB_rcpp(const arma::vec & A,
                        const arma::vec & B) {

  const arma::uword l_A = A.n_elem;
  const arma::uword l_B = B.n_elem;
  arma::vec result(l_A);

  for (size_t i = 0; i < l_A; ++i) {

    for (size_t j = 0; j < l_B; ++j) {

      if (A(i) >= B(j)) {

        result(i) += 1;

      } else {

        break;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::vec rankAinB_W_rcpp(const arma::vec & A,
                          const arma::vec & B,
                          const arma::vec & W_B) {

  const arma::uword l_A = A.n_elem;
  const arma::uword l_B = B.n_elem;
  arma::vec result(l_A);

  for (size_t i = 0; i < l_A; ++i) {

    for (size_t j = 0; j < l_B; ++j) {

      if (A(i) >= B(j)) {

        result(i) += W_B(j);

      } else {

        break;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::mat outer_times_rcpp(const arma::vec & a,
                           const arma::vec & b) {

  const arma::uword l_a = a.n_elem;
  const arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i) {

    for (size_t j = 0; j < l_b; ++j) {

      outer(i, j) = a(i) * b(j);
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::mat outer_minus_rcpp(const arma::vec & a,
                           const arma::vec & b) {

  const arma::uword l_a = a.n_elem;
  const arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i) {

    for (size_t j = 0; j < l_b; ++j) {

      outer(i, j) = a(i) - b(j);
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::mat outer_plus_rcpp(const arma::vec & a,
                          const arma::vec & b) {

  const arma::uword l_a = a.n_elem;
  const arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i) {

    for (size_t j = 0; j < l_b; ++j) {

      outer(i, j) = a(i) + b(j);
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::mat outer_leq_rcpp(const arma::vec & a,
                         const arma::vec & b) {

  const arma::uword l_a = a.n_elem;
  const arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i) {

    for (size_t j = 0; j < l_b; ++j) {

      outer(i, j) = (a(i) <= b(j));
    }
  }

  return outer;
}

// [[Rcpp::export]]
arma::mat outer_geq_rcpp(const arma::vec & a,
                         const arma::vec & b) {

  const arma::uword l_a = a.n_elem;
  const arma::uword l_b = b.n_elem;
  arma::mat outer(l_a, l_b);

  for (size_t i = 0; i < l_a; ++i) {

    for (size_t j = 0; j < l_b; ++j) {

      outer(i, j) = (a(i) >= b(j));
    }
  }

  return outer;
}



