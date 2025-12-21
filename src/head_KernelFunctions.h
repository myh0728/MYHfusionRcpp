// 1. Header Guards 開始
//這三行是為了告訴編譯器：「如果已經讀過這個檔，就不要再讀一次」
#ifndef KERNELFUNCTION_H  // 命名慣例：全大寫，用 _ 代替點
#define KERNELFUNCTION_H

// 2. 引入必要的函式庫
#include <RcppArmadillo.h>

// ==========================================
//  定義常數 (放在函數上面)
// ==========================================

constexpr double INV_SQRT_2PI = 0.3989422804014327;

// =========================================================
//  定義 Kernel 函數 (Inline Kernels)
// =========================================================

inline double K2_Ep_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  return 0.75 * (1.0 - u * u);
}

inline double K2_Ep_diff1_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  return -1.5 * u;
}

inline double K2_Ep_diff2_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  return -1.5;
}

inline double K2_Bw_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  double t = 1.0 - u * u;
  return 0.9375 * t * t;
}

inline double K2_Bw_diff1_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  return -3.75 * u * (1.0 - u * u);
}

inline double K2_Bw_diff2_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  return -3.75 * (1.0 - 3.0 * u * u);
}

inline double K4_Bw_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  double u2 = u * u;
  double t = 1.0 - u2;
  return 1.640625 * (1.0 - 3.0 * u2) * t * t;
}

inline double K4_Bw_diff1_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  double u2 = u * u;
  return 3.28125 * u * (1.0 - u2) * (9.0 * u2 - 5.0);
}

inline double K4_Bw_diff2_inline(const double & u) {

  if (std::abs(u) >= 1.0) {

    return 0.0;
  }

  double u2 = u * u;
  return 3.28125 * (-5.0 + u2 * (42.0 - 45.0 * u2));
}

inline double K2_G_inline(const double & u) {

  if (std::abs(u) > 6.0) {

    return 0.0;
  }

  double u2 = u * u;
  return INV_SQRT_2PI * std::exp(-0.5 * u2);
}

inline double K2_G_diff1_inline(const double & u) {

  if (std::abs(u) > 6.0) {

    return 0.0;
  }

  double u2 = u * u;
  double phi = INV_SQRT_2PI * std::exp(-0.5 * u2);
  return -u * phi;
}

inline double K2_G_diff2_inline(const double & u) {

  if (std::abs(u) > 6.0) {

    return 0.0;
  }

  double u2 = u * u;
  double phi = INV_SQRT_2PI * std::exp(-0.5 * u2);
  return (u2 - 1.0) * phi;
}

#endif
