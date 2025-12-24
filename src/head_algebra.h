// 1. Header Guards 開始
//這三行是為了告訴編譯器：「如果已經讀過這個檔，就不要再讀一次」
#ifndef ALGEBRA_H  // 命名慣例：全大寫，用 _ 代替點
#define ALGEBRA_H

// 2. 引入必要的函式庫
#include <RcppArmadillo.h>
#include <cmath>

// 輔助函數：數值穩定的 log(1 + exp(x))
inline double log1pexp(double x) {

  if (x > 35.0) return x;
  if (x < -10.0) return std::exp(x);

  return std::log1p(std::exp(x));
}

// 輔助函數：數值穩定的 log(1 + exp(-x))
inline double log1pexp_neg(double x) {

  if (x > 10.0) return std::exp(-x);
  if (x < -35.0) return -x;

  return std::log1p(std::exp(-x));
}

#endif
