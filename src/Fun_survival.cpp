#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include "head_KernelFunctions.h"

// [[Rcpp::export]]
List KME_exact_rcpp(
    const arma::vec & time_last,   // [注意] 必須預先由小到大排序
    const arma::uvec & is_event,   // [注意] 對應排序後的 time_last
    const arma::vec & time_event   // [注意] 必須是唯一值且由小到大排序
) {

  const arma::uword n_n = time_last.n_elem;
  const arma::uword n_event = time_event.n_elem;

  // 準備輸出容器
  arma::vec S_out(n_event);          // Survival Probability
  arma::vec hazard_out(n_event);     // Discrete Hazard (d/n)
  arma::vec cum_hazard_out(n_event); // Cumulative Hazard (Nelson-Aalen)

  // 狀態變數
  double current_S = 1.0;
  double current_cum_H = 0.0; // Nelson-Aalen 從 0 開始累加

  // 定義指針 j: 指向 time_last (原始樣本)
  arma::uword j = 0;

  // 單層迴圈：直接遍歷每一個事件時間點
  for (arma::uword k = 0; k < n_event; ++k) {
    double t_e = time_event(k);

    // 1. 更新 Risk Set 指針 (j)
    // 快速滑動：跳過所有時間小於當前 t_e 的人
    while (j < n_n && time_last(j) < t_e) {
      j++;
    }

    // 計算 Risk Set 大小 (分母)
    double n_risk = (double)(n_n - j);

    // 初始化當前的 hazard 為 0 (以防 risk set 為 0 的邊界情況)
    double h_val = 0.0;

    if (n_risk > 0) {
      // 2. 計算 Events 數量 (分子)
      double n_events = 0.0;
      arma::uword temp_j = j;

      while (temp_j < n_n && time_last(temp_j) == t_e) {
        if (is_event(temp_j) == 1) {
          n_events++;
        }
        temp_j++;
      }

      // 3. 計算各項指標
      // Hazard (d / n)
      h_val = n_events / n_risk;

      // Update Kaplan-Meier S(t) = S(t-1) * (1 - h)
      current_S *= (1.0 - h_val);

      // Update Nelson-Aalen H(t) = H(t-1) + h
      current_cum_H += h_val;
    }

    // 4. 儲存結果
    hazard_out(k) = h_val;
    S_out(k) = current_S;
    cum_hazard_out(k) = current_cum_H;
  }

  return List::create(
    Named("time") = time_event,
    Named("survival") = S_out,
    Named("hazard") = hazard_out,
    Named("cum_hazard") = cum_hazard_out
  );
}

// [[Rcpp::export]]
List KME_at_times_rcpp(
    const arma::vec & time_event,                  // [已排序] 事件時間點
    const arma::vec & time_points,                 // [已排序] 查詢時間點
    Rcpp::Nullable<Rcpp::NumericVector> S_event = R_NilValue, // [可選] 對應的 S(t)
    Rcpp::Nullable<Rcpp::NumericVector> H_event = R_NilValue  // [可選] 對應的 H(t)
) {

  const arma::uword n_event = time_event.n_elem;
  const arma::uword n_points = time_points.n_elem;

  // 1. 檢查並準備 Survival 相關變數
  bool calc_S = false;
  arma::vec S_vec;      // 存傳入的 S_event
  arma::vec S_out;      // 存輸出的結果
  double current_S = 1.0; // 初始值

  if (S_event.isNotNull()) {
    calc_S = true;
    S_vec = as<arma::vec>(S_event);
    if (S_vec.n_elem != n_event) stop("Length of S_event must match time_event.");
    S_out.set_size(n_points);
  }

  // 2. 檢查並準備 Cumulative Hazard 相關變數
  bool calc_H = false;
  arma::vec H_vec;      // 存傳入的 H_event
  arma::vec H_out;      // 存輸出的結果
  double current_H = 0.0; // 初始值

  if (H_event.isNotNull()) {
    calc_H = true;
    H_vec = as<arma::vec>(H_event);
    if (H_vec.n_elem != n_event) stop("Length of H_event must match time_event.");
    H_out.set_size(n_points);
  }

  // 若兩者都沒傳，直接回傳空 List
  if (!calc_S && !calc_H) {
    return List::create();
  }

  // 3. 雙指針演算法 (Two-Pointer Algorithm)
  arma::uword k = 0; // 指向 event

  for (arma::uword m = 0; m < n_points; ++m) {
    double t_query = time_points(m);

    // 【核心邏輯】
    // 只要「事件時間」小於等於「查詢時間」，就更新狀態
    while (k < n_event && time_event(k) <= t_query) {
      if (calc_S) current_S = S_vec(k);
      if (calc_H) current_H = H_vec(k);
      k++;
    }

    // 儲存當下查詢結果
    if (calc_S) S_out(m) = current_S;
    if (calc_H) H_out(m) = current_H;
  }

  // 4. 建立輸出 List
  List ret;
  if (calc_S) ret["survival"] = S_out;
  if (calc_H) ret["cum_hazard"] = H_out;

  return ret;
}

// [[Rcpp::export]]
List SKME_exact_rcpp(
    const arma::vec & time_last,   // [已排序]
    const arma::uvec & is_event,   // [已排序]
    const arma::vec & time_event,  // [已排序] 唯一事件時間
    const arma::mat & X,           // 樣本共變數
    const arma::mat & x,           // 目標點
    const arma::vec & h,           // Bandwidth
    std::string kernel = "K2_Ep"   // Kernel 名稱
) {

  const arma::uword n_n = X.n_rows;
  const arma::uword n_p = X.n_cols;
  const arma::uword n_k = x.n_rows;
  const arma::uword n_event = time_event.n_elem;

  // 取得使用者選擇的 Kernel 函數指標
  KernelFunc k_func = get_kernel_func(kernel);

  // 準備輸出容器 (初始化為 0)
  arma::mat mat_S(n_k, n_event, arma::fill::zeros);
  arma::mat mat_H(n_k, n_event, arma::fill::zeros);
  arma::mat mat_CumH(n_k, n_event, arma::fill::zeros);

  // 預先計算 1/h，減少除法運算
  arma::vec inv_h = 1.0 / h;

  // =========================================================
  // OpenMP 平行化策略：針對目標點 k 進行平行
  // 使用 schedule(dynamic) 因為不同 k 的計算時間大致相同，但
  // 若 k 很多，static 也可以。這裡用 dynamic 稍微彈性一點。
  // =========================================================
#pragma omp parallel for schedule(dynamic)
  for (arma::uword k = 0; k < n_k; ++k) {

    // [變數私有化] 每個執行緒必須有自己的權重向量
    arma::vec current_weights(n_n);
    double total_weight_sum = 0.0;

    // 1. 計算 Kernel Weights
    for (arma::uword i = 0; i < n_n; ++i) {

      double Kik_h = 1.0;

      // [記憶體優化] 直接讀取元素，不建立 Xrow 向量
      for (arma::uword l = 0; l < n_p; ++l) {
        double u = (X(i, l) - x(k, l)) * inv_h(l); // 改用乘法加速
        Kik_h *= k_func(u) * inv_h(l);
      }

      current_weights(i) = Kik_h;
      total_weight_sum += Kik_h;
    }

    // 2. 執行 Weighted Kaplan-Meier
    // 這些狀態變數宣告在迴圈內，自動成為 Thread-private
    double current_risk_weight = total_weight_sum;
    double current_S = 1.0;
    double current_cum_H = 0.0;
    arma::uword j = 0; // 用於追蹤 time_last 的指標

    for (arma::uword m = 0; m < n_event; ++m) {

      double t_e = time_event(m);

      // (A) 更新 Risk Set: 移除時間 < t_e 的樣本
      // 因為 time_last 是排序過的，可以使用線性掃描 (Two-pointer style)
      while (j < n_n && time_last(j) < t_e) {
        current_risk_weight -= current_weights(j);
        j++;
      }

      // 數值穩定性修正
      if (current_risk_weight < 1e-12) current_risk_weight = 0.0;

      double h_val = 0.0;

      // 只有當分母 > 0 時才計算
      if (current_risk_weight > 0) {

        // (B) 計算 Events (at time t_e)
        double weighted_events = 0.0;
        arma::uword temp_j = j;

        // 繼續往後掃描所有等於 t_e 的樣本
        while (temp_j < n_n && time_last(temp_j) == t_e) {
          if (is_event(temp_j) == 1) {
            weighted_events += current_weights(temp_j);
          }
          temp_j++;
        }

        // (C) 計算 Hazard
        h_val = weighted_events / current_risk_weight;

        // 強制限制 Hazard 範圍
        if (h_val < 0.0) h_val = 0.0;
        if (h_val > 1.0) h_val = 1.0;

        // 更新 S(t) 和 CumH(t)
        current_S *= (1.0 - h_val);
        current_cum_H += h_val;
      }

      // (D) 儲存結果
      // 雖然 mat_H 是共享的，但不同執行緒寫入不同的 row(k)，所以安全
      mat_H(k, m) = h_val;
      mat_S(k, m) = current_S;
      mat_CumH(k, m) = current_cum_H;
    }
  }

  return List::create(
    Named("time") = time_event,
    Named("survival") = mat_S,
    Named("hazard") = mat_H,
    Named("cum_hazard") = mat_CumH
  );
}














