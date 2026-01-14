SurvData_process <- function(data = NULL, X.name = NULL, Y.name = NULL, D.name = NULL,
                             X = NULL, Y = NULL, D = NULL)
{
  # 1. 資料前處理
  if (!is.null(data)) {
    # --- 從 data frame 讀取 ---
    Y <- as.vector(data[, Y.name])
    D <- as.vector(data[, D.name])

    # 判斷是否輸入了 X.name
    if (!is.null(X.name)) {
      X <- as.matrix(data[, X.name])
    } else {
      X <- NULL
    }

  } else {
    # --- 直接讀取向量/矩陣 ---
    Y <- as.vector(Y)
    D <- as.vector(D)

    # X 維持原本輸入 (可能是 NULL，也可能有值)
    if (!is.null(X)) {
      X <- as.matrix(X)
    }
  }

  # 2. 排序 (根據時間 Y 由小到大)
  ord <- order(Y)
  time_last <- Y[ord]
  is_event <- D[ord]

  # 因為 time_last 已經排序過，unique 取出的也會是排序過的唯一時間點
  time_event <- unique(time_last[is_event == 1])

  # 3. 處理共變數 (僅當 X 存在時)
  if (!is.null(X)) {
    # drop = FALSE 確保即使只有一個變數，也會維持矩陣格式 (n x 1)
    covariate <- X[ord, , drop = FALSE]
  } else {
    covariate <- NULL
  }

  results <- list(time_last = time_last,
                  is_event = is_event,
                  time_event = time_event,
                  covariate = covariate)

  return(results)
}

KaplanMeier <- function(data = NULL, Y.name = NULL, D.name = NULL,
                        X = NULL, Y = NULL, D = NULL,
                        time.points = NULL)
{
  surv.data <- SurvData_process(data = data, Y.name = Y.name, D.name = D.name,
                                Y = Y, D = D)

  results.exact <- KME_exact_rcpp(time_last = surv.data$time_last,
                                  is_event = surv.data$is_event,
                                  time_event = surv.data$time_event)

  results <- list(time = as.vector(results.exact$time),
                  survival = as.vector(results.exact$survival),
                  hazard = as.vector(results.exact$hazard),
                  cum_hazard = as.vector(results.exact$cum_hazard))

  if (!is.null(time.points))
  {
    time.points <- sort(unique(as.vector(time.points)))

    results.time.points <- KME_at_times_rcpp(time_event = surv.data$time_event,
                                             time_points = time.points,
                                             S_event = results$survival,
                                             H_event = results$cum_hazard)

    results$survival_tps <- as.vector(results.time.points$survival)
    results$cum_hazard_tps <- as.vector(results.time.points$cum_hazard)
  }

  return(results)
}


