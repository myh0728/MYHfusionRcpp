Panel.SemiSI.Frailty <- function(X, Y, kernel = "K2.Biweight", bandwidth = NULL,
                                 bandwidth.initial = 1, beta.initial = NULL,
                                 link.eval = NULL)
{
  number_n <- dim(Y)[1]
  number_t <- dim(Y)[2]
  number_p <- dim(X)[3]

  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  # 如果使用者輸入的名稱在對應表中，就轉換；否則維持原樣
  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  if (is.null(beta.initial)) {

    beta.free.initial <- rep(0, number_p - 1)

  } else {

    beta.free.initial <- beta.initial[-1] / beta.initial[1]
  }

  results <- list()

  # === Step 1 ===

  if (is.null(bandwidth))
  {
    CVcriterion <- function(par.free)
    {
      beta <- c(1, par.free[1:(number_p - 1)])
      h <- exp(par.free[number_p])
      SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
      value <- criterion_panel_CV_rcpp(X = SI, Y = Y, h = h, kernel = cpp_kernel)

      return(value)
    }

    estimation.CV <- nlminb(start = c(beta.free.initial, log(bandwidth.initial)),
                            objective = CVcriterion)
    beta.free.initial <- estimation.CV$par[1:(number_p - 1)]
    bandwidth <- exp(estimation.CV$par[number_p])

    results$coef.CV <- c(1, beta.free.initial)
    results$bandwidth <- bandwidth
    results$details.CV <- estimation.CV
  }

  # === Step 2 ===

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SS_rcpp(X = SI, Y = Y, h = bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SS <- c(1, estimation$par)
  results$details.SS <- estimation

  # === Step 3.1 ===

  effV.inv.CVbeta <- get_Vinv_panel_frailty_rcpp(
    X = apply(aperm(X, c(3, 1, 2)) * results$coef.CV, c(2, 3), sum), Y = Y,
    h = results$bandwidth, kernel = cpp_kernel
  )

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.CVbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SSeff.CV <- c(1, estimation$par)
  results$details.SSeff.CV <- estimation

  # === Step 3.2 ===

  effV.inv.SSbeta <- get_Vinv_panel_frailty_rcpp(
    X = apply(aperm(X, c(3, 1, 2)) * results$coef.SS, c(2, 3), sum), Y = Y,
    h = results$bandwidth, kernel = cpp_kernel
  )

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.SSbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SSeff.SS <- c(1, estimation$par)
  results$details.SSeff.SS <- estimation

  # === estimate the link function ===

  if (!is.null(link.eval))
  {
    link.eval <- as.vector(link.eval)

    SI.hat.CV <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.CV, c(2, 3), sum)
    link.hat.CV <- get_linkF_panel_rcpp(X = SI.hat.CV, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    SI.hat.SS <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.SS, c(2, 3), sum)
    link.hat.SS <- get_linkF_panel_rcpp(X = SI.hat.SS, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    results$link.CV <- as.vector(link.hat.CV)
    results$link.SS <- as.vector(link.hat.SS)
  }

  # =============

  return(results)
}

Panel.SemiSI.Exchangable <- function(X, Y, kernel = "K2.Biweight", bandwidth = NULL,
                                     bandwidth.initial = 1, beta.initial = NULL,
                                     link.eval = NULL)
{
  number_n <- dim(Y)[1]
  number_t <- dim(Y)[2]
  number_p <- dim(X)[3]

  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  # 如果使用者輸入的名稱在對應表中，就轉換；否則維持原樣
  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  if (is.null(beta.initial)) {

    beta.free.initial <- rep(0, number_p - 1)

  } else {

    beta.free.initial <- beta.initial[-1] / beta.initial[1]
  }

  results <- list()

  # === Step 1 ===

  if (is.null(bandwidth))
  {
    CVcriterion <- function(par.free)
    {
      beta <- c(1, par.free[1:(number_p - 1)])
      h <- exp(par.free[number_p])
      SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
      value <- criterion_panel_CV_rcpp(X = SI, Y = Y, h = h, kernel = cpp_kernel)

      return(value)
    }

    estimation.CV <- nlminb(start = c(beta.free.initial, log(bandwidth.initial)),
                            objective = CVcriterion)
    beta.free.initial <- estimation.CV$par[1:(number_p - 1)]
    bandwidth <- exp(estimation.CV$par[number_p])

    results$coef.CV <- c(1, beta.free.initial)
    results$bandwidth <- bandwidth
    results$details.CV <- estimation.CV
  }

  # === Step 2 ===

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SS_rcpp(X = SI, Y = Y, h = bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SS <- c(1, estimation$par)
  results$details.SS <- estimation

  # === Step 3.1 ===

  effV.inv.CVbeta <- get_Vinv_panel_exchangable_rcpp(
    X = apply(aperm(X, c(3, 1, 2)) * results$coef.CV, c(2, 3), sum), Y = Y,
    h = results$bandwidth, kernel = cpp_kernel
  )

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.CVbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SSeff.CV <- c(1, estimation$par)
  results$details.SSeff.CV <- estimation

  # === Step 3.2 ===

  effV.inv.SSbeta <- get_Vinv_panel_exchangable_rcpp(
    X = apply(aperm(X, c(3, 1, 2)) * results$coef.SS, c(2, 3), sum), Y = Y,
    h = results$bandwidth, kernel = cpp_kernel
  )

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.SSbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SSeff.SS <- c(1, estimation$par)
  results$details.SSeff.SS <- estimation

  # === estimate the link function ===

  if (!is.null(link.eval))
  {
    link.eval <- as.vector(link.eval)

    SI.hat.CV <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.CV, c(2, 3), sum)
    link.hat.CV <- get_linkF_panel_rcpp(X = SI.hat.CV, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    SI.hat.SS <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.SS, c(2, 3), sum)
    link.hat.SS <- get_linkF_panel_rcpp(X = SI.hat.SS, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    results$link.CV <- as.vector(link.hat.CV)
    results$link.SS <- as.vector(link.hat.SS)
  }

  # =============

  return(results)
}

Panel.SemiSI.AR1 <- function(X, Y, kernel = "K2.Biweight", bandwidth = NULL,
                             bandwidth.initial = 1, beta.initial = NULL,
                             link.eval = NULL)
{
  number_n <- dim(Y)[1]
  number_t <- dim(Y)[2]
  number_p <- dim(X)[3]

  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  # 如果使用者輸入的名稱在對應表中，就轉換；否則維持原樣
  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  if (is.null(beta.initial)) {

    beta.free.initial <- rep(0, number_p - 1)

  } else {

    beta.free.initial <- beta.initial[-1] / beta.initial[1]
  }

  results <- list()

  # === Step 1 ===

  if (is.null(bandwidth))
  {
    CVcriterion <- function(par.free)
    {
      beta <- c(1, par.free[1:(number_p - 1)])
      h <- exp(par.free[number_p])
      SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
      value <- criterion_panel_CV_rcpp(X = SI, Y = Y, h = h, kernel = cpp_kernel)

      return(value)
    }

    estimation.CV <- nlminb(start = c(beta.free.initial, log(bandwidth.initial)),
                            objective = CVcriterion)
    beta.free.initial <- estimation.CV$par[1:(number_p - 1)]
    bandwidth <- exp(estimation.CV$par[number_p])

    results$coef.CV <- c(1, beta.free.initial)
    results$bandwidth <- bandwidth
    results$details.CV <- estimation.CV
  }

  # === Step 2 ===

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SS_rcpp(X = SI, Y = Y, h = bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SS <- c(1, estimation$par)
  results$details.SS <- estimation

  # === Step 3.1 ===

  effV.inv.CVbeta <- get_Vinv_panel_AR1_rcpp(
    X = apply(aperm(X, c(3, 1, 2)) * results$coef.CV, c(2, 3), sum), Y = Y,
    h = results$bandwidth, kernel = cpp_kernel
  )

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.CVbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SSeff.CV <- c(1, estimation$par)
  results$details.SSeff.CV <- estimation

  # === Step 3.2 ===

  effV.inv.SSbeta <- get_Vinv_panel_AR1_rcpp(
    X = apply(aperm(X, c(3, 1, 2)) * results$coef.SS, c(2, 3), sum), Y = Y,
    h = results$bandwidth, kernel = cpp_kernel
  )

  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.SSbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SSeff.SS <- c(1, estimation$par)
  results$details.SSeff.SS <- estimation

  # === estimate the link function ===

  if (!is.null(link.eval))
  {
    link.eval <- as.vector(link.eval)

    SI.hat.CV <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.CV, c(2, 3), sum)
    link.hat.CV <- get_linkF_panel_rcpp(X = SI.hat.CV, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    SI.hat.SS <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.SS, c(2, 3), sum)
    link.hat.SS <- get_linkF_panel_rcpp(X = SI.hat.SS, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    results$link.CV <- as.vector(link.hat.CV)
    results$link.SS <- as.vector(link.hat.SS)
  }

  # =============

  return(results)
}

Panel.SemiSI <- function(X, Y,
                         correlation = c("frailty", "exchangeable", "AR1"),
                         kernel = "K2.Biweight", bandwidth = NULL,
                         bandwidth.initial = 1, beta.initial = NULL,
                         link.eval = NULL)
{
  # 確保使用者輸入正確的 correlation 選項，預設為第一個 "frailty"
  correlation <- match.arg(correlation)

  number_n <- dim(Y)[1]
  number_t <- dim(Y)[2]
  number_p <- dim(X)[3]

  kernel_map <- c("K2.Epanechnikov" = "K2_Ep",
                  "Epanechnikov"    = "K2_Ep",
                  "K2.Gaussian"     = "K2_G",
                  "Gaussian"        = "K2_G",
                  "K2.Biweight"     = "K2_Bw",
                  "K4.Biweight"     = "K4_Bw")

  # 如果使用者輸入的名稱在對應表中，就轉換；否則維持原樣
  cpp_kernel <- if (kernel %in% names(kernel_map)) kernel_map[[kernel]] else kernel

  if (is.null(beta.initial)) {
    beta.free.initial <- rep(0, number_p - 1)
  } else {
    beta.free.initial <- beta.initial[-1] / beta.initial[1]
  }

  results <- list()

  # === Step 1: CV estimation ===
  if (is.null(bandwidth))
  {
    CVcriterion <- function(par.free)
    {
      beta <- c(1, par.free[1:(number_p - 1)])
      h <- exp(par.free[number_p])
      SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
      value <- criterion_panel_CV_rcpp(X = SI, Y = Y, h = h, kernel = cpp_kernel)

      return(value)
    }

    estimation.CV <- nlminb(start = c(beta.free.initial, log(bandwidth.initial)),
                            objective = CVcriterion)
    beta.free.initial <- estimation.CV$par[1:(number_p - 1)]
    bandwidth <- exp(estimation.CV$par[number_p])

    results$coef.CV <- c(1, beta.free.initial)
    results$bandwidth <- bandwidth
    results$details.CV <- estimation.CV
  }

  # === Step 2: Sum of Squares (SS) estimation ===
  SScriterion <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SS_rcpp(X = SI, Y = Y, h = bandwidth, kernel = cpp_kernel)

    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion)

  results$coef.SS <- c(1, estimation$par)
  results$details.SS <- estimation

  # === Step 3.1: Efficient SS estimation using CV beta ===
  # 根據 correlation 參數切換對應的 get_Vinv 函數
  SI_CV <- apply(aperm(X, c(3, 1, 2)) * results$coef.CV, c(2, 3), sum)

  if (correlation == "frailty") {
    effV.inv.CVbeta <- get_Vinv_panel_frailty_rcpp(X = SI_CV, Y = Y, h = results$bandwidth, kernel = cpp_kernel)
  } else if (correlation == "exchangeable") {
    effV.inv.CVbeta <- get_Vinv_panel_exchangable_rcpp(X = SI_CV, Y = Y, h = results$bandwidth, kernel = cpp_kernel)
  } else if (correlation == "AR1") {
    effV.inv.CVbeta <- get_Vinv_panel_AR1_rcpp(X = SI_CV, Y = Y, h = results$bandwidth, kernel = cpp_kernel)
  }

  SScriterion_eff_CV <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.CVbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)
    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion_eff_CV)

  results$coef.SSeff.CV <- c(1, estimation$par)
  results$details.SSeff.CV <- estimation

  # === Step 3.2: Efficient SS estimation using SS beta ===
  # 根據 correlation 參數切換對應的 get_Vinv 函數
  SI_SS <- apply(aperm(X, c(3, 1, 2)) * results$coef.SS, c(2, 3), sum)

  if (correlation == "frailty") {
    effV.inv.SSbeta <- get_Vinv_panel_frailty_rcpp(X = SI_SS, Y = Y, h = results$bandwidth, kernel = cpp_kernel)
  } else if (correlation == "exchangeable") {
    effV.inv.SSbeta <- get_Vinv_panel_exchangable_rcpp(X = SI_SS, Y = Y, h = results$bandwidth, kernel = cpp_kernel)
  } else if (correlation == "AR1") {
    effV.inv.SSbeta <- get_Vinv_panel_AR1_rcpp(X = SI_SS, Y = Y, h = results$bandwidth, kernel = cpp_kernel)
  }

  SScriterion_eff_SS <- function(par.free)
  {
    beta <- c(1, par.free)
    SI <- apply(aperm(X, c(3, 1, 2)) * beta, c(2, 3), sum)
    value <- criterion_panel_SSeff_rcpp(X = SI, Y = Y, V_inv = effV.inv.SSbeta,
                                        h = results$bandwidth, kernel = cpp_kernel)
    return(value)
  }

  estimation <- nlminb(start = beta.free.initial,
                       objective = SScriterion_eff_SS)

  results$coef.SSeff.SS <- c(1, estimation$par)
  results$details.SSeff.SS <- estimation

  # === Step 4: estimate the link function ===
  if (!is.null(link.eval))
  {
    link.eval <- as.vector(link.eval)

    SI.hat.CV <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.CV, c(2, 3), sum)
    link.hat.CV <- get_linkF_panel_rcpp(X = SI.hat.CV, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    SI.hat.SS <- apply(aperm(X, c(3, 1, 2)) * results$coef.SSeff.SS, c(2, 3), sum)
    link.hat.SS <- get_linkF_panel_rcpp(X = SI.hat.SS, Y = Y, u_points = link.eval,
                                        h = results$bandwidth, kernel = cpp_kernel)

    results$link.CV <- as.vector(link.hat.CV)
    results$link.SS <- as.vector(link.hat.SS)
  }

  # =============
  return(results)
}
