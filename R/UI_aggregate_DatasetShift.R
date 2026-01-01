GLMcombineAD.DatasetShift <- function(data = NULL,
                                      X.name = NULL,
                                      Y.name = NULL,
                                      X = NULL,
                                      Y = NULL,
                                      shift = "NS",       # "NS", "PPS", "CS"
                                      distribution = "normal", # "normal", "Bernoulli"
                                      method = "EL",
                                      do.SE = TRUE,
                                      initial = NULL,
                                      initial.DRM = NULL,
                                      info.EX = list(phi = NULL, ext.size = NULL, ext.var = NULL),
                                      info.EY = list(phi = NULL, ext.size = NULL, ext.var = NULL),
                                      info.EXsubY = list(phi = NULL, y.pts = NULL, ext.size = NULL, ext.var = NULL),
                                      info.EYsubX = list(phi = NULL, inclusion = NULL, ext.size = NULL, ext.var = NULL),
                                      iter.max = 10,
                                      step.rate = 2,
                                      step.max = 5,
                                      tol = 1e-5,
                                      eps.inv = 1e-7,
                                      diff.size = 1e-5) {

  # --- Data Preparation ---
  if (!is.null(data)) {
    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])
  } else {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  # Initial MLE Estimate
  MLE.initial <- GLM.MLE(X = X, Y = Y, distribution = distribution, do.SE = do.SE, initial = initial)
  results <- list(MLE.initial = MLE.initial)

  # Helper to set initial parameters
  get_initial_params <- function(dist) {

    # 1. 提取基礎參數
    out_list <- list(
      alpha = MLE.initial$alpha,
      beta  = MLE.initial$beta
    )

    # 2. 額外參數處理
    if (dist == "normal") out_list$sigma <- MLE.initial$sigma
    if (dist == "Gamma")  out_list$nu    <- MLE.initial$nu

    return(out_list)
  }

  # 執行設定
  init_pars <- get_initial_params(distribution)
  alpha.initial <- init_pars$alpha
  beta.initial  <- init_pars$beta
  if (distribution == "normal") sigma.initial <- init_pars$sigma
  if (distribution == "Gamma")  nu.initial    <- init_pars$nu

  # ==============================================================================
  # DISTRIBUTION: NORMAL
  # ==============================================================================
  if (distribution == "normal") {

    par.names <- c("alpha", paste0("beta", 1:number_p), "sigma")

    # ---------------------------------------------------------
    # Shift: No Shift (NS)
    # ---------------------------------------------------------
    if (shift == "NS") {

      # Case: EY
      if (!is.null(info.EY$phi)) {

        number_m <- 1
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_AD_EY_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EY$phi,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            kappa = if(!is.null(info.EY$ext.size)) info.EY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EY$ext.var)) info.EY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EY <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EY$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_AD_EY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EY$phi, eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial)),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- par.names

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EY$ext.var)) {

              ext.var.inv <- 1 / info.EY$ext.var

            } else {

              ext.var  <- grad_psi_AD_EY_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = info.EY$phi,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- 1 / ext.var
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_AD_EY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = phi.par, eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 * ext.var.inv / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EY$phi),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names)

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if ((method == "fast")|is.null(info.EY$ext.size)) {

            par.hat <- c(par.hat, info.EY$phi)
            names(par.hat) <- c(par.names, phi.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_AD_EY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            phi = par.hat[phi.names],
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_AD_EY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            kappa = if(!is.null(info.EY$ext.size)) info.EY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EY$ext.var)) info.EY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EY$Cov.coef <- asy.Cov.par / number_n
        }
      }

      # Case: EXsubY
      if (!is.null(info.EXsubY$phi)) {

        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- number_k * number_p
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_AD_EXsubY_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            kappa = if(!is.null(info.EXsubY$ext.size)) info.EXsubY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EXsubY$ext.var)) info.EXsubY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EXsubY <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EXsubY$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_AD_EXsubY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial)),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- par.names

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EXsubY$ext.var)) {

              ext.var.inv <- inv_sympd_rcpp(info.EXsubY$ext.var)

            } else {

              ext.var  <- grad_psi_AD_EXsubY_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = info.EXsubY$phi,
                y_pts = info.EXsubY$y.pts,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- inv_sympd_rcpp(ext.var)
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- t(matrix(par[(number_p + 3):(number_p + 2 + number_m)],
                                  nrow = number_p,
                                  ncol = number_k))
              phi.diff <- as.vector(t(info.EXsubY$phi - phi.par))

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_AD_EXsubY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EXsubY$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                           as.vector(t(info.EXsubY$phi))),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names)

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if ((method == "fast")|is.null(info.EXsubY$ext.size)) {

            par.hat <- c(par.hat, t(info.EXsubY$phi))
            names(par.hat) <- c(par.names, phi.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_AD_EXsubY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = t(matrix(par.hat[phi.names], nrow = number_p, ncol = number_k)),
            y_pts = info.EXsubY$y.pts,
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_AD_EXsubY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = t(matrix(par.hat[phi.names], nrow = number_p, ncol = number_k)),
            y_pts = info.EXsubY$y.pts,
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            kappa = if(!is.null(info.EXsubY$ext.size)) info.EXsubY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EXsubY$ext.var)) info.EXsubY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EXsubY$Cov.coef <- asy.Cov.par / number_n
        }
      }

      # Case: EYsubX
      if (!is.null(info.EYsubX$phi)) {

        number_m <- length(info.EYsubX$phi)
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_AD_EYsubX_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EYsubX$phi,
            inclusion = info.EYsubX$inclusion,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            kappa = if(!is.null(info.EYsubX$ext.size)) info.EYsubX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EYsubX$ext.var)) info.EYsubX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EYsubX <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EYsubX$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_AD_EYsubX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial)),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- par.names

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EYsubX$ext.var)) {

              ext.var.inv <- inv_sympd_rcpp(info.EYsubX$ext.var)

            } else {

              ext.var  <- grad_psi_AD_EYsubX_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = info.EYsub$phi,
                inclusion = info.EYsubX$inclusion,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- inv_sympd_rcpp(ext.var)
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              phi.diff <- info.EYsubX$phi - phi.par

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_AD_EYsubX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = phi.par, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EYsubX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EYsubX$phi),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names)

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if ((method == "fast")|is.null(info.EYsubX$ext.size)) {

            par.hat <- c(par.hat, info.EYsubX$phi)
            names(par.hat) <- c(par.names, phi.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_AD_EYsubX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            phi = par.hat[phi.names],
            inclusion = info.EYsubX$inclusion,
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_AD_EYsubX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            inclusion = info.EYsubX$inclusion,
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            kappa = if(!is.null(info.EYsubX$ext.size)) info.EYsubX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EYsubX$ext.var)) info.EYsubX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EYsubX$Cov.coef <- asy.Cov.par / number_n
        }
      }
    }

    # ---------------------------------------------------------
    # Shift: Covariate Shift (CS)
    # ---------------------------------------------------------
    if (shift == "CS") {

      if (is.null(initial.DRM)) initial.DRM <- rep(0, number_p)
      CS.names <- paste0("CSbeta", 1:number_p)

      # Case: EY
      if (!is.null(info.EY$phi)) {

        number_m <- 1
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          SS <- function(b) {

            meanPsi <- mean_psi_ADCS_EY_normal_rcpp(X = X,
                                                    alpha = alpha.initial,
                                                    beta = beta.initial,
                                                    phi = info.EY$phi,
                                                    CS_beta = b)
            sum(meanPsi ^ 2)
          }
          initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_ADCS_EY_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EY$phi,
            CS_beta = initial.DRM,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, CS.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, CS.names, drop = FALSE],
            kappa = if(!is.null(info.EY$ext.size)) info.EY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EY$ext.var)) info.EY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EY <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EY$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              CS_beta <- par[(number_p + 3):(number_p * 2 + 2)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADCS_EY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EY$phi, CS_beta = CS_beta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, CS.names)

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EY$ext.var)) {

              ext.var.inv <- 1 / info.EY$ext.var

            } else {

              nll <- function(par) {

                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                CS_beta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  SolveLagrange_ADCS_EY_normal_rcpp(
                    X = X, alpha = alpha, beta = beta,
                    phi = phi.par, CS_beta = CS_beta,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EY$phi, initial.DRM),
                                   objective = nll)
              par.hat <- estimation$par
              par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
              names(par.hat) <- c(par.names, phi.names, CS.names)
              alpha.initial <- par.hat["alpha"]
              beta.initial <- par.hat[paste0("beta", 1:number_p)]
              sigma.initial <- par.hat["sigma"]
              phi.initial <- par.hat[phi.names]
              initial.DRM <- par.hat[CS.names]

              ext.var  <- grad_psi_ADCS_EY_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = phi.initial, CS_beta = initial.DRM,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- 1 / ext.var
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              CS_beta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADCS_EY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = phi.par, CS_beta = CS_beta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 * ext.var.inv / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EY$phi, initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names, CS.names)

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if (method == "fast") {

            par.hat <- c(par.hat, info.EY$phi, initial.DRM)
            names(par.hat) <- c(par.names, phi.names, CS.names)
          } else if (is.null(info.EY$ext.size)) {

            par.hat <- c(par.hat[par.names], info.EY$phi, par.hat[CS.names])
            names(par.hat) <- c(par.names, phi.names, CS.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_ADCS_EY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            phi = par.hat[phi.names],
            CS_beta = par.hat[CS.names],
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_ADCS_EY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            CS_beta = par.hat[CS.names],
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, CS.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, CS.names, drop = FALSE],
            kappa = if(!is.null(info.EY$ext.size)) info.EY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EY$ext.var)) info.EY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EY$Cov.coef <- asy.Cov.par / number_n
        }
      }

      # Case: EXsubY
      if (!is.null(info.EXsubY$phi)) {

        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- number_k * number_p
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          SS <- function(b) {

            meanPsi <- mean_psi_ADCS_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              CS_beta = b,
              y_pts = info.EXsubY$y.pts)
            sum(meanPsi ^ 2)
          }
          initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_ADCS_EXsubY_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EXsubY$phi, CS_beta = initial.DRM, y_pts = info.EXsubY$y.pts,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, CS.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, CS.names, drop = FALSE],
            kappa = if(!is.null(info.EXsubY$ext.size)) info.EXsubY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EXsubY$ext.var)) info.EXsubY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EXsubY <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EXsubY$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              CS_beta <- par[(number_p + 3):(number_p * 2 + 2)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADCS_EXsubY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, CS_beta = CS_beta, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, CS.names)

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EXsubY$ext.var)) {

              ext.var.inv <- inv_sympd_rcpp(info.EXsubY$ext.var)

            } else {

              nll <- function(par) {

                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- t(matrix(par[(number_p + 3):(number_p + 2 + number_m)],
                                    nrow = number_p,
                                    ncol = number_k))
                CS_beta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]
                phi.diff <- as.vector(t(info.EXsubY$phi - phi.par))

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  SolveLagrange_ADCS_EXsubY_normal_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = phi.par, CS_beta = CS_beta, y_pts = info.EXsubY$y.pts,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EXsubY$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                             as.vector(t(info.EXsubY$phi)), initial.DRM),
                                   objective = nll)
              par.hat <- estimation$par
              par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
              names(par.hat) <- c(par.names, phi.names, CS.names)
              alpha.initial <- par.hat["alpha"]
              beta.initial <- par.hat[paste0("beta", 1:number_p)]
              sigma.initial <- par.hat["sigma"]
              phi.initial <- par.hat[phi.names]
              initial.DRM <- par.hat[CS.names]

              ext.var  <- grad_psi_ADCS_EXsubY_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = t(matrix(phi.initial, nrow = number_p, ncol = number_k)),
                CS_beta = initial.DRM,
                y_pts = info.EXsubY$y.pts,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- inv_sympd_rcpp(ext.var)
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- t(matrix(par[(number_p + 3):(number_p + 2 + number_m)],
                                  nrow = number_p,
                                  ncol = number_k))
              CS_beta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]
              phi.diff <- as.vector(t(info.EXsubY$phi - phi.par))

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADCS_EXsubY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, CS_beta = CS_beta, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EXsubY$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                           as.vector(t(info.EXsubY$phi)), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names, CS.names)

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if (method == "fast") {

            par.hat <- c(par.hat, info.EXsubY$phi, initial.DRM)
            names(par.hat) <- c(par.names, phi.names, CS.names)
          } else if (is.null(info.EXsubY$ext.size)) {

            par.hat <- c(par.hat[par.names], info.EXsubY$phi, par.hat[CS.names])
            names(par.hat) <- c(par.names, phi.names, CS.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_ADCS_EXsubY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = t(matrix(par.hat[phi.names], nrow = number_p, ncol = number_k)),
            CS_beta = par.hat[CS.names],
            y_pts = info.EXsubY$y.pts,
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_ADCS_EXsubY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = t(matrix(par.hat[phi.names], nrow = number_p, ncol = number_k)),
            CS_beta = par.hat[CS.names],
            y_pts = info.EXsubY$y.pts,
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, CS.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, CS.names, drop = FALSE],
            kappa = if(!is.null(info.EXsubY$ext.size)) info.EXsubY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EXsubY$ext.var)) info.EXsubY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EXsubY$Cov.coef <- asy.Cov.par / number_n
        }
      }

      # Case: EYsubX
      if (!is.null(info.EYsubX$phi)) {

        number_m <- length(info.EYsubX$phi)
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          SS <- function(b) {

            meanPsi <- mean_psi_ADCS_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EYsubX$phi,
              CS_beta = b,
              inclusion = info.EYsubX$inclusion)
            sum(meanPsi ^ 2)
          }
          initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_ADCS_EYsubX_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EYsubX$phi, CS_beta = initial.DRM, inclusion = info.EYsubX$inclusion,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, CS.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, CS.names, drop = FALSE],
            kappa = if(!is.null(info.EYsubX$ext.size)) info.EYsubX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EYsubX$ext.var)) info.EYsubX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EYsubX <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EYsubX$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              CS_beta <- par[(number_p + 3):(number_p * 2 + 2)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADCS_EYsubX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EYsubX$phi, CS_beta = CS_beta, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, CS.names)

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EYsubX$ext.var)) {

              ext.var.inv <- inv_sympd_rcpp(info.EYsubX$ext.var)

            } else {

              nll <- function(par) {

                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                CS_beta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]
                phi.diff <- info.EYsubX$phi - phi.par

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  SolveLagrange_ADCS_EYsubX_normal_rcpp(
                    X = X, alpha = alpha, beta = beta,
                    phi = phi.par, CS_beta = CS_beta, inclusion = info.EYsubX$inclusion,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EYsubX$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                             as.vector(t(info.EYsubX$phi)), initial.DRM),
                                   objective = nll)
              par.hat <- estimation$par
              par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
              names(par.hat) <- c(par.names, phi.names, CS.names)
              alpha.initial <- par.hat["alpha"]
              beta.initial <- par.hat[paste0("beta", 1:number_p)]
              sigma.initial <- par.hat["sigma"]
              phi.initial <- par.hat[phi.names]
              initial.DRM <- par.hat[CS.names]

              ext.var  <- grad_psi_ADCS_EYsubX_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = phi.initial,
                CS_beta = initial.DRM,
                inclusion = info.EYsubX$inclusion,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- inv_sympd_rcpp(ext.var)
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              CS_beta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]
              phi.diff <- as.vector(t(info.EYsubX$phi - phi.par))

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADCS_EYsubX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = phi.par, CS_beta = CS_beta, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EYsubX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                           as.vector(t(info.EYsubX$phi)), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names, CS.names)

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if (method == "fast") {

            par.hat <- c(par.hat, info.EYsubX$phi, initial.DRM)
            names(par.hat) <- c(par.names, phi.names, CS.names)
          } else if (is.null(info.EYsubX$ext.size)) {

            par.hat <- c(par.hat[par.names], info.EYsubX$phi, par.hat[CS.names])
            names(par.hat) <- c(par.names, phi.names, CS.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_ADCS_EYsubX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            phi = par.hat[phi.names],
            CS_beta = par.hat[CS.names],
            inclusion = info.EYsubX$inclusion,
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_ADCS_EYsubX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            CS_beta = par.hat[CS.names],
            inclusion = info.EYsubX$inclusion,
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, CS.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, CS.names, drop = FALSE],
            kappa = if(!is.null(info.EYsubX$ext.size)) info.EYsubX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EYsubX$ext.var)) info.EYsubX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EYsubX$Cov.coef <- asy.Cov.par / number_n
        }
      }
    }

    # ---------------------------------------------------------
    # Shift: Prior Probability Shift (PPS)
    # ---------------------------------------------------------
    if (shift == "PPS") {

      if (is.null(initial.DRM)) initial.DRM <- 0
      PPS.names <- "PPSbeta"

      # Case: EX
      if (!is.null(info.EX$phi)) {

        number_m <- number_p
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          SS <- function(b) {

            meanPsi <- mean_psi_ADPPS_EX_normal_rcpp(X = X,
                                                     alpha = alpha.initial,
                                                     beta = beta.initial,
                                                     sigma = sigma.initial,
                                                     phi = info.EX$phi,
                                                     PPS_beta = b)
            sum(meanPsi ^ 2)
          }
          initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_ADPPS_EX_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EX$phi,
            PPS_beta = initial.DRM,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EX$ext.size)) info.EX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EX$ext.var)) info.EX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EX <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EX$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              PPS_beta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EX$phi, PPS_beta = PPS_beta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, PPS.names)

            results$EX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EX$ext.var)) {

              ext.var.inv <- inv_sympd_rcpp(info.EX$ext.var)

            } else {

              nll <- function(par) {

                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                PPS_beta <- par[number_p + 3 + number_m]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  SolveLagrange_ADPPS_EX_normal_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = phi.par, PPS_beta = PPS_beta,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EX$ext.size * sum((info.EX$phi - phi.par) ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EX$phi, initial.DRM),
                                   objective = nll)
              par.hat <- estimation$par
              par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
              names(par.hat) <- c(par.names, phi.names, PPS.names)
              alpha.initial <- par.hat["alpha"]
              beta.initial <- par.hat[paste0("beta", 1:number_p)]
              sigma.initial <- par.hat["sigma"]
              phi.initial <- par.hat[phi.names]
              initial.DRM <- par.hat[PPS.names]

              ext.var  <- grad_psi_ADPPS_EX_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = phi.initial, PPS_beta = initial.DRM,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- inv_sympd_rcpp(ext.var)
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              PPS_beta <- par[number_p + 3 + number_m]
              phi.diff <- info.EX$phi - phi.par

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, PPS_beta = PPS_beta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EX$phi, initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names, PPS.names)

            results$EX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if (method == "fast") {

            par.hat <- c(par.hat, info.EX$phi, initial.DRM)
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          } else if (is.null(info.EX$ext.size)) {

            par.hat <- c(par.hat[par.names], info.EX$phi, par.hat[PPS.names])
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_ADPPS_EX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            PPS_beta = par.hat[PPS.names],
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_ADPPS_EX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            PPS_beta = par.hat[PPS.names],
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EX$ext.size)) info.EX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EX$ext.var)) info.EX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EX$Cov.coef <- asy.Cov.par / number_n
        }
      }

      # Case: EY
      if (!is.null(info.EY$phi)) {

        number_m <- 1
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          SS <- function(b) {

            meanPsi <- mean_psi_ADPPS_EY_normal_rcpp(X = X,
                                                    alpha = alpha.initial,
                                                    beta = beta.initial,
                                                    sigma = sigma.initial,
                                                    phi = info.EY$phi,
                                                    PPS_beta = b)
            sum(meanPsi ^ 2)
          }
          initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_ADPPS_EY_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EY$phi,
            PPS_beta = initial.DRM,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EY$ext.size)) info.EY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EY$ext.var)) info.EY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EY <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EY$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              PPS_beta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EY$phi, PPS_beta = PPS_beta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, PPS.names)

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EY$ext.var)) {

              ext.var.inv <- 1 / info.EY$ext.var

            } else {

              nll <- function(par) {

                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                PPS_beta <- par[number_p + 3 + number_m]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  SolveLagrange_ADPPS_EY_normal_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = phi.par, PPS_beta = PPS_beta,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EY$phi, initial.DRM),
                                   objective = nll)
              par.hat <- estimation$par
              par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
              names(par.hat) <- c(par.names, phi.names, PPS.names)
              alpha.initial <- par.hat["alpha"]
              beta.initial <- par.hat[paste0("beta", 1:number_p)]
              sigma.initial <- par.hat["sigma"]
              phi.initial <- par.hat[phi.names]
              initial.DRM <- par.hat[PPS.names]

              ext.var  <- grad_psi_ADPPS_EY_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = phi.initial, PPS_beta = initial.DRM,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- 1 / ext.var
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              PPS_beta <- par[number_p + 3 + number_m]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, PPS_beta = PPS_beta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 * ext.var.inv / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), info.EY$phi, initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names, PPS.names)

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if (method == "fast") {

            par.hat <- c(par.hat, info.EY$phi, initial.DRM)
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          } else if (is.null(info.EY$ext.size)) {

            par.hat <- c(par.hat[par.names], info.EY$phi, par.hat[PPS.names])
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_ADPPS_EY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            PPS_beta = par.hat[PPS.names],
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_ADPPS_EY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            PPS_beta = par.hat[PPS.names],
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EY$ext.size)) info.EY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EY$ext.var)) info.EY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EY$Cov.coef <- asy.Cov.par / number_n
        }
      }

      # Case: EXsubY
      if (!is.null(info.EXsubY$phi)) {

        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- number_k * number_p
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          SS <- function(b) {

            meanPsi <- mean_psi_ADPPS_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              PPS_beta = b,
              y_pts = info.EXsubY$y.pts)
            sum(meanPsi ^ 2)
          }
          initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_ADPPS_EXsubY_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EXsubY$phi, PPS_beta = initial.DRM, y_pts = info.EXsubY$y.pts,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EXsubY$ext.size)) info.EXsubY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EXsubY$ext.var)) info.EXsubY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EXsubY <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EXsubY$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              PPS_beta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EXsubY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, PPS_beta = PPS_beta, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, PPS.names)

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EXsubY$ext.var)) {

              ext.var.inv <- inv_sympd_rcpp(info.EXsubY$ext.var)

            } else {

              nll <- function(par) {

                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- t(matrix(par[(number_p + 3):(number_p + 2 + number_m)],
                                    nrow = number_p,
                                    ncol = number_k))
                PPS_beta <- par[number_p + 3 + number_m]
                phi.diff <- as.vector(t(info.EXsubY$phi - phi.par))

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  SolveLagrange_ADPPS_EXsubY_normal_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = phi.par, PPS_beta = PPS_beta, y_pts = info.EXsubY$y.pts,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EXsubY$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                             as.vector(t(info.EXsubY$phi)), initial.DRM),
                                   objective = nll)
              par.hat <- estimation$par
              par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
              names(par.hat) <- c(par.names, phi.names, PPS.names)
              alpha.initial <- par.hat["alpha"]
              beta.initial <- par.hat[paste0("beta", 1:number_p)]
              sigma.initial <- par.hat["sigma"]
              phi.initial <- par.hat[phi.names]
              initial.DRM <- par.hat[PPS.names]

              ext.var  <- grad_psi_ADPPS_EXsubY_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = t(matrix(phi.initial, nrow = number_p, ncol = number_k)),
                PPS_beta = initial.DRM,
                y_pts = info.EXsubY$y.pts,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- inv_sympd_rcpp(ext.var)
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- t(matrix(par[(number_p + 3):(number_p + 2 + number_m)],
                                  nrow = number_p,
                                  ncol = number_k))
              PPS_beta <- par[number_p + 3 + number_m]
              phi.diff <- as.vector(t(info.EXsubY$phi - phi.par))

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EXsubY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, PPS_beta = PPS_beta, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EXsubY$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                           as.vector(t(info.EXsubY$phi)), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names, PPS.names)

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if (method == "fast") {

            par.hat <- c(par.hat, info.EXsubY$phi, initial.DRM)
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          } else if (is.null(info.EXsubY$ext.size)) {

            par.hat <- c(par.hat[par.names], info.EXsubY$phi, par.hat[PPS.names])
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_ADPPS_EXsubY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = t(matrix(par.hat[phi.names], nrow = number_p, ncol = number_k)),
            PPS_beta = par.hat[PPS.names],
            y_pts = info.EXsubY$y.pts,
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_ADPPS_EXsubY_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = t(matrix(par.hat[phi.names], nrow = number_p, ncol = number_k)),
            PPS_beta = par.hat[PPS.names],
            y_pts = info.EXsubY$y.pts,
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EXsubY$ext.size)) info.EXsubY$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EXsubY$ext.var)) info.EXsubY$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EXsubY$Cov.coef <- asy.Cov.par / number_n
        }
      }

      # Case: EYsubX
      if (!is.null(info.EYsubX$phi)) {

        number_m <- length(info.EYsubX$phi)
        eta.names <- paste0("eta", 1:number_m)
        phi.names <- paste0("phi", 1:number_m)

        if (method == "fast") {

          SS <- function(b) {

            meanPsi <- mean_psi_ADPPS_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              PPS_beta = b,
              inclusion = info.EYsubX$inclusion)
            sum(meanPsi ^ 2)
          }
          initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial
          )

          if (rcond_rcpp(MLE.score$hessian) > eps.inv) {
            MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score$hessian)
          } else {
            MLE.score.H.inv <- ginv.sym.eigen(MLE.score$hessian, eps.inv)
          }

          AD.relevant <- grad_psi_ADPPS_EYsubX_normal_rcpp(
            X = X,
            alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
            phi = info.EYsubX$phi, PPS_beta = initial.DRM, inclusion = info.EYsubX$inclusion,
            eta = rep(0, number_m),
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          ADvar <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE], # 取出對參數的微分部分
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EYsubX$ext.size)) info.EYsubX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EYsubX$ext.var)) info.EYsubX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv
          )

          J.V <- ADvar$J.V
          J.V.inv <- ADvar$J.V.inv

          update_term <- MLE.score.H.inv %*%
            J.V[par.names, eta.names, drop = FALSE] %*%
            J.V.inv[eta.names, eta.names, drop = FALSE] %*%
            (-AD.relevant$mean_psi)

          par.hat <- as.vector(c(alpha.initial, beta.initial, sigma.initial) + update_term)
          names(par.hat) <- par.names

          results$EYsubX <- list(
            alpha = par.hat["alpha"],
            beta  = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            parameter = par.hat
          )
        }

        if (method == "EL") {

          if (is.null(info.EYsubX$ext.size)) {

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              PPS_beta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EYsubX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, PPS_beta = PPS_beta, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, PPS.names)

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          } else {

            if (!is.null(info.EYsubX$ext.var)) {

              ext.var.inv <- inv_sympd_rcpp(info.EYsubX$ext.var)

            } else {

              nll <- function(par) {

                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                PPS_beta <- par[number_p + 3 + number_m]
                phi.diff <- info.EYsubX$phi - phi.par

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  SolveLagrange_ADPPS_EYsubX_normal_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = phi.par, PPS_beta = PPS_beta, inclusion = info.EYsubX$inclusion,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EYsubX$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                             as.vector(t(info.EYsubX$phi)), initial.DRM),
                                   objective = nll)
              par.hat <- estimation$par
              par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
              names(par.hat) <- c(par.names, phi.names, PPS.names)
              alpha.initial <- par.hat["alpha"]
              beta.initial <- par.hat[paste0("beta", 1:number_p)]
              sigma.initial <- par.hat["sigma"]
              phi.initial <- par.hat[phi.names]
              initial.DRM <- par.hat[PPS.names]

              ext.var  <- grad_psi_ADPPS_EYsubX_normal_rcpp(
                X = X,
                alpha = alpha.initial, beta = beta.initial, sigma = sigma.initial,
                phi = phi.initial,
                PPS_beta = initial.DRM,
                inclusion = info.EYsubX$inclusion,
                eta = rep(0, number_m),
                diff_size = diff.size
              )$ADvar

              ext.var.inv <- inv_sympd_rcpp(ext.var)
            }

            nll <- function(par) {

              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              PPS_beta <- par[number_p + 3 + number_m]
              phi.diff <- as.vector(t(info.EYsubX$phi - phi.par))

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_ADPPS_EYsubX_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, PPS_beta = PPS_beta, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EYsubX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial, log(sigma.initial),
                                           as.vector(t(info.EYsubX$phi)), initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c(par.names, phi.names, PPS.names)

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta  = par.hat[paste0("beta", 1:number_p)],
              sigma = par.hat["sigma"],
              parameter = par.hat
            )
          }
        }

        if (do.SE) {

          if (method == "fast") {

            par.hat <- c(par.hat, info.EYsubX$phi, initial.DRM)
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          } else if (is.null(info.EYsubX$ext.size)) {

            par.hat <- c(par.hat[par.names], info.EYsubX$phi, par.hat[PPS.names])
            names(par.hat) <- c(par.names, phi.names, PPS.names)
          }

          MLE.score <- diff_lL_normal_rcpp(
            X = X, Y = Y,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"]
          )

          eta.hat <- SolveLagrange_ADPPS_EYsubX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            PPS_beta = par.hat[PPS.names],
            inclusion = info.EYsubX$inclusion,
            eta_initial = rep(0, number_m),
            iter_max = iter.max, step_rate = step.rate,
            step_max = step.max, tol = tol, eps_inv = eps.inv)$eta

          AD.relevant <- grad_psi_ADPPS_EYsubX_normal_rcpp(
            X = X,
            alpha = par.hat["alpha"],
            beta = par.hat[paste0("beta", 1:number_p)],
            sigma = par.hat["sigma"],
            phi = par.hat[phi.names],
            PPS_beta = par.hat[PPS.names],
            inclusion = info.EYsubX$inclusion,
            eta = eta.hat,
            diff_size = diff.size
          )

          AD.diff <- AD.relevant$mean_psi_gradient
          colnames(AD.diff) <- c(par.names, phi.names, PPS.names)

          asy.Cov.par <- GLMcombineADavar(
            name.par = par.names,
            MLE.hessian = MLE.score$hessian,
            Psi.square = AD.relevant$mean_psi_outer,
            Psi.diff.beta = AD.diff[, par.names, drop = FALSE],
            Psi.diff.phi = AD.diff[, phi.names, drop = FALSE],
            Psi.diff.theta = AD.diff[, PPS.names, drop = FALSE],
            kappa = if(!is.null(info.EYsubX$ext.size)) info.EYsubX$ext.size / number_n else NULL,
            avar.phi = if (!is.null(info.EYsubX$ext.var)) info.EYsubX$ext.var else AD.relevant$ADvar,
            eps.inv = eps.inv)$asy.Cov.par

          results$EYsubX$Cov.coef <- asy.Cov.par / number_n
        }
      }
    }
  }

  # ==============================================================================
  # DISTRIBUTION: BERNOULLI (Logistic)
  # ==============================================================================
  if (distribution == "Bernoulli") {

    # ---------------------------------------------------------
    # Shift: No Shift (NS)
    # ---------------------------------------------------------
    if (shift == "NS") {

      # Case: EY
      if (!is.null(info.EY$phi)) {}

      # Case: EXsubY
      if (!is.null(info.EXsubY$phi)) {}

      # Case: EYsubX
      if (!is.null(info.EYsubX$phi)) {}
    }

    # ---------------------------------------------------------
    # Shift: Covariate Shift (CS) (Filled in based on C++ signatures)
    # ---------------------------------------------------------
    if (shift == "CS") {
      if (is.null(initial.DRM)) initial.DRM <- rep(0, number_p)

      # Case: EY
      if (!is.null(info.EY$phi)) {}

      # Case: EXsubY
      if (!is.null(info.EXsubY$phi)) {}

      # Case: EYsubX
      if (!is.null(info.EYsubX$phi)) {}
    }

    # ---------------------------------------------------------
    # Shift: Prior Probability Shift (PPS) (Filled in based on C++ signatures)
    # ---------------------------------------------------------
    if (shift == "PPS") {
      if (is.null(initial.DRM)) initial.DRM <- 0

      # Case: EX
      if (!is.null(info.EX$phi)) {}

      # Case: EY
      if (!is.null(info.EY$phi)) {}

      # Case: EYsubX
      if (!is.null(info.EYsubX$phi)) {}
    }
  }

  return(results)
}








GLMcombineAD.DatasetShift0 <- function(
    data = NULL, X.name = NULL, Y.name = NULL,
    X = NULL, Y = NULL, shift = "NS",            # "NS": no shift, "PPS": prior probability shift, "CS": covariate shift
    distribution = "normal",                     # "normal", "Gamma", "Bernoulli")
    method = "EL", do.SE = TRUE,                 # "EL", "fast"
    initial = NULL, initial.DRM = NULL,
    info.EX = list(phi = NULL,
                   ext.size = NULL,
                   ext.var = NULL),
    info.EY = list(phi = NULL,
                   ext.size = NULL,
                   ext.var = NULL),
    info.EXsubY = list(phi = NULL,
                       y.pts = NULL,
                       ext.size = NULL,
                       ext.var = NULL),
    info.EYsubX = list(phi = NULL,
                       inclusion = NULL,
                       ext.size = NULL,
                       ext.var = NULL),
    iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-5, eps.inv = 1e-7, diff.size = 1e-5)
{
  if (!is.null(data))
  {
    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])
  }else
  {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  MLE.initial <- GLM.MLE(
    X = X, Y = Y,
    distribution = distribution,
    do.SE = do.SE,
    initial = initial)

  results <- list(MLE.initial = MLE.initial)

  if (distribution == "normal")
  {
    if (is.null(initial))
    {
      alpha.initial <- MLE.initial$alpha
      beta.initial <- MLE.initial$beta
      sigma.initial <- MLE.initial$sigma
    }else
    {
      alpha.initial <- initial[1]
      beta.initial <- initial[2:(number_p + 1)]
      sigma.initial <- initial[number_p + 2]
    }

    if (shift == "NS")
    {
      if (!is.null(info.EY$phi))
      {
        number_m <- 1

        if (is.null(info.EY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                SolveLagrange_AD_EY_normal_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EY$phi, eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial)),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma")

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- grad_psi_AD_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              eta = rep(0, number_m),
              diff_size = diff.size)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$mean_psi_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.square <- AD.score$mean_psi_outer
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma")

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- AD_EY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EY$phi)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta)$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  AD_EY_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EY$phi, eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                             log(sigma.initial), info.EY$phi),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                  paste("beta", 1:number_p, sep = ""),
                                  "sigma",
                                  paste("phi", 1:number_m, sep = ""))
              eta.tilde <- AD_EY_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADvar_EY_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                eta = eta.tilde)
              avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)
              ext.var.inv <- 1 / avar.phi
            }else
            {
              avar.phi <- info.EY$ext.var
              ext.var.inv <- 1 / avar.phi

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             info.EY)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                AD_EY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EY$phi, eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 * ext.var.inv / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")]),
              objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- AD_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi)
            eta.tilde <- AD_EY_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADvar_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              eta = eta.tilde)
            avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, info.EY$phi)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- AD_EY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = par.hat[paste("phi", 1:number_m, sep = "")])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EXsubY$phi))
      {
        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- length(info.EXsubY$phi)

        if (is.null(info.EXsubY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                AD_EXsubY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial)),
                                 objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma")

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- AD_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv%*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma")

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- AD_EXsubY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EXsubY$phi,
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta)$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EXsubY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  AD_EXsubY_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EXsubY$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             as.vector(t(info.EXsubY$phi))),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""))
              eta.tilde <- AD_EXsubY_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                y_pts = info.EXsubY$y.pts,
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADvar_EXsubY_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                y_pts = info.EXsubY$y.pts,
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
                ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EXsubY$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             as.vector(t(info.EXsubY$phi)))
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                AD_EXsubY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EXsubY$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- AD_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              y_pts = info.EXsubY$y.pts)
            eta.tilde <- AD_EXsubY_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              y_pts = info.EXsubY$y.pts,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADvar_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              y_pts = info.EXsubY$y.pts,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, t(info.EXsubY$phi))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- AD_EXsubY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = t(matrix(par.hat[paste("phi", 1:number_m, sep = "")],
                             nrow = number_p, ncol = number_k)),
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EYsubX$phi))
      {
        number_m <- length(info.EYsubX$phi)

        if (is.null(info.EYsubX$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                AD_EYsubX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial)),
                                 objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma")

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- AD_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma")

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- AD_EYsubX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta)$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EYsubX$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                phi.diff <- as.vector(info.EYsubX$phi - phi.par)

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  AD_EYsubX_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EYsubX$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                             log(sigma.initial), info.EYsubX$phi),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""))

              eta.tilde <- AD_EYsubX_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                inclusion = info.EYsubX$inclusion,
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADvar_EYsubX_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                inclusion = info.EYsubX$inclusion,
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
                ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EYsubX$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             info.EYsubX$phi)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              phi.diff <- as.vector(info.EYsubX$phi - phi.par)

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                AD_EYsubX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EYsubX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- AD_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion)
            eta.tilde <- AD_EYsubX_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADvar_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, t(info.EYsubX$phi))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- AD_EYsubX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = par.hat[paste("phi", 1:number_m, sep = "")],
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }
    }

    if (shift == "CS")
    {
      if (is.null(initial.DRM))
      {
        initial.DRM <- rep(0, number_p)
      }

      if (!is.null(info.EY$phi))
      {
        number_m <- 1

        if (is.null(info.EY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              theta <- par[(number_p + 3):(number_p * 2 + 2)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADCS_EY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EY$phi, CS_beta = theta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial),
                                           initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("theta", 1:number_p, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADCS_EY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EY$phi,
                CS_beta = b)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADCS_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              CS_beta = initial.DRM)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("theta", 1:number_p, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADCS_EY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EY$phi,
              CS_beta = par.hat[paste("theta", 1:number_p, sep = "")])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                paste("beta", 1:number_p, sep = ""),
                "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta)$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                theta <- par[
                  (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  ADCS_EY_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = phi.par, CS_beta = theta,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             info.EY$phi,
                                             initial.DRM),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    paste("theta", 1:number_p, sep = ""))
              eta.tilde <- ADCS_EY_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                CS_beta = par.tilde[paste("theta", 1:number_p, sep = "")],
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADCSvar_EY_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                CS_beta = par.tilde[paste("theta", 1:number_p, sep = "")],
                eta = eta.tilde)
              avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)
              ext.var.inv <- 1 / avar.phi
            }else
            {
              avar.phi <- info.EY$ext.var
              ext.var.inv <- 1 / avar.phi

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             info.EY$phi,
                             initial.DRM)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    paste("theta", 1:number_p, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              theta <- par[
                (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADCS_EY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, CS_beta = theta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 * ext.var.inv / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")],
                        par.tilde[paste("theta", 1:number_p, sep = "")]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                paste("theta", 1:number_p, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADCS_EY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EY$phi,
                CS_beta = b)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADCS_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              CS_beta = initial.DRM)
            eta.tilde <- ADCS_EY_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              CS_beta = initial.DRM,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADCSvar_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              CS_beta = initial.DRM,
              eta = eta.tilde)
            avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, info.EY$phi, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                paste("theta", 1:number_p, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)

          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADCS_EY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = par.hat[paste("phi", 1:number_m, sep = "")],
              CS_beta = par.hat[paste("theta", 1:number_p, sep = "")])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EXsubY$phi))
      {
        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- length(info.EXsubY$phi)

        if (is.null(info.EXsubY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              theta <- par[(number_p + 3):(number_p * 2 + 2)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADCS_EXsubY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, CS_beta = theta,
                  y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial),
                                           initial.DRM),
                                 objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("theta", 1:number_p, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADCS_EXsubY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EXsubY$phi,
                CS_beta = b,
                y_pts = info.EXsubY$y.pts)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADCS_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              CS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[ , c("alpha",
                                         paste("beta", 1:number_p, sep = ""),
                                         "sigma")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("theta", 1:number_p, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADCS_EXsubY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EXsubY$phi,
              CS_beta = par.hat[paste("theta", 1:number_p, sep = "")],
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[ , c("alpha",
                                         paste("beta", 1:number_p, sep = ""),
                                         "sigma")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta)$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EXsubY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)
                theta <- par[(number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  ADCS_EXsubY_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EXsubY$phi, CS_beta = theta,
                    y_pts = info.EXsubY$y.pts,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EXsubY$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             as.vector(t(info.EXsubY$phi)),
                                             initial.DRM),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    paste("theta", 1:number_p, sep = ""))
              eta.tilde <- ADCS_EXsubY_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                CS_beta = par.tilde[paste("theta", 1:number_p, sep = "")],
                y_pts = info.EXsubY$y.pts,
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADCSvar_EXsubY_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                CS_beta = par.tilde[paste("theta", 1:number_p, sep = "")],
                y_pts = info.EXsubY$y.pts,
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EXsubY$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             as.vector(t(info.EXsubY$phi)),
                             initial.DRM)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    paste("theta", 1:number_p, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)
              theta <- par[(number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADCS_EXsubY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, CS_beta = theta,
                  y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EXsubY$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")],
                        par.tilde[paste("theta", 1:number_p, sep = "")]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                paste("theta", 1:number_p, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADCS_EXsubY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EXsubY$phi,
                CS_beta = b,
                y_pts = info.EXsubY$y.pts)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADCS_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              CS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts)
            eta.tilde <- ADCS_EXsubY_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              CS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADCSvar_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              CS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
            ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, t(info.EXsubY$phi), initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                paste("theta", 1:number_p, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADCS_EXsubY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = t(matrix(par.hat[paste("phi", 1:number_m, sep = "")],
                             nrow = number_p, ncol = number_k)),
              CS_beta = par.hat[paste("theta", 1:number_p, sep = "")],
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EYsubX$phi))
      {
        number_m <- length(info.EYsubX$phi)

        if (is.null(info.EYsubX$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              theta <- par[(number_p + 3):(number_p * 2 + 2)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADCS_EYsubX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, CS_beta = theta,
                  inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial),
                                           initial.DRM),
                                 objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("theta", 1:number_p, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADCS_EYsubX_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EYsubX$phi,
                CS_beta = b,
                inclusion = info.EYsubX$inclusion)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADCS_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              CS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("theta", 1:number_p, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADCS_EYsubX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EYsubX$phi,
              CS_beta = par.hat[paste("theta", 1:number_p, sep = "")],
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta)$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EYsubX$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                phi.diff <- as.vector(info.EYsubX$phi - phi.par)
                theta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  ADCS_EYsubX_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EYsubX$phi, CS_beta = theta,
                    inclusion = info.EYsubX$inclusion,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EYsubX$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             info.EYsubX$phi,
                                             initial.DRM),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    paste("theta", 1:number_p, sep = ""))

              eta.tilde <- ADCS_EYsubX_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                CS_beta = par.tilde[paste("theta", 1:number_p, sep = "")],
                inclusion = info.EYsubX$inclusion,
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADCSvar_EYsubX_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                CS_beta = par.tilde[paste("theta", 1:number_p, sep = "")],
                inclusion = info.EYsubX$inclusion,
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EYsubX$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             info.EYsubX$phi,
                             initial.DRM)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    paste("theta", 1:number_p, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              phi.diff <- as.vector(info.EYsubX$phi - phi.par)
              theta <- par[(number_p + 3 + number_m):(number_p * 2 + 2 + number_m)]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADCS_EYsubX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, CS_beta = theta,
                  inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EYsubX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")],
                        par.tilde[paste("theta", 1:number_p, sep = "")]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                paste("theta", 1:number_p, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADCS_EYsubX_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EYsubX$phi,
                CS_beta = b,
                inclusion = info.EYsubX$inclusion)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADCS_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              CS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion)
            eta.tilde <- ADCS_EYsubX_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              CS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADCSvar_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              CS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
            ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, info.EYsubX$phi, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                paste("theta", 1:number_p, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADCS_EYsubX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = par.hat[paste("phi", 1:number_m, sep = "")],
              CS_beta = par.hat[paste("theta", 1:number_p, sep = "")],
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   paste("theta", 1:number_p, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, paste("theta", 1:number_p, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }
    }

    if (shift == "PPS")
    {
      if (is.null(initial.DRM))
      {
        initial.DRM <- 0
      }

      if (!is.null(info.EX$phi))
      {
        number_m <- number_p

        if (is.null(info.EX$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              theta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EX$phi, PPS_beta = theta,
                  eta_initial = rep(0, number_p),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial),
                                           initial.DRM),
                                 objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EX_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EX$phi,
                PPS_beta = b)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EX$phi,
              PPS_beta = initial.DRM)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_p, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta,
              eps.inv = eps.inv)
            J.V  <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EX$phi,
              PPS_beta = par.hat["theta"])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_p, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta)$asy.Cov.par

            results$EX$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EX$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_p)]
                phi.diff <- as.vector(info.EX$phi - phi.par)
                theta <- par[number_p + 3 + number_p]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  ADPPS_EX_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EX$phi, PPS_beta = theta,
                    eta_initial = rep(0, number_p),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EX$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             info.EX$phi,
                                             initial.DRM),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_p, sep = ""),
                                    "theta")

              eta.tilde <- ADPPS_EX_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_p, sep = "")],
                PPS_beta = par.tilde["theta"],
                eta_initial = rep(0, number_p),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADPPSvar_EX_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_p, sep = "")],
                PPS_beta = par.tilde["theta"],
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_p)
              ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_p))
            }else
            {
              avar.phi <- info.EX$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_p))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             info.EX$phi,
                             initial.DRM)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_p, sep = ""),
                                    "theta")
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_p)]
              phi.diff <- as.vector(info.EX$phi - phi.par)
              theta <- par[number_p + 3 + number_p]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EX$phi, PPS_beta = theta,
                  eta_initial = rep(0, number_p),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_p, sep = "")],
                        par.tilde["theta"]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_p, sep = ""),
                                "theta")

            results$EX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EX_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EX$phi,
                PPS_beta = b)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EX$phi,
              PPS_beta = initial.DRM)
            eta.tilde <- ADPPS_EX_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EX$phi,
              PPS_beta = initial.DRM,
              eta_initial = rep(0, number_p),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADPPSvar_EX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EX$phi,
              PPS_beta = initial.DRM,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_p)
            ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_p, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_p, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, info.EX$phi, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                "theta")

            results$EX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = par.hat[paste("phi", 1:number_p, sep = "")],
              PPS_beta = par.hat["theta"])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_p, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_p, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EX$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EY$phi))
      {
        number_m <- 1

        if (is.null(info.EY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              theta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EY$phi, PPS_beta = theta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial),
                                           initial.DRM),
                                 objective = nll)
            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EY$phi,
                PPS_beta = b)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              PPS_beta = initial.DRM)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EY$phi,
              PPS_beta = par.hat["theta"])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta)$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                theta <- par[number_p + 3 + number_m]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  ADPPS_EY_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = phi.par, PPS_beta = theta,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             info.EY$phi,
                                             initial.DRM),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    "theta")
              eta.tilde <- ADPPS_EY_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                PPS_beta = par.tilde["theta"],
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADPPSvar_EY_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                PPS_beta = par.tilde["theta"],
                eta = eta.tilde)
              avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)
              ext.var.inv <- 1 / avar.phi
            }else
            {
              avar.phi <- info.EY$ext.var
              ext.var.inv <- 1 / avar.phi

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             info.EY$phi,
                             initial.DRM)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    "theta")
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              theta <- par[number_p + 3 + number_m]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, PPS_beta = theta,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 * ext.var.inv / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")],
                        par.tilde["theta"]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                "theta")

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EY$phi,
                PPS_beta = b)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              PPS_beta = initial.DRM)
            eta.tilde <- ADPPS_EY_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              PPS_beta = initial.DRM,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADPPSvar_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi,
              PPS_beta = initial.DRM,
              eta = eta.tilde)
            avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, info.EY$phi, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                "theta")

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = par.hat[paste("phi", 1:number_m, sep = "")],
              PPS_beta = par.hat["theta"])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EXsubY$phi))
      {
        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- length(info.EXsubY$phi)

        if (is.null(info.EXsubY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              theta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EXsubY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, PPS_beta = theta,
                  y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial),
                                           initial.DRM),
                                 objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EXsubY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EXsubY$phi,
                PPS_beta = b,
                y_pts = info.EXsubY$y.pts)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              PPS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[ , c("alpha",
                                         paste("beta", 1:number_p, sep = ""),
                                         "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EXsubY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EXsubY$phi,
              PPS_beta = par.hat["theta"],
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[ , c("alpha",
                                         paste("beta", 1:number_p, sep = ""),
                                         "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta)$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EXsubY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)
                theta <- par[number_p + 3 + number_m]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  ADPPS_EXsubY_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EXsubY$phi, PPS_beta = theta,
                    y_pts = info.EXsubY$y.pts,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EXsubY$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             as.vector(t(info.EXsubY$phi)),
                                             initial.DRM),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    "theta")
              eta.tilde <- ADPPS_EXsubY_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                PPS_beta = par.tilde["theta"],
                y_pts = info.EXsubY$y.pts,
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADPPSvar_EXsubY_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                PPS_beta = par.tilde["theta"],
                y_pts = info.EXsubY$y.pts,
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EXsubY$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             as.vector(t(info.EXsubY$phi)),
                             initial.DRM)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    "theta")
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)
              theta <- par[number_p + 3 + number_m]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EXsubY_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, PPS_beta = theta,
                  y_pts = info.EXsubY$y.pts,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EXsubY$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")],
                        par.tilde["theta"]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                "theta")

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EXsubY_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EXsubY$phi,
                PPS_beta = b,
                y_pts = info.EXsubY$y.pts)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              PPS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts)
            eta.tilde <- ADPPS_EXsubY_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              PPS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADPPSvar_EXsubY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EXsubY$phi,
              PPS_beta = initial.DRM,
              y_pts = info.EXsubY$y.pts,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
            ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, t(info.EXsubY$phi), initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                "theta")

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EXsubY_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = t(matrix(par.hat[paste("phi", 1:number_m, sep = "")],
                             nrow = number_p, ncol = number_k)),
              PPS_beta = par.hat["theta"],
              y_pts = info.EXsubY$y.pts)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EYsubX$phi))
      {
        number_m <- length(info.EYsubX$phi)

        if (is.null(info.EYsubX$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              theta <- par[number_p + 3]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EYsubX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, PPS_beta = theta,
                  inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial,
                                           log(sigma.initial),
                                           initial.DRM),
                                 objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EYsubX_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EYsubX$phi,
                PPS_beta = b,
                inclusion = info.EYsubX$inclusion)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              PPS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                "theta")

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EYsubX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = info.EYsubX$phi,
              PPS_beta = par.hat["theta"],
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.theta = AD.diff.theta)$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EYsubX$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                sigma <- exp(par[number_p + 2])
                phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
                phi.diff <- as.vector(info.EYsubX$phi - phi.par)
                theta <- par[number_p + 3 + number_m]

                ll <- lL_normal_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta, sigma = sigma) -
                  ADPPS_EYsubX_normal_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta, sigma = sigma,
                    phi = info.EYsubX$phi, PPS_beta = theta,
                    inclusion = info.EYsubX$inclusion,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EYsubX$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             log(sigma.initial),
                                             info.EYsubX$phi,
                                             initial.DRM),
                                   objective = nll)

              par.tilde <- estimation$par
              par.tilde[number_p + 2] <- exp(par.tilde[number_p + 2])
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    "theta")

              eta.tilde <- ADPPS_EYsubX_normal_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                PPS_beta = par.tilde["theta"],
                inclusion = info.EYsubX$inclusion,
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADPPSvar_EYsubX_normal_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                sigma = par.tilde["sigma"],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                PPS_beta = par.tilde["theta"],
                inclusion = info.EYsubX$inclusion,
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EYsubX$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             sigma.initial,
                             info.EYsubX$phi,
                             initial.DRM)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    "sigma",
                                    paste("phi", 1:number_m, sep = ""),
                                    "theta")
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              sigma <- exp(par[number_p + 2])
              phi.par <- par[(number_p + 3):(number_p + 2 + number_m)]
              phi.diff <- as.vector(info.EYsubX$phi - phi.par)
              theta <- par[number_p + 3 + number_m]

              ll <- lL_normal_rcpp(X = X, Y = Y,
                                   alpha = alpha, beta = beta, sigma = sigma) -
                ADPPS_EYsubX_normal_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, PPS_beta = theta,
                  inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EYsubX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        log(par.tilde["sigma"]),
                        par.tilde[paste("phi", 1:number_m, sep = "")],
                        par.tilde["theta"]),
              objective = nll)

            par.hat <- estimation$par
            par.hat[number_p + 2] <- exp(par.hat[number_p + 2])
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                "theta")

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            SS <- function(b)
            {
              Psi <- ADPPS_EYsubX_normal_rcpp(
                X = X,
                alpha = alpha.initial,
                beta = beta.initial,
                sigma = sigma.initial,
                phi = info.EYsubX$phi,
                PPS_beta = b,
                inclusion = info.EYsubX$inclusion)$score
              ss <- sum(Psi ^ 2)
              return(ss)
            }
            initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par

            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial)
            AD.score <- ADPPS_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              PPS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion)
            eta.tilde <- ADPPS_EYsubX_normal_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              PPS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADPPSvar_EYsubX_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EYsubX$phi,
              PPS_beta = initial.DRM,
              inclusion = info.EYsubX$inclusion,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
            ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial, sigma.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = ""),
                             "sigma"),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 2,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, info.EYsubX$phi, initial.DRM)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                "sigma",
                                paste("phi", 1:number_m, sep = ""),
                                "theta")

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_normal_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"])
            AD.score <- ADPPS_EYsubX_normal_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              sigma = par.hat["sigma"],
              phi = par.hat[paste("phi", 1:number_m, sep = "")],
              PPS_beta = par.hat["theta"],
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   "sigma",
                                   paste("phi", 1:number_m, sep = ""),
                                   "theta")
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""),
                                        "sigma")]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.diff.theta <- AD.diff[, "theta"]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = ""),
                           "sigma"),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              Psi.diff.theta = AD.diff.theta,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }
    }

  }

  if (distribution == "Bernoulli")
  {
    if (is.null(initial))
    {
      alpha.initial <- MLE.initial$alpha
      beta.initial <- MLE.initial$beta
    }else
    {
      alpha.initial <- initial[1]
      beta.initial <- initial[2:(number_p + 1)]
    }

    if (shift == "NS")
    {
      if (!is.null(info.EY$phi))
      {
        number_m <- 1

        if (is.null(info.EY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]

              ll <- lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta) -
                AD_EY_logistic_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EY$phi, eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial),
                                 objective = nll)
            par.hat <- estimation$par
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial)
            AD.score <- AD_EY_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EY$phi)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = "")),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 1,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")])
            AD.score <- AD_EY_logistic_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              phi = info.EY$phi)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta)$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                phi.par <- par[(number_p + 2):(number_p + 1 + number_m)]

                ll <- lL_logistic_rcpp(X = X, Y = Y,
                                       alpha = alpha, beta = beta) -
                  AD_EY_logistic_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta,
                    phi = info.EY$phi, eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                             info.EY$phi),
                                   objective = nll)

              par.tilde <- estimation$par
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    paste("phi", 1:number_m, sep = ""))
              eta.tilde <- AD_EY_logistic_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADvar_EY_logistic_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                eta = eta.tilde)
              avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)
              ext.var.inv <- 1 / avar.phi
            }else
            {
              avar.phi <- info.EY$ext.var
              ext.var.inv <- 1 / avar.phi

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             info.EY)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    paste("phi", 1:number_m, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              phi.par <- par[(number_p + 2):(number_p + 1 + number_m)]

              ll <- lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta) -
                AD_EY_logistic_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EY$phi, eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EY$ext.size * (info.EY$phi - phi.par) ^ 2 * ext.var.inv / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        par.tilde[paste("phi", 1:number_m, sep = "")]),
              objective = nll)
            par.hat <- estimation$par
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                paste("phi", 1:number_m, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial)
            AD.score <- AD_EY_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EY$phi)
            eta.tilde <- AD_EY_logistic_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EY$phi,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADvar_EY_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EY$phi,
              eta = eta.tilde)
            avar.phi <- avar.phi.AB$var / (avar.phi.AB$diff ^ 2)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = "")),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 1,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, info.EY$phi)
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                paste("phi", 1:number_m, sep = ""))

            results$EY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")])
            AD.score <- AD_EY_logistic_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              phi = par.hat[paste("phi", 1:number_m, sep = "")])

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EXsubY$phi))
      {
        number_k <- 2
        number_m <- length(info.EXsubY$phi)

        if (is.null(info.EXsubY$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]

              ll <- lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta) -
                AD_EXsubY_logistic_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EXsubY$phi,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial,
                                           beta.initial),
                                 objective = nll)

            par.hat <- estimation$par
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial)
            AD.score <- AD_EXsubY_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EXsubY$phi)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = "")),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 1,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")])
            AD.score <- AD_EXsubY_logistic_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              phi = info.EXsubY$phi)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta)$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EXsubY$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                phi.par <- par[(number_p + 2):(number_p + 1 + number_m)]
                phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)

                ll <- lL_logistic_rcpp(X = X, Y = Y,
                                       alpha = alpha, beta = beta) -
                  AD_EXsubY_logistic_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta,
                    phi = info.EXsubY$phi,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EXsubY$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial,
                                             beta.initial,
                                             as.vector(t(info.EXsubY$phi))),
                                   objective = nll)

              par.tilde <- estimation$par
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    paste("phi", 1:number_m, sep = ""))
              eta.tilde <- AD_EXsubY_logistic_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADvar_EXsubY_logistic_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                phi = t(matrix(par.tilde[paste("phi", 1:number_m, sep = "")],
                               nrow = number_p, ncol = number_k)),
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EXsubY$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             as.vector(t(info.EXsubY$phi)))
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    paste("phi", 1:number_m, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              phi.par <- par[(number_p + 2):(number_p + 1 + number_m)]
              phi.diff <- as.vector(t(info.EXsubY$phi) - phi.par)

              ll <- lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta) -
                AD_EXsubY_logistic_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EXsubY$phi,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EXsubY$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        par.tilde[paste("phi", 1:number_m, sep = "")]),
              objective = nll)

            par.hat <- estimation$par
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                paste("phi", 1:number_m, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial)
            AD.score <- AD_EXsubY_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EXsubY$phi)
            eta.tilde <- AD_EXsubY_logistic_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EXsubY$phi,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADvar_EXsubY_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EXsubY$phi,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
            ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = "")),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 1,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, t(info.EXsubY$phi))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                paste("phi", 1:number_m, sep = ""))

            results$EXsubY <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")])
            AD.score <- AD_EXsubY_logistic_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              phi = t(matrix(par.hat[paste("phi", 1:number_m, sep = "")],
                             nrow = number_p, ncol = number_k)))

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EXsubY$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EXsubY$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }

      if (!is.null(info.EYsubX$phi))
      {
        number_m <- length(info.EYsubX$phi)

        if (is.null(info.EYsubX$ext.size))
        {
          if (method == "EL")
          {
            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]

              ll <- lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta) -
                AD_EYsubX_logistic_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial),
                                 objective = nll)

            par.hat <- estimation$par
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial)
            AD.score <- AD_EYsubX_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = "")),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 1,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")])
            AD.score <- AD_EYsubX_logistic_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta)$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }else
        {
          if (method == "EL")
          {
            if (is.null(info.EYsubX$ext.var))
            {
              nll <- function(par)
              {
                alpha <- par[1]
                beta <- par[2:(number_p + 1)]
                phi.par <- par[(number_p + 2):(number_p + 1 + number_m)]
                phi.diff <- as.vector(info.EYsubX$phi - phi.par)

                ll <- lL_logistic_rcpp(X = X, Y = Y,
                                       alpha = alpha, beta = beta) -
                  AD_EYsubX_logistic_SolveLagrange_rcpp(
                    X = X, alpha = alpha, beta = beta,
                    phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
                    eta_initial = rep(0, number_m),
                    iter_max = iter.max, step_rate = step.rate,
                    step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                  info.EYsubX$ext.size * sum(phi.diff ^ 2) / 2

                return(-ll)
              }

              estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                             info.EYsubX$phi),
                                   objective = nll)

              par.tilde <- estimation$par
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    paste("phi", 1:number_m, sep = ""))

              eta.tilde <- AD_EYsubX_logistic_SolveLagrange_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                inclusion = info.EYsubX$inclusion,
                eta_initial = rep(0, number_m),
                iter_max = iter.max, step_rate = step.rate,
                step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
              avar.phi.AB <- ADvar_EYsubX_logistic_rcpp(
                X = X,
                alpha = par.tilde["alpha"],
                beta = par.tilde[paste("beta", 1:number_p, sep = "")],
                phi = par.tilde[paste("phi", 1:number_m, sep = "")],
                inclusion = info.EYsubX$inclusion,
                eta = eta.tilde)
              H.inv <- inv_sympd_rcpp(
                t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
              ) %*% t(avar.phi.AB$diff)
              avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
            }else
            {
              avar.phi <- info.EYsubX$ext.var
              ext.var.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))

              par.tilde <- c(alpha.initial,
                             beta.initial,
                             info.EYsubX$phi)
              names(par.tilde) <- c("alpha",
                                    paste("beta", 1:number_p, sep = ""),
                                    paste("phi", 1:number_m, sep = ""))
            }

            nll <- function(par)
            {
              alpha <- par[1]
              beta <- par[2:(number_p + 1)]
              phi.par <- par[(number_p + 2):(number_p + 1 + number_m)]
              phi.diff <- as.vector(info.EYsubX$phi - phi.par)

              ll <- lL_logistic_rcpp(X = X, Y = Y,
                                     alpha = alpha, beta = beta) -
                AD_EYsubX_logistic_SolveLagrange_rcpp(
                  X = X, alpha = alpha, beta = beta,
                  phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
                  eta_initial = rep(0, number_m),
                  iter_max = iter.max, step_rate = step.rate,
                  step_max = step.max, tol = tol, eps_inv = eps.inv)$value -
                info.EYsubX$ext.size * sum(t(ext.var.inv * phi.diff) * phi.diff) / 2

              return(-ll)
            }

            estimation <- nlminb(
              start = c(par.tilde["alpha"],
                        par.tilde[paste("beta", 1:number_p, sep = "")],
                        par.tilde[paste("phi", 1:number_m, sep = "")]),
              objective = nll)

            par.hat <- estimation$par
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                paste("phi", 1:number_m, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (method == "fast")
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = alpha.initial,
              beta = beta.initial)
            AD.score <- AD_EYsubX_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion)
            eta.tilde <- AD_EYsubX_logistic_SolveLagrange_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion,
              eta_initial = rep(0, number_m),
              iter_max = iter.max, step_rate = step.rate,
              step_max = step.max, tol = tol, eps_inv = eps.inv)$eta
            avar.phi.AB <- ADvar_EYsubX_logistic_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              phi = info.EYsubX$phi,
              inclusion = info.EYsubX$inclusion,
              eta = eta.tilde)
            H.inv <- inv_sympd_rcpp(
              t(avar.phi.AB$diff) %*% avar.phi.AB$diff + eps.inv * diag(number_m)
            ) %*% t(avar.phi.AB$diff)
            avar.phi <- H.inv %*% avar.phi.AB$var %*% H.inv

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.inv <- ginv.sym.eigen(MLE.score.H)
            }
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            ADvar <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi),
              eps.inv = eps.inv)
            J.V <- ADvar$J.V
            J.V.inv <- ADvar$J.V.inv
            dimnames(J.V.inv) <- dimnames(J.V)

            par.hat <- as.vector(
              c(alpha.initial, beta.initial) + MLE.score.H.inv %*%
                matrix(J.V[c("alpha",
                             paste("beta", 1:number_p, sep = "")),
                           paste("eta", 1:number_m, sep = "")],
                       nrow = number_p + 1,
                       ncol = number_m) %*%
                J.V.inv[paste("eta", 1:number_m, sep = ""),
                        paste("eta", 1:number_m, sep = "")] %*%
                (-AD.score$score))
            par.hat <- c(par.hat, t(info.EYsubX$phi))
            names(par.hat) <- c("alpha",
                                paste("beta", 1:number_p, sep = ""),
                                paste("phi", 1:number_m, sep = ""))

            results$EYsubX <- list(
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              parameter = par.hat)
          }

          if (do.SE)
          {
            MLE.score <- diff_lL_logistic_rcpp(
              X = X, Y = Y,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")])
            AD.score <- AD_EYsubX_logistic_rcpp(
              X = X,
              alpha = par.hat["alpha"],
              beta = par.hat[paste("beta", 1:number_p, sep = "")],
              phi = par.hat[paste("phi", 1:number_m, sep = "")],
              inclusion = info.EYsubX$inclusion)

            MLE.score.H <- MLE.score$hessian
            AD.diff <- AD.score$score_gradient
            colnames(AD.diff) <- c("alpha",
                                   paste("beta", 1:number_p, sep = ""),
                                   paste("phi", 1:number_m, sep = ""))
            AD.diff.beta <- AD.diff[, c("alpha",
                                        paste("beta", 1:number_p, sep = ""))]
            AD.diff.phi <- AD.diff[, paste("phi", 1:number_m, sep = "")]
            AD.square <- AD.score$score_square
            asy.Cov.par <- GLMcombineADavar(
              name.par = c("alpha",
                           paste("beta", 1:number_p, sep = "")),
              MLE.hessian = MLE.score.H,
              Psi.square = AD.square,
              Psi.diff.beta = AD.diff.beta,
              Psi.diff.phi = AD.diff.phi,
              kappa = info.EYsubX$ext.size / number_n,
              avar.phi = as.matrix(avar.phi))$asy.Cov.par

            results$EYsubX$Cov.coef <- asy.Cov.par / number_n
          }
        }
      }
    }

    if (shift == "CS")
    {
      if (!is.null(info.EY$phi))
      {

      }

      if (!is.null(info.EXsubY$phi))
      {

      }

      if (!is.null(info.EYsubX$phi))
      {

      }
    }

    if (shift == "PPS")
    {
      if (!is.null(info.EX$phi))
      {

      }

      if (!is.null(info.EY$phi))
      {

      }

      if (!is.null(info.EXsubY$phi))
      {

      }

      if (!is.null(info.EYsubX$phi))
      {

      }
    }

  }

  return(results)
}




