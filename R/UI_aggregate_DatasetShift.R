GLMcombineAD.DatasetShift <- function(
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
    iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-5, eps.inv = 1e-7)
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
                AD_EY_normal_SolveLagrange_rcpp(
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
            AD.score <- AD_EY_normal_rcpp(
              X = X,
              alpha = alpha.initial,
              beta = beta.initial,
              sigma = sigma.initial,
              phi = info.EY$phi)

            MLE.score.H <- MLE.score$hessian
            if (rcond_rcpp(MLE.score.H) > eps.inv)
            {
              MLE.score.H.inv <- -inv_sympd_rcpp(-MLE.score.H)
            }else
            {
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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

  if (distribution == "Gamma")
  {
    if (is.null(initial))
    {
      alpha.initial <- MLE.initial$alpha
      beta.initial <- MLE.initial$beta
      nu.initial <- MLE.initial$nu
    }else
    {
      alpha.initial <- initial[1]
      beta.initial <- initial[2:(number_p + 1)]
      nu.initial <- initial[number_p + 2]
    }

    if (shift == "NS")
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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
              MLE.score.H.eigen <- eigen_rcpp(MLE.score.H)
              MLE.score.H.eigen.value <- as.vector(MLE.score.H.eigen$value)
              MLE.score.H.eigen.vector <- MLE.score.H.eigen$vector
              MLE.score.H.eigen.value[abs(MLE.score.H.eigen.value) < eps.inv] <- 0
              MLE.score.H.eigen.value.inv <- (MLE.score.H.eigen.value !=0) /
                (MLE.score.H.eigen.value + (MLE.score.H.eigen.value == 0))
              MLE.score.H.inv <- MLE.score.H.eigen.vector %*%
                diag(MLE.score.H.eigen.value.inv,
                     nrow = length(MLE.score.H.eigen.value),
                     ncol = length(MLE.score.H.eigen.value)) %*%
                t( MLE.score.H.eigen.vector)
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




# CombineAD.DatasetShift <- function(
#     data = NULL, X.name = NULL, Y.name = NULL,
#     X = NULL, Y = NULL, shift = "NS",            # "NS": no shift, "PPS": prior probability shift, "CS": covariate shift
#     distribution = "normal",                     # "normal", "Gamma", "Bernoulli"
#     ext.sample.size = NULL,                      # NULL: no uncertainty
#     ext.var = NULL,
#     method = "fast",                             # "EL", "fast"
#     do.SE = TRUE,
#     initial = NULL, initial.DRM = NULL,
#     info.EX = list(phi = NULL),
#     info.EY = list(phi = NULL),
#     info.EXsubY = list(phi = NULL,
#                        y.pts = NULL),
#     info.EYsubX = list(phi = NULL,
#                        inclusion = NULL),
#     iter.max = 10, step.rate = 2, step.max = 5, tol = 1e-5)
# {
#   if (!is.null(data))
#   {
#     X <- as.matrix(data[, X.name])
#     Y <- as.matrix(data[, Y.name])
#   }else
#   {
#     X <- as.matrix(X)
#     Y <- as.matrix(Y)
#   }
#
#   number_n <- dim(X)[1]
#   number_p <- dim(X)[2]
#
#   MLE.initial <- GLM.MLE(
#     X = X, Y = Y,
#     distribution = distribution,
#     do.SE = do.SE,
#     initial = initial)
#
#   if (method == "EL")
#   {
#     if (distribution == "normal")
#     {
#       if (is.null(initial))
#       {
#         alpha.initial <- MLE.initial$alpha
#         beta.initial <- MLE.initial$beta
#         sigma.initial <- MLE.initial$sigma
#       }else
#       {
#         alpha.initial <- initial[1]
#         beta.initial <- initial[2:(number_p + 1)]
#         sigma.initial <- initial[number_p + 2]
#       }
#
#       if (shift == "NS")
#       {
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 AD_EY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EY$phi, eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m, number_l + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 AD_EY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EY$phi, eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EY$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = theta.hat[paste("phi", 1:number_m, sep = "")])
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2, number_l + number_m * 2)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 AD_EXsubY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EXsubY$phi,
#                 y_pts = info.EXsubY$y.pts)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m, number_l + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 AD_EXsubY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EXsubY$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                             dim(info.EXsubY$phi)),
#                 y_pts = info.EXsubY$y.pts)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2, number_l + number_m * 2)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 AD_EYsubX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 inclusion = info.EYsubX$inclusion)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m, number_l + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 AD_EYsubX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EYsubX$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#                 inclusion = info.EYsubX$inclusion)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2, number_l + number_m * 2)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#       }
#
#       if (shift == "PPS")
#       {
#         if (is.null(initial.DRM))
#         {
#           initial.DRM <- 0
#         }
#
#         if (is.null(info.EX$phi))
#         {
#           results.EX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[number_p + 3]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EX$phi, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   "theta")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EX$phi,
#                 PPS_beta = theta.hat["theta"])
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EX$phi - phi.par)
#               theta <- theta.beta.phi[number_p + 3 + number_m]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = phi.par, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EX$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#                 PPS_beta = theta.hat["theta"])
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""),"sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""),"sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[number_p + 3]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EY$phi, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   "theta")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi,
#                 PPS_beta = theta.hat["theta"])
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               theta <- theta.beta.phi[number_p + 3 + number_m]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = phi.par, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EY$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#                 PPS_beta = theta.hat["theta"])
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[number_p + 3]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EXsubY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EXsubY$phi, PPS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   "theta")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 PPS_beta = theta.hat["theta"],
#                 phi = info.EXsubY$phi,
#                 y_pts = info.EXsubY$y.pts)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#               theta <- theta.beta.phi[number_p + 3 + number_m]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EXsubY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = array(phi.par, dim(info.EXsubY$phi)),
#                   PPS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EXsubY$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                             dim(info.EXsubY$phi)),
#                 PPS_beta = theta.hat["theta"],
#                 y_pts = info.EXsubY$y.pts)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[number_p + 3]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EYsubX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EYsubX$phi, PPS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   "theta")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 PPS_beta = theta.hat["theta"],
#                 inclusion = info.EYsubX$inclusion)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#               theta <- theta.beta.phi[number_p + 3 + number_m]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADPPS_EYsubX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = phi.par, PPS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EYsubX$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#                 PPS_beta = theta.hat["theta"],
#                 inclusion = info.EYsubX$inclusion)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#       }
#
#       if (shift == "CS")
#       {
#         if (is.null(initial.DRM))
#         {
#           initial.DRM <- rep(0, number_p)
#         }
#
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[(number_p + 3):(number_p * 2 + 2)]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADCS_EY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EY$phi, CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi,
#                 CS_beta = theta.hat[paste("theta", 1:number_p, sep = "")])
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_p + number_m, number_l + number_p + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""))] <-
#                 AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l,
#                                               (number_l + number_m + 1):(number_l + number_m + number_p))])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               theta <- theta.beta.phi[
#                 (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADCS_EY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = phi.par, CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EY$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#                 CS_beta = theta.hat[paste("theta", 1:number_p, sep = "")])
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + number_p, number_l + number_m * 2 + number_p)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[(number_p + 3):(number_p * 2 + 2)]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADCS_EXsubY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EXsubY$phi, CS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 CS_beta = theta.hat[paste("theta", 1:number_p, sep = "")],
#                 phi = info.EXsubY$phi,
#                 y_pts = info.EXsubY$y.pts)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_p + number_m, number_l + number_p + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""))] <-
#                 AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l,
#                                               (number_l + number_m + 1):(number_l + number_m + number_p))])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#               theta <- theta.beta.phi[
#                 (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADCS_EXsubY_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = array(phi.par, dim(info.EXsubY$phi)),
#                   CS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EXsubY$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                             dim(info.EXsubY$phi)),
#                 CS_beta = theta.hat[paste("theta", 1:number_p, sep = "")],
#                 y_pts = info.EXsubY$y.pts)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + number_p,
#                 number_l + number_m * 2 + number_p)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               sigma <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[(number_p + 3):(number_p * 2 + 2)]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADCS_EYsubX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = info.EYsubX$phi, CS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 CS_beta = theta.hat[paste("theta", 1:number_p, sep = "")],
#                 inclusion = info.EYsubX$inclusion)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_p + number_m, number_l + number_p + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""))] <-
#                 AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l,
#                                               (number_l + number_m + 1):(number_l + number_m + number_p))])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               sigma <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#               theta <- theta.beta.phi[
#                 (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]
#
#               ll <- lL_normal_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, sigma = sigma) -
#                 ADCS_EYsubX_normal_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, sigma = sigma,
#                   phi = phi.par, CS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(sigma.initial), info.EYsubX$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#                 CS_beta = theta.hat[paste("theta", 1:number_p, sep = "")],
#                 inclusion = info.EYsubX$inclusion)
#
#               number_l <- number_p + 2
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + number_p,
#                 number_l + number_m * 2 + number_p)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV,
#                                           diag(number_l + number_m * 2 + number_p))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#       }
#     }
#
#     if (distribution == "Gamma")
#     {
#       if (is.null(initial))
#       {
#         alpha.initial <- MLE.initial$alpha
#         beta.initial <- MLE.initial$beta
#         nu.initial <- MLE.initial$nu
#       }else
#       {
#         alpha.initial <- initial[1]
#         beta.initial <- initial[2:(number_p + 1)]
#         nu.initial <- initial[number_p + 2]
#       }
#
#       if (shift == "NS")
#       {
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 AD_EY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EY$phi, eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 AD_EY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EY$phi, eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EY$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 AD_EXsubY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                    alpha = alpha, beta = beta, nu = nu) -
#                 AD_EXsubY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EXsubY$phi, y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EXsubY$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 AD_EYsubX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 AD_EYsubX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EYsubX$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#       }
#
#       if (shift == "PPS")
#       {
#         if (is.null(initial.DRM))
#         {
#           initial.DRM <- -1e-3
#         }
#
#         if (is.null(info.EX$phi))
#         {
#           results.EX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#               theta <- -exp(theta.beta[number_p + 3])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EX$phi, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3] <- -exp(theta.hat[number_p + 3])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   "theta")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EX$phi - phi.par)
#               theta <- -exp(theta.beta.phi[number_p + 3 + number_m])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = phi.par, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EX$phi,
#                                            log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3 + number_m] <- -exp(theta.hat[number_p + 3 + number_m])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#               theta <- -exp(theta.beta[number_p + 3])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EY$phi, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3] <- -exp(theta.hat[number_p + 3])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   "theta")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               theta <- -exp(theta.beta.phi[number_p + 3 + number_m])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = phi.par, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EY$phi,
#                                            log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3 + number_m] <- -exp(theta.hat[number_p + 3 + number_m])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#               theta <- -exp(theta.beta[number_p + 3])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EXsubY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EXsubY$phi, PPS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3] <- -exp(theta.hat[number_p + 3])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   "theta")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#               theta <- -exp(theta.beta.phi[number_p + 3 + number_m])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EXsubY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = array(phi.par, dim(info.EXsubY$phi)),
#                   PPS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EXsubY$phi,
#                                            log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3 + number_m] <- -exp(theta.hat[number_p + 3 + number_m])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#               theta <- -exp(theta.beta[number_p + 3])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EYsubX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EYsubX$phi, PPS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3] <- -exp(theta.hat[number_p + 3])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   "theta")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#               theta <- -exp(theta.beta.phi[number_p + 3 + number_m])
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADPPS_EYsubX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = phi.par, PPS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EYsubX$phi,
#                                            log(-initial.DRM)),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             theta.hat[number_p + 3 + number_m] <- -exp(theta.hat[number_p + 3 + number_m])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#       }
#
#       if (shift == "CS")
#       {
#         if (is.null(initial.DRM))
#         {
#           initial.DRM <- rep(0, number_p)
#         }
#
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[(number_p + 3):(number_p * 2 + 2)]
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADCS_EY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EY$phi, CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               theta <- theta.beta.phi[
#                 (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADCS_EY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = phi.par, CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EY$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[(number_p + 3):(number_p * 2 + 2)]
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADCS_EXsubY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EXsubY$phi, CS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#               theta <- theta.beta.phi[
#                 (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADCS_EXsubY_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = array(phi.par, dim(info.EXsubY$phi)),
#                   CS_beta = theta,
#                   y_pts = info.EXsubY$y.pts,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EXsubY$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               nu <- exp(theta.beta[number_p + 2])
#               theta <- theta.beta[(number_p + 3):(number_p * 2 + 2)]
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADCS_EYsubX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = info.EYsubX$phi, CS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               nu <- exp(theta.beta.phi[number_p + 2])
#               phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#               theta <- theta.beta.phi[
#                 (number_p + 3 + number_m):(number_p + 2 + number_m + number_p)]
#
#               ll <- lL_Gamma_rcpp(X = X, Y = Y,
#                                   alpha = alpha, beta = beta, nu = nu) -
#                 ADCS_EYsubX_Gamma_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta, nu = nu,
#                   phi = phi.par, CS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            log(nu.initial), info.EYsubX$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             theta.hat[number_p + 2] <- exp(theta.hat[number_p + 2])
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "nu",
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               nu = theta.hat["nu"],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#       }
#     }
#
#     if (distribution == "Bernoulli")
#     {
#       if (is.null(initial))
#       {
#         alpha.initial <- MLE.initial$alpha
#         beta.initial <- MLE.initial$beta
#       }else
#       {
#         alpha.initial <- initial[1]
#         beta.initial <- initial[2:(number_p + 1)]
#       }
#
#       if (shift == "NS")
#       {
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EY$phi, eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EY$phi, eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EY$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EXsubY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EXsubY$phi,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EXsubY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EXsubY$phi,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EXsubY$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EYsubX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EYsubX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EYsubX$phi, inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EYsubX$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#       }
#
#       if (shift == "PPS")
#       {
#         if (is.null(initial.DRM))
#         {
#           initial.DRM <- 0
#         }
#
#         if (is.null(info.EX$phi))
#         {
#           results.EX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               theta <- theta.beta[number_p + 2]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADPPS_EX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EX$phi, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "theta")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               phi.diff <- as.vector(info.EX$phi - phi.par)
#               theta <- theta.beta.phi[number_p + 2 + number_m]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADPPS_EX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = phi.par, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EX$phi, initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               theta <- theta.beta[number_p + 2]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADPPS_EY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EY$phi, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "theta")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               theta <- theta.beta.phi[number_p + 2 + number_m]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADPPS_EY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = phi.par, PPS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EY$phi, initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EXsubY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EXsubY$phi,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 AD_EXsubY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EXsubY$phi,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EXsubY$phi),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               theta <- theta.beta[number_p + 2]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADPPS_EYsubX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EYsubX$phi, PPS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "theta")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#               theta <- theta.beta.phi[number_p + 2 + number_m]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADPPS_EYsubX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = phi.par, PPS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EYsubX$phi, initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""),
#                                   "theta")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat["theta"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#       }
#
#       if (shift == "CS")
#       {
#         if (is.null(initial.DRM))
#         {
#           initial.DRM <- rep(0, number_p)
#         }
#
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               theta <- theta.beta[(number_p + 2):(number_p * 2 + 1)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADCS_EY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EY$phi, CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               theta <- theta.beta.phi[
#                 (number_p + 2 + number_m):(number_p + 1 + number_m + number_p)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADCS_EY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = phi.par, CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * (info.EY$phi - phi.par) ^ 2 * inv.ext.var / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EY$phi, initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               theta <- theta.beta[(number_p + 2):(number_p * 2 + 1)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADCS_EXsubY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EXsubY$phi, CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               phi.diff <- as.vector(info.EXsubY$phi - phi.par)
#               theta <- theta.beta.phi[
#                 (number_p + 2 + number_m):(number_p + 1 + number_m + number_p)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADCS_EXsubY_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = array(phi.par, dim(info.EXsubY$phi)),
#                   CS_beta = theta,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EXsubY$phi, initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = array(theta.hat[paste("phi", 1:number_m, sep = "")],
#                           dim(info.EXsubY$phi)),
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             nll <- function(theta.beta)
#             {
#               alpha <- theta.beta[1]
#               beta <- theta.beta[2:(number_p + 1)]
#               theta <- theta.beta[(number_p + 2):(number_p * 2 + 1)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADCS_EYsubX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = info.EYsubX$phi, CS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             nll <- function(theta.beta.phi)
#             {
#               alpha <- theta.beta.phi[1]
#               beta <- theta.beta.phi[2:(number_p + 1)]
#               phi.par <- theta.beta.phi[(number_p + 2):(number_p + 1 + number_m)]
#               phi.diff <- as.vector(info.EYsubX$phi - phi.par)
#               theta <- theta.beta.phi[
#                 (number_p + 2 + number_m):(number_p + 1 + number_m + number_p)]
#
#               ll <- lL_logistic_rcpp(X = X, Y = Y,
#                                      alpha = alpha, beta = beta) -
#                 ADCS_EYsubX_logistic_SolveLagrange_rcpp(
#                   X = X, alpha = alpha, beta = beta,
#                   phi = phi.par, CS_beta = theta,
#                   inclusion = info.EYsubX$inclusion,
#                   eta_initial = rep(0, number_m),
#                   iter_max = iter.max, step_rate = step.rate,
#                   step_max = step.max, tol = tol)$value -
#                 ext.sample.size * sum(t(inv.ext.var * phi.diff) * phi.diff) / 2
#
#               return(-ll)
#             }
#
#             estimation <- nlminb(start = c(alpha.initial, beta.initial,
#                                            info.EYsubX$phi,
#                                            initial.DRM),
#                                  objective = nll)
#             theta.hat <- estimation$par
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   paste("phi", 1:number_m, sep = ""),
#                                   paste("theta", 1:number_p, sep = ""))
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               phi = theta.hat[paste("phi", 1:number_m, sep = "")],
#               theta = theta.hat[paste("theta", 1:number_p, sep = "")],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#
#             }
#           }
#         }
#       }
#     }
#   }
#
#   if (method == "fast")
#   {
#     if (distribution == "normal")
#     {
#       alpha.initial <- MLE.initial$alpha
#       beta.initial <- MLE.initial$beta
#       sigma.initial <- MLE.initial$sigma
#       number_l <- number_p + 2
#
#       if (shift == "NS")
#       {
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#           AD.score <- AD_EY_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EY$phi)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m, number_l + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <-
#               AD.score$score_gradient[, 1:number_l]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m, number_l + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2, number_l + number_m * 2)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2, number_l + number_m * 2)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#           AD.score <- AD_EXsubY_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EXsubY$phi,
#             y_pts = info.EXsubY$y.pts)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m, number_l + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EXsubY$phi,
#                 y_pts = info.EXsubY$y.pts)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m, number_l + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2, number_l + number_m * 2)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EXsubY$phi,
#                 y_pts = info.EXsubY$y.pts)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2, number_l + number_m * 2)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#           AD.score <- AD_EYsubX_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EYsubX$phi,
#             inclusion = info.EYsubX$inclusion)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m, number_l + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 inclusion = info.EYsubX$inclusion)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m, number_l + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- AD.score$score_gradient[, 1:number_l]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient[, 1:number_l])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2, number_l + number_m * 2)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- AD_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 inclusion = info.EYsubX$inclusion)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2, number_l + number_m * 2)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#       }
#
#       if (shift == "PPS")
#       {
#         if (is.null(info.EX$phi))
#         {
#           results.EX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EX$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#
#           SS <- function(b)
#           {
#             Psi <- ADPPS_EX_normal_rcpp(
#               X = X,
#               alpha = alpha.initial,
#               beta = beta.initial,
#               sigma = sigma.initial,
#               phi = info.EX$phi,
#               PPS_beta = b)$score
#             ss <- sum(Psi ^ 2)
#             return(ss)
#           }
#
#           if (is.null(initial.DRM))
#           {
#             initial.DRM <- nlminb(start = 0, objective = SS)$par
#           }else
#           {
#             initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par
#           }
#
#           AD.score <- ADPPS_EX_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EX$phi,
#             PPS_beta = initial.DRM)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + 1 + number_m, number_l + 1 + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#               AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#               paste("eta", 1:number_m, sep = "")] <-
#               t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EX$phi,
#                 PPS_beta = initial.DRM)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta"),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EX$phi,
#                 PPS_beta = initial.DRM)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#
#           SS <- function(b)
#           {
#             Psi <- ADPPS_EY_normal_rcpp(
#               X = X,
#               alpha = alpha.initial,
#               beta = beta.initial,
#               sigma = sigma.initial,
#               phi = info.EY$phi,
#               PPS_beta = b)$score
#             ss <- sum(Psi ^ 2)
#             return(ss)
#           }
#
#           if (is.null(initial.DRM))
#           {
#             initial.DRM <- nlminb(start = 0, objective = SS)$par
#           }else
#           {
#             initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par
#           }
#
#           AD.score <- ADPPS_EY_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EY$phi,
#             PPS_beta = initial.DRM)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + 1 + number_m, number_l + 1 + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#               AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#               paste("eta", 1:number_m, sep = "")] <-
#               t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi,
#                 PPS_beta = initial.DRM)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta"),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi,
#                 PPS_beta = initial.DRM)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#
#           SS <- function(b)
#           {
#             Psi <- ADPPS_EXsubY_normal_rcpp(
#               X = X,
#               alpha = alpha.initial,
#               beta = beta.initial,
#               sigma = sigma.initial,
#               phi = info.EXsubY$phi,
#               PPS_beta = b,
#               y_pts = info.EXsubY$y.pts)$score
#             ss <- sum(Psi ^ 2)
#             return(ss)
#           }
#
#           if (is.null(initial.DRM))
#           {
#             initial.DRM <- nlminb(start = 0, objective = SS)$par
#           }else
#           {
#             initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par
#           }
#
#           AD.score <- ADPPS_EXsubY_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EXsubY$phi,
#             PPS_beta = initial.DRM,
#             y_pts = info.EXsubY$y.pts)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m + 1, number_l + number_m + 1)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#               AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#               paste("eta", 1:number_m, sep = "")] <-
#               t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EXsubY$phi,
#                 PPS_beta = initial.DRM,
#                 y_pts = info.EXsubY$y.pts)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta"),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EXsubY$phi,
#                 PPS_beta = initial.DRM,
#                 y_pts = info.EXsubY$y.pts)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#
#           SS <- function(b)
#           {
#             Psi <- ADPPS_EYsubX_normal_rcpp(
#               X = X,
#               alpha = alpha.initial,
#               beta = beta.initial,
#               sigma = sigma.initial,
#               phi = info.EYsubX$phi,
#               PPS_beta = b,
#               inclusion = info.EYsubX$inclusion)$score
#             ss <- sum(Psi ^ 2)
#             return(ss)
#           }
#
#           if (is.null(initial.DRM))
#           {
#             initial.DRM <- nlminb(start = 0, objective = SS)$par
#           }else
#           {
#             initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par
#           }
#
#           AD.score <- ADPPS_EYsubX_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EYsubX$phi,
#             PPS_beta = initial.DRM,
#             inclusion = info.EYsubX$inclusion)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + 1 + number_m, number_l + 1 + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#               AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#               paste("eta", 1:number_m, sep = "")] <-
#               t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 PPS_beta = initial.DRM,
#                 inclusion = info.EYsubX$inclusion)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + 1 + number_m, number_l + 1 + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta")] <-
#                 AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma", "theta"),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l, number_l + number_m + 1)])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + 1 + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta",
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta",
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""), "theta"),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADPPS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 PPS_beta = initial.DRM,
#                 inclusion = info.EYsubX$inclusion)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + 1, number_l + number_m * 2 + 1)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta",
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta")] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""), "theta"),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + 1))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#       }
#
#       if (shift == "CS")
#       {
#         results.EX <- list(message = "No additional information.")
#
#         if (is.null(info.EY$phi))
#         {
#           results.EY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EY$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#
#           SS <- function(b)
#           {
#             Psi <- ADCS_EY_normal_rcpp(
#               X = X,
#               alpha = alpha.initial,
#               beta = beta.initial,
#               sigma = sigma.initial,
#               phi = info.EY$phi,
#               CS_beta = b)$score
#             ss <- sum(Psi ^ 2)
#             return(ss)
#           }
#
#           if (is.null(initial.DRM))
#           {
#             initial.DRM <- nlminb(start = rep(0, number_p),
#                                   objective = SS)$par
#           }else
#           {
#             initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par
#           }
#
#           AD.score <- ADCS_EY_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EY$phi,
#             CS_beta = initial.DRM)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + number_p + number_m, number_l + number_p + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""))] <-
#               AD.score$score_gradient[, c(1:number_l,
#                                           (number_l + number_m + 1):(number_l + number_m + number_p))]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <-
#               t(AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi,
#                 CS_beta = initial.DRM)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_p + number_m, number_l + number_p + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""))] <-
#                 AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l,
#                                               (number_l + number_m + 1):(number_l + number_m + number_p))])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- 1
#             }else
#             {
#               inv.ext.var <- 1 / ext.var
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2 + number_p,
#               number_l + number_m * 2 + number_p)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EY$phi,
#                 CS_beta = initial.DRM)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + number_p,
#                 number_l + number_m * 2 + number_p)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EXsubY$phi))
#         {
#           results.EXsubY <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EXsubY$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#
#           SS <- function(b)
#           {
#             Psi <- ADCS_EXsubY_normal_rcpp(
#               X = X,
#               alpha = alpha.initial,
#               beta = beta.initial,
#               sigma = sigma.initial,
#               phi = info.EXsubY$phi,
#               CS_beta = b,
#               y_pts = info.EXsubY$y.pts)$score
#             ss <- sum(Psi ^ 2)
#             return(ss)
#           }
#
#           if (is.null(initial.DRM))
#           {
#             initial.DRM <- nlminb(start = rep(0, number_p),
#                                   objective = SS)$par
#           }else
#           {
#             initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par
#           }
#
#           AD.score <- ADCS_EXsubY_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EXsubY$phi,
#             CS_beta = initial.DRM,
#             y_pts = info.EXsubY$y.pts)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m + number_p, number_l + number_m + number_p)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""))] <-
#               AD.score$score_gradient[, c(1:number_l,
#                                           (number_l + number_m + 1):(number_l + number_m + number_p))]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <-
#               t(AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EXsubY$phi,
#                 CS_beta = initial.DRM,
#                 y_pts = info.EXsubY$y.pts)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_p + number_m, number_l + number_p + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""))] <-
#                 AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l,
#                                               (number_l + number_m + 1):(number_l + number_m + number_p))])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2 + number_p,
#               number_l + number_m * 2 + number_p)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EXsubY <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EXsubY_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EXsubY$phi,
#                 CS_beta = initial.DRM,
#                 y_pts = info.EXsubY$y.pts)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + number_p,
#                 number_l + number_m * 2 + number_p)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EXsubY$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#
#         if (is.null(info.EYsubX$phi))
#         {
#           results.EYsubX <- list(message = "No additional information.")
#         }else
#         {
#           number_m <- length(info.EYsubX$phi)
#           MLE.score <- diff_lL_normal_rcpp(
#             X = X, Y = Y,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial)
#           MLE.score.H <- MLE.score$hessian
#           MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#
#           SS <- function(b)
#           {
#             Psi <- ADCS_EYsubX_normal_rcpp(
#               X = X,
#               alpha = alpha.initial,
#               beta = beta.initial,
#               sigma = sigma.initial,
#               phi = info.EYsubX$phi,
#               CS_beta = b,
#               inclusion = info.EYsubX$inclusion)$score
#             ss <- sum(Psi ^ 2)
#             return(ss)
#           }
#
#           if (is.null(initial.DRM))
#           {
#             initial.DRM <- nlminb(start = rep(0, number_p),
#                                   objective = SS)$par
#           }else
#           {
#             initial.DRM <- nlminb(start = initial.DRM, objective = SS)$par
#           }
#
#           AD.score <- ADCS_EYsubX_normal_rcpp(
#             X = X,
#             alpha = alpha.initial,
#             beta = beta.initial,
#             sigma = sigma.initial,
#             phi = info.EYsubX$phi,
#             CS_beta = initial.DRM,
#             inclusion = info.EYsubX$inclusion)
#
#           if (is.null(ext.sample.size))
#           {
#             asy.Cov.JV <- matrix(
#               0, number_l + number_p + number_m, number_l + number_p + number_m)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = ""))] <-
#               AD.score$score_gradient[, c(1:number_l,
#                                           (number_l + number_m + 1):(number_l + number_m + number_p))]
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("theta", 1:number_p, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <-
#               t(AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))])
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 CS_beta = initial.DRM,
#                 inclusion = info.EYsubX$inclusion)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_p + number_m, number_l + number_p + number_m)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = ""))] <-
#                 AD.score$score_gradient[, c(1:number_l,
#                                             (number_l + number_m + 1):(number_l + number_m + number_p))]
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <-
#                 t(AD.score$score_gradient[, c(1:number_l,
#                                               (number_l + number_m + 1):(number_l + number_m + number_p))])
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_p + number_m))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }else
#           {
#             if (is.null(ext.var))
#             {
#               inv.ext.var <- diag(number_m)
#             }else
#             {
#               inv.ext.var <- inv_sympd_rcpp(ext.var)
#             }
#
#             asy.Cov.JV <- matrix(
#               0, number_l + number_m * 2 + number_p,
#               number_l + number_m * 2 + number_p)
#             dimnames(asy.Cov.JV) <- list(
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""),
#                 paste("eta", 1:number_m, sep = "")))
#
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#             asy.Cov.JV[
#               paste("phi", 1:number_m, sep = ""),
#               paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#             asy.Cov.JV[
#               paste("eta", 1:number_m, sep = ""),
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#             asy.Cov.JV[
#               c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("theta", 1:number_p, sep = "")),
#               paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#             asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#             dimnames(asy.Cov.invJV) <- dimnames(asy.Cov.JV)
#             theta.hat <- c(alpha.initial, beta.initial, sigma.initial) -
#               MLE.score.invH %*%
#               matrix(t(AD.score$score_gradient[, 1:number_l]),
#                      nrow = number_l, ncol = number_m) %*%
#               asy.Cov.invJV[paste("eta", 1:number_m, sep = ""),
#                             paste("eta", 1:number_m, sep = "")] %*% AD.score$score
#             theta.hat <- as.vector(theta.hat)
#             names(theta.hat) <- c("alpha",
#                                   paste("beta", 1:number_p, sep = ""),
#                                   "sigma")
#
#             results.EYsubX <- list(
#               alpha = theta.hat["alpha"],
#               beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#               sigma = theta.hat["sigma"],
#               parameter = theta.hat)
#
#             if (do.SE)
#             {
#               MLE.score <- diff_lL_normal_rcpp(
#                 X = X, Y = Y,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"])
#               MLE.score.H <- MLE.score$hessian
#               MLE.score.invH <- -inv_sympd_rcpp(-MLE.score$hessian)
#               AD.score <- ADCS_EYsubX_normal_rcpp(
#                 X = X,
#                 alpha = theta.hat["alpha"],
#                 beta = theta.hat[paste("beta", 1:number_p, sep = "")],
#                 sigma = theta.hat["sigma"],
#                 phi = info.EYsubX$phi,
#                 CS_beta = initial.DRM,
#                 inclusion = info.EYsubX$inclusion)
#
#               asy.Cov.SigmaS <- asy.Cov.JV <- matrix(
#                 0, number_l + number_m * 2 + number_p, number_l + number_m * 2 + number_p)
#               dimnames(asy.Cov.SigmaS) <- dimnames(asy.Cov.JV) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""),
#                   paste("eta", 1:number_m, sep = "")))
#
#               asy.Cov.SigmaS[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.SigmaS[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.SigmaS[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- AD.score$score_square
#
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")] <- -MLE.score.H
#               asy.Cov.JV[
#                 paste("phi", 1:number_m, sep = ""),
#                 paste("phi", 1:number_m, sep = "")] <- ext.sample.size / number_n * inv.ext.var
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 paste("eta", 1:number_m, sep = "")] <- -AD.score$score_square
#               asy.Cov.JV[
#                 paste("eta", 1:number_m, sep = ""),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = ""))] <- AD.score$score_gradient
#               asy.Cov.JV[
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma",
#                   paste("phi", 1:number_m, sep = ""),
#                   paste("theta", 1:number_p, sep = "")),
#                 paste("eta", 1:number_m, sep = "")] <- t(AD.score$score_gradient)
#               asy.Cov.invJV <- solve_rcpp(asy.Cov.JV, diag(number_l + number_m * 2 + number_p))
#
#               asy.Cov <- asy.Cov.invJV %*% asy.Cov.SigmaS %*% asy.Cov.invJV
#               asy.Cov <- asy.Cov[1:number_l, 1:number_l]
#               dimnames(asy.Cov) <- list(
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
#                 c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))
#
#               results.EYsubX$Cov.coef <- asy.Cov / number_n
#             }
#           }
#         }
#       }
#     }
#
#     if (distribution == "Gamma")
#     {
#
#     }
#
#     if (distribution == "Bernoulli")
#     {
#
#     }
#   }
#
#   results <- list(MLE.initial = MLE.initial,
#                   EX = results.EX,
#                   EY = results.EY,
#                   EXsubY = results.EXsubY,
#                   EYsubX = results.EYsubX)
#
#   return(results)
# }
