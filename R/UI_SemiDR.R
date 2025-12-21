SemiSI <- function(data = NULL, X.name = NULL, Y.name = NULL,
                   X = NULL, Y = NULL,
                   type = "mean", kernel = "K2.Biweight",
                   bandwidth = NULL, bandwidth.initial = 1,
                   beta.initial = NULL,
                   wi.boot = NULL)
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

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (number_p < 2)
  {
    warnings("There should be at least two covariates.")
  }

  if (is.null(beta.initial))
  {
    beta.initial <- rep(0, number_p - 1)
  }else
  {
    beta.initial <- as.vector(beta.initial)
    beta.initial <- beta.initial[-1] / beta.initial[1]
  }

  if (type == "mean")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNW_K2Ep_rcpp(Y = Y, X = X %*% beta, h = h)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNW_K2Ep_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNW_K2Ep_w_rcpp(Y = Y, X = X %*% beta, h = h, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNW_K2Ep_w_rcpp(Y = Y, X = X %*% beta, h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNW_K2Bw_rcpp(Y = Y, X = X %*% beta, h = h)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNW_K2Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNW_K2Bw_w_rcpp(Y = Y, X = X %*% beta, h = h, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNW_K2Bw_w_rcpp(Y = Y, X = X %*% beta, h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNW_K4Bw_rcpp(Y = Y, X = X %*% beta, h = h)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNW_K4Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNW_K4Bw_w_rcpp(Y = Y, X = X %*% beta, h = h, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNW_K4Bw_w_rcpp(Y = Y, X = X %*% beta, h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }
    }
  }

  if (type == "distribution")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNWcdf_K2Ep_rcpp(Y = Y, X = X %*% beta, h = h)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNWcdf_K2Ep_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNWcdf_K2Ep_w_rcpp(Y = Y, X = X %*% beta, h = h, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNWcdf_K2Ep_w_rcpp(Y = Y, X = X %*% beta, h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNWcdf_K2Bw_rcpp(Y = Y, X = X %*% beta, h = h)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNWcdf_K2Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNWcdf_K2Bw_w_rcpp(Y = Y, X = X %*% beta, h = h, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNWcdf_K2Bw_w_rcpp(Y = Y, X = X %*% beta, h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNWcdf_K4Bw_rcpp(Y = Y, X = X %*% beta, h = h)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNWcdf_K4Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- c(1, theta[1:(number_p - 1)])
            h <- exp(theta[number_p])
            cv_value <- CVKNWcdf_K4Bw_w_rcpp(Y = Y, X = X %*% beta, h = h, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(beta.initial, log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- LOOCV.esti$par[1:(number_p - 1)]
          bandwidth <- exp(LOOCV.esti$par[number_p])

          results <- list(coef.CV = c(1, beta.CV),
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- beta.CV
        }

        SemiLS <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNWcdf_K4Bw_w_rcpp(Y = Y, X = X %*% beta, h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = beta.initial,
                              objective = SemiLS)

        results$coef <- c(1, SemiLS.esti$par)
        results$LS.value <- SemiLS.esti$objective
      }
    }
  }

  return(results)
}

SemiSI.adtLASSO.ADMM <- function(data = NULL, X.name = NULL, Y.name = NULL,
                                 X = NULL, Y = NULL,
                                 type = "mean", kernel = "K2.Biweight",
                                 beta.initial, bandwidth,
                                 weight.beta, lambda.L1.penalty,
                                 rho.augLagrange.penalty = 1,
                                 iter.max = 50, tol = 1e-7)
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

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (number_p < 2)
  {
    warnings("There should be at least two covariates.")
  }

  beta.free.initial <- beta.initial[-1] / beta.initial[1]
  weight.free.beta <- lambda.L1.penalty / (weight.beta[-1] / weight.beta[1])

  if (type == "mean")
  {
    if (kernel == "K2.Epanechnikov")
    {
      beta_free_iter <- beta.free.initial
      gamma_free_iter <- beta.free.initial
      u_free_iter <- rep(0, length(beta.free.initial))

      for (iter in 1:iter.max)
      {
        SemiLS.penalized <- function(theta)
        {
          beta <- c(1, theta)
          ls_value <- LSKNW_K2Ep_rcpp(Y = Y, X = X %*% beta, h = bandwidth) * number_n +
            (rho.augLagrange.penalty / 2) *
            sum((theta - gamma_free_iter + u_free_iter) ^ 2)

          return(ls_value)
        }

        SemiLS.penalized.esti <- nlminb(start = beta_free_iter,
                                        objective = SemiLS.penalized)

        beta_free_iter_new <- SemiLS.penalized.esti$par

        gamma_free_iter_new <- sign(beta_free_iter + u_free_iter) * pmax(
          abs(beta_free_iter + u_free_iter) - weight.free.beta / rho.augLagrange.penalty,
          0)
        u_free_iter_new <- u_free_iter + beta_free_iter_new - gamma_free_iter_new

        res_primal <- sqrt(sum((beta_free_iter_new - gamma_free_iter_new) ^ 2))
        res_dual <- rho.augLagrange.penalty *
          sqrt(sum((gamma_free_iter_new - gamma_free_iter) ^ 2))

        if ((res_primal > tol) | (res_dual > tol))
        {
          beta_free_iter <- beta_free_iter_new
          gamma_free_iter <- gamma_free_iter_new
          u_free_iter <- u_free_iter_new

          print(c(paste("Iteration", iter, " par = ", sep = ""), beta_free_iter))
        }else
        {
          break
        }
      }

      results <- list(coef = c(1, beta_free_iter),
                      iteration = iter)
    }

    if (kernel == "K2.Biweight")
    {

    }

    if (kernel == "K4.Biweight")
    {

    }
  }

  if (type == "distribution")
  {
    if (kernel == "K2.Epanechnikov")
    {

    }

    if (kernel == "K2.Biweight")
    {

    }

    if (kernel == "K4.Biweight")
    {

    }
  }


}

SemiMI <- function(data = NULL, X.name = NULL, Y.name = NULL,
                   X = NULL, Y = NULL,
                   type = "mean", dimension = 2, kernel = "K2.Biweight",
                   bandwidth = NULL, bandwidth.initial = 1,
                   beta.initial = NULL,
                   wi.boot = NULL)
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

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (number_p <= dimension)
  {
    warnings("The number of covariates should be larger than the number of indices.")
  }

  number_k <- (number_p - dimension) * dimension

  if (is.null(beta.initial))
  {
    beta.initial <- matrix(0,
                           nrow = number_p - dimension,
                           ncol = dimension)
  }else
  {
    beta.initial <- matrix(beta.initial,
                           nrow = number_p,
                           ncol = dimension)
    beta.initial <- matrix(beta.initial[(dimension + 1):number_p, ],
                           nrow = number_p - dimension,
                           ncol = dimension)
  }

  if (type == "mean")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNW_K2Ep_rcpp(Y = Y, X = X_beta, h = h * sd_X_beta)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNW_K2Ep_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNW_K2Ep_w_rcpp(Y = Y, X = X_beta,
                                          h = h * sd_X_beta, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNW_K2Ep_w_rcpp(Y = Y, X = X %*% beta,
                                        h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNW_K2Bw_rcpp(Y = Y, X = X_beta, h = h * sd_X_beta)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNW_K2Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNW_K2Bw_w_rcpp(Y = Y, X = X_beta,
                                          h = h * sd_X_beta, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNW_K2Bw_w_rcpp(Y = Y, X = X %*% beta,
                                        h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNW_K4Bw_rcpp(Y = Y, X = X_beta, h = h * sd_X_beta)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNW_K4Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNW_K4Bw_w_rcpp(Y = Y, X = X_beta,
                                          h = h * sd_X_beta, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNW_K4Bw_w_rcpp(Y = Y, X = X %*% beta,
                                        h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }
    }
  }

  if (type == "distribution")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNWcdf_K2Ep_rcpp(Y = Y, X = X_beta, h = h * sd_X_beta)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNWcdf_K2Ep_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNWcdf_K2Ep_w_rcpp(Y = Y, X = X_beta,
                                             h = h * sd_X_beta, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNWcdf_K2Ep_w_rcpp(Y = Y, X = X %*% beta,
                                           h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNWcdf_K2Bw_rcpp(Y = Y, X = X_beta, h = h * sd_X_beta)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNWcdf_K2Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNWcdf_K2Bw_w_rcpp(Y = Y, X = X_beta,
                                             h = h * sd_X_beta, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNWcdf_K2Bw_w_rcpp(Y = Y, X = X %*% beta,
                                           h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(wi.boot))
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNWcdf_K4Bw_rcpp(Y = Y, X = X_beta, h = h * sd_X_beta)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNWcdf_K4Bw_rcpp(Y = Y, X = X %*% beta, h = bandwidth)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }else
      {
        if (is.null(bandwidth))
        {
          LOOCV <- function(theta)
          {
            beta <- rbind(diag(dimension),
                          matrix(theta[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension))
            h <- exp(theta[number_k + 1])
            X_beta <- X %*% beta
            sd_X_beta <- apply(X_beta, 2, sd)
            cv_value <- CVKNWcdf_K4Bw_w_rcpp(Y = Y, X = X_beta,
                                             h = h * sd_X_beta, w = wi.boot)

            return(cv_value)
          }

          LOOCV.esti <- nlminb(start = c(as.vector(beta.initial),
                                         log(bandwidth.initial)),
                               objective = LOOCV)
          beta.CV <- rbind(diag(dimension),
                           matrix(LOOCV.esti$par[1:number_k],
                                  nrow = number_p - dimension,
                                  ncol = dimension))
          bandwidth <- exp(LOOCV.esti$par[number_k + 1]) *
            apply(X %*% beta.CV, 2, sd)

          results <- list(coef.CV = beta.CV,
                          bandwidth = bandwidth,
                          CV.value = LOOCV.esti$objective)

          beta.initial <- matrix(LOOCV.esti$par[1:number_k],
                                 nrow = number_p - dimension,
                                 ncol = dimension)
        }else
        {
          bandwidth <- rep(bandwidth, length = dimension)
        }

        SemiLS <- function(theta)
        {
          beta <- rbind(diag(dimension),
                        matrix(theta,
                               nrow = number_p - dimension,
                               ncol = dimension))
          ls_value <- LSKNWcdf_K4Bw_w_rcpp(Y = Y, X = X %*% beta,
                                           h = bandwidth, w = wi.boot)

          return(ls_value)
        }

        SemiLS.esti <- nlminb(start = as.vector(beta.initial),
                              objective = SemiLS)

        results$coef <- rbind(diag(dimension),
                              matrix(SemiLS.esti$par,
                                     nrow = number_p - dimension,
                                     ncol = dimension))
        results$LS.value <- SemiLS.esti$objective
      }
    }
  }

  return(results)
}
