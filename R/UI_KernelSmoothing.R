KernelSmooth <- function(data = NULL, X.name = NULL, Y.name = NULL,
                         X = NULL, Y = NULL, x = NULL, y = NULL,
                         type = "mean", kernel = "K2.Biweight",
                         bandwidth = NULL, bandwidth.initial = 1,
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

  if (is.null(x))
  {
    x <- X
  }else
  {
    x <- as.matrix(x)
  }

  if (is.null(y))
  {
    y <- sort(unique(as.vector(Y)))
  }else
  {
    y <- as.vector(y)
  }

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_l <- length(y)

  if (type == "mean")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNW_K2Ep_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        mhat <- KNW_K2Ep_rcpp(Y = Y, X = X, x = x, h = bandwidth)
      }else
      {
        mhat <- KNW_K2Ep_w_rcpp(Y = Y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(mean = mhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNW_K2Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        mhat <- KNW_K2Bw_rcpp(Y = Y, X = X, x = x, h = bandwidth)
      }else
      {
        mhat <- KNW_K2Bw_w_rcpp(Y = Y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(mean = mhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNW_K4Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        mhat <- KNW_K4Bw_rcpp(Y = Y, X = X, x = x, h = bandwidth)
      }else
      {
        mhat <- KNW_K4Bw_w_rcpp(Y = Y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(mhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(mean = mhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }
  }

  if (type == "distribution")
  {
    if (kernel == "K2.Epanechnikov")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNWcdf_K2Ep_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        Fhat <- KNWcdf_K2Ep_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth)
      }else
      {
        Fhat <- KNWcdf_K2Ep_w_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(Fhat) <- y
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(distribution = Fhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K2.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNWcdf_K2Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        Fhat <- KNWcdf_K2Bw_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth)
      }else
      {
        Fhat <- KNWcdf_K2Bw_w_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(Fhat) <- y
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(distribution = Fhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }

    if (kernel == "K4.Biweight")
    {
      if (is.null(bandwidth))
      {
        LOOCV <- function(h.log)
        {
          cv_value <- CVKNWcdf_K4Bw_rcpp(Y = Y, X = X, h = rep(exp(h.log), number_p))

          return(cv_value)
        }

        LOOCV.esti <- nlminb(start = log(bandwidth.initial),
                             objective = LOOCV)
        bandwidth <- rep(exp(LOOCV.esti$par), number_p)
      }else
      {
        bandwidth <- rep(bandwidth, length.out = number_p)
      }

      if (is.null(wi.boot))
      {
        Fhat <- KNWcdf_K4Bw_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth)
      }else
      {
        Fhat <- KNWcdf_K4Bw_w_rcpp(Y = Y, y = y, X = X, x = x, h = bandwidth, w = wi.boot)
      }

      rownames(Fhat) <- rownames(x) <- paste("subject", 1:number_k, sep = ".")
      colnames(Fhat) <- y
      colnames(x) <- names(bandwidth) <- paste("covariate", 1:number_p, sep = ".")
      results <- list(distribution = Fhat,
                      eval.covariate = x,
                      bandwidth = bandwidth)
    }
  }

  return(results)
}

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








