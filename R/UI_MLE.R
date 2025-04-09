GLM.MLE <- function(data = NULL, X.name = NULL, Y.name = NULL,
                    X = NULL, Y = NULL,
                    distribution = "normal", N.binomial = NULL,
                    initial = NULL, wi.boot = NULL, do.SE = TRUE,
                    X.future = NULL)
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

  if (distribution == "normal")
  {
    if (is.null(initial))
    {
      initial.trans <- c(0, rep(0, number_p), 0)
    }else
    {
      initial.trans <- initial
      initial.trans[number_p + 2] <- log(initial[number_p + 2])
    }

    if (is.null(wi.boot))
    {
      nlL.run <- function(theta.trans)
      {
        value <- -lL_normal_rcpp(X = X, Y = Y,
                                 alpha = theta.trans[1],
                                 beta = theta.trans[2:(number_p + 1)],
                                 sigma = exp(theta.trans[number_p + 2]))

        return(value)
      }
    }else
    {
      nlL.run <- function(theta.trans)
      {
        value <- -lL_normal_w_rcpp(X = X, Y = Y,
                                   alpha = theta.trans[1],
                                   beta = theta.trans[2:(number_p + 1)],
                                   sigma = exp(theta.trans[number_p + 2]),
                                   w = wi.boot)

        return(value)
      }
    }

    theta.trans.hat <- nlminb(start = initial.trans,
                              objective = nlL.run)$par
    theta.hat <- theta.trans.hat
    theta.hat[number_p + 2] <- exp(theta.trans.hat[number_p + 2])
    names(theta.hat) <- c("alpha", paste("beta", 1:number_p, sep = ""), "sigma")

    results <- list(alpha = theta.hat[1],
                    beta = theta.hat[2:(number_p + 1)],
                    sigma = theta.hat[number_p + 2],
                    parameter = theta.hat)

    if (do.SE)
    {
      asy.Cov <- inv_sympd_rcpp(-diff_lL_normal_rcpp(
        X = X, Y = Y,
        alpha = theta.hat[1],
        beta = theta.hat[2:(number_p + 1)],
        sigma = theta.hat[number_p + 2])$hessian)
      dimnames(asy.Cov) <- list(
        c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"),
        c("alpha", paste("beta", 1:number_p, sep = ""), "sigma"))

      results$Cov.coef <- asy.Cov / number_n
    }

    if (!is.null(X.future))
    {
      X.future <- as.matrix(X.future)

      Y.predict <- theta.hat[1] + X.future %*% theta.hat[2:(number_p + 1)]

      results$Y.predict <- as.vector(Y.predict)
    }
  }

  if (distribution == "Bernoulli")
  {
    if (is.null(initial))
    {
      initial <- c(0, rep(0, number_p))
    }

    if (is.null(wi.boot))
    {
      nlL.run <- function(theta)
      {
        value <- -lL_logistic_rcpp(X = X, Y = Y,
                                   alpha = theta[1],
                                   beta = theta[2:(number_p + 1)])

        return(value)
      }
    }else
    {
      nlL.run <- function(theta)
      {
        value <- -lL_logistic_w_rcpp(X = X, Y = Y,
                                     alpha = theta[1],
                                     beta = theta[2:(number_p + 1)],
                                     w = wi.boot)

        return(value)
      }
    }

    theta.hat <- nlminb(start = initial,
                        objective = nlL.run)$par
    names(theta.hat) <- c("alpha", paste("beta", 1:number_p, sep = ""))

    results <- list(alpha = theta.hat[1],
                    beta = theta.hat[2:(number_p + 1)],
                    parameter = theta.hat)

    if (do.SE)
    {
      asy.Cov <- inv_sympd_rcpp(-diff_lL_logistic_rcpp(
        X = X, Y = Y,
        alpha = theta.hat[1],
        beta = theta.hat[2:(number_p + 1)])$hessian)
      dimnames(asy.Cov) <- list(
        c("alpha", paste("beta", 1:number_p, sep = "")),
        c("alpha", paste("beta", 1:number_p, sep = "")))

      results$Cov.coef <- asy.Cov / number_n
    }

    if (!is.null(X.future))
    {
      X.future <- as.matrix(X.future)

      SI.future <- as.vector(
        theta.hat[1] + X.future %*% theta.hat[2:(number_p + 1)])
      P1.future <- 1 / (exp(-SI.future) + 1)
      Y.predict <- (P1.future >= 0.5) * 1

      results$Y.predict <- as.vector(Y.predict)
      results$Y.posterior <- as.vector(P1.future)
    }
  }

  if (distribution == "Gamma")
  {
    if (is.null(initial))
    {
      initial.trans <- c(0, rep(0, number_p), 0)
    }else
    {
      initial.trans <- initial
      initial.trans[number_p + 2] <- log(initial[number_p + 2])
    }

    if (is.null(wi.boot))
    {
      nlL.run <- function(theta.trans)
      {
        value <- -lL_Gamma_rcpp(
          X = X, Y = Y,
          alpha = theta.trans[1],
          beta = theta.trans[2:(number_p + 1)],
          nu = exp(theta.trans[number_p + 2]))

        return(value)
      }
    }else
    {
      nlL.run <- function(theta.trans)
      {
        value <- -lL_Gamma_w_rcpp(
          X = X, Y = Y,
          alpha = theta.trans[1],
          beta = theta.trans[2:(number_p + 1)],
          nu = exp(theta.trans[number_p + 2]),
          w = wi.boot)

        return(value)
      }
    }

    theta.trans.hat <- nlminb(start = initial.trans,
                              objective = nlL.run)$par
    theta.hat <- theta.trans.hat
    theta.hat[number_p + 2] <- exp(theta.trans.hat[number_p + 2])
    names(theta.hat) <- c("alpha", paste("beta", 1:number_p, sep = ""), "shape")

    results <- list(alpha = theta.hat[1],
                    beta = theta.hat[2:(number_p + 1)],
                    nu = theta.hat[number_p + 2],
                    parameter = theta.hat)
    if (do.SE)
    {
      asy.Cov <- inv_sympd_rcpp(-diff_lL_Gamma_rcpp(
        X = X, Y = Y,
        alpha = theta.hat[1],
        beta = theta.hat[2:(number_p + 1)],
        nu = theta.hat[number_p + 2])$hessian)
      dimnames(asy.Cov) <- list(
        c("alpha", paste("beta", 1:number_p, sep = ""), "shape"),
        c("alpha", paste("beta", 1:number_p, sep = ""), "shape"))

      results$Cov.coef <- asy.Cov / number_n
    }

    if (!is.null(X.future))
    {
      X.future <- as.matrix(X.future)

      Y.predict <- theta.hat[number_p + 2] / exp(
        theta.hat[1] + X.future %*% theta.hat[2:(number_p + 1)])

      results$Y.predict <- as.vector(Y.predict)
    }
  }

  if (distribution == "binomial")
  {
    if (is.null(N.binomial))
    {
      N.binomial <- max(Y)
    }

    if (is.null(initial))
    {
      initial <- c(0, rep(0, number_p))
    }

    if (is.null(wi.boot))
    {
      nlL.run <- function(theta)
      {
        value <- -lL_binomial_rcpp(X = X, Y = Y,
                                   alpha = theta[1],
                                   beta = theta[2:(number_p + 1)],
                                   N = N.binomial)

        return(value)
      }
    }else
    {
      nlL.run <- function(theta)
      {
        value <- -lL_binomial_w_rcpp(X = X, Y = Y,
                                     alpha = theta[1],
                                     beta = theta[2:(number_p + 1)],
                                     N = N.binomial,
                                     w = wi.boot)

        return(value)
      }
    }

    theta.hat <- nlminb(start = initial,
                        objective = nlL.run)$par
    names(theta.hat) <- c("alpha", paste("beta", 1:number_p, sep = ""))

    results <- list(alpha = theta.hat[1],
                    beta = theta.hat[2:(number_p + 1)],
                    parameter = theta.hat)

    if (do.SE)
    {
      asy.Cov <- inv_sympd_rcpp(-diff_lL_binomial_rcpp(
        X = X, Y = Y,
        alpha = theta.hat[1],
        beta = theta.hat[2:(number_p + 1)],
        N = N.binomial)$hessian)
      dimnames(asy.Cov) <- list(
        c("alpha", paste("beta", 1:number_p, sep = "")),
        c("alpha", paste("beta", 1:number_p, sep = "")))

      results$Cov.coef <- asy.Cov / number_n
    }

    if (!is.null(X.future))
    {
      X.future <- as.matrix(X.future)

      SI.future <- as.vector(
        theta.hat[1] + X.future %*% theta.hat[2:(number_p + 1)])
      Y.predict <- N.binomial / (exp(-SI.future) + 1)

      results$Y.predict <- as.vector(Y.predict)
    }
  }

  if (distribution == "Poisson")
  {
    if (is.null(initial))
    {
      initial <- c(0, rep(0, number_p))
    }

    if (is.null(wi.boot))
    {
      nlL.run <- function(theta)
      {
        value <- -lL_Poisson_rcpp(X = X, Y = Y,
                                  alpha = theta[1],
                                  beta = theta[2:(number_p + 1)])

        return(value)
      }
    }else
    {
      nlL.run <- function(theta)
      {
        value <- -lL_Poisson_w_rcpp(X = X, Y = Y,
                                    alpha = theta[1],
                                    beta = theta[2:(number_p + 1)],
                                    w = wi.boot)

        return(value)
      }
    }

    theta.hat <- nlminb(start = initial,
                        objective = nlL.run)$par
    names(theta.hat) <- c("alpha", paste("beta", 1:number_p, sep = ""))

    results <- list(alpha = theta.hat[1],
                    beta = theta.hat[2:(number_p + 1)],
                    parameter = theta.hat)

    if (do.SE)
    {
      asy.Cov <- inv_sympd_rcpp(-diff_lL_Poisson_rcpp(
        X = X, Y = Y,
        alpha = theta.hat[1],
        beta = theta.hat[2:(number_p + 1)])$hessian)
      dimnames(asy.Cov) <- list(
        c("alpha", paste("beta", 1:number_p, sep = "")),
        c("alpha", paste("beta", 1:number_p, sep = "")))

      results$Cov.coef <- asy.Cov / number_n
    }

    if (!is.null(X.future))
    {
      X.future <- as.matrix(X.future)

      SI.future <- as.vector(
        theta.hat[1] + X.future %*% theta.hat[2:(number_p + 1)])
      Y.predict <- exp(SI.future)

      results$Y.predict <- as.vector(Y.predict)
    }
  }

  return(results)
}
