GLMcombineADavar <- function(name.par,
                             MLE.hessian,
                             Psi.square,
                             Psi.diff.beta,
                             Psi.diff.phi = NULL,
                             Psi.diff.theta = NULL,
                             kappa = NULL,
                             avar.phi = NULL,
                             avar.working = NULL,
                             eps.inv = 1e-7)
{
  number_l <- dim(MLE.hessian)[1]
  number_m <- dim(Psi.square)[1]
  M.name <- name.par
  if (!is.null(kappa))
  {
    M.name <- c(M.name, paste("phi", 1:number_m, sep = ""))
  }
  if (!is.null(Psi.diff.theta))
  {
    number_theta <- length(Psi.diff.theta) / number_m
    M.name <- c(M.name, paste("theta", 1:number_theta, sep = ""))
  }
  M.name <- c(M.name, paste("eta", 1:number_m, sep = ""))
  number_total <- length(M.name)

  avar.S <- matrix(0, nrow = number_total, ncol = number_total)
  J.V <- matrix(0, nrow = number_total, ncol = number_total)
  dimnames(avar.S) <- list(M.name, M.name)
  dimnames(J.V) <- list(M.name, M.name)

  avar.S[name.par, name.par] <- -MLE.hessian
  avar.S[paste("eta", 1:number_m, sep = ""),
         paste("eta", 1:number_m, sep = "")] <- Psi.square

  J.V[name.par, name.par] <- -MLE.hessian
  J.V[paste("eta", 1:number_m, sep = ""),
      paste("eta", 1:number_m, sep = "")] <- -Psi.square
  J.V[paste("eta", 1:number_m, sep = ""),
      name.par] <- as.vector(Psi.diff.beta)
  J.V[name.par,
      paste("eta", 1:number_m, sep = "")] <- as.vector(t(Psi.diff.beta))

  if (!is.null(kappa))
  {
    if (is.null(avar.working))
    {
      avar.inv <- inv_sympd_rcpp(avar.phi + eps.inv * diag(number_m))
      avar.S[paste("phi", 1:number_m, sep = ""),
             paste("phi", 1:number_m, sep = "")] <- kappa * avar.inv
      J.V[paste("phi", 1:number_m, sep = ""),
          paste("phi", 1:number_m, sep = "")] <- kappa * avar.inv
    }else
    {
      avar.working.inv <- inv_sympd_rcpp(avar.working + eps.inv * diag(number_m))
      avar.S[paste("phi", 1:number_m, sep = ""),
             paste("phi", 1:number_m, sep = "")] <- kappa * avar.working.inv %*%
        avar.phi %*% avar.working.inv
      J.V[paste("phi", 1:number_m, sep = ""),
          paste("phi", 1:number_m, sep = "")] <- kappa * avar.working.inv
    }

    J.V[paste("eta", 1:number_m, sep = ""),
        paste("phi", 1:number_m, sep = "")] <- as.vector(Psi.diff.phi)
    J.V[paste("phi", 1:number_m, sep = ""),
        paste("eta", 1:number_m, sep = "")] <- as.vector(t(Psi.diff.phi))
  }

  if (!is.null(Psi.diff.theta))
  {
    J.V[paste("eta", 1:number_m, sep = ""),
        paste("theta", 1:number_theta, sep = "")] <- as.vector(Psi.diff.theta)
    J.V[paste("theta", 1:number_theta, sep = ""),
        paste("eta", 1:number_m, sep = "")] <- as.vector(t(Psi.diff.theta))
  }

  J.V.inv <- solve_rcpp(J.V, diag(number_total))
  asy.Cov <- J.V.inv %*% avar.S %*% J.V.inv
  dimnames(asy.Cov) <- list(M.name, M.name)
  asy.Cov.par <- asy.Cov[name.par, name.par]

  results <- list(avar.S = avar.S,
                  J.V = J.V,
                  J.V.inv = J.V.inv,
                  asy.Cov = asy.Cov,
                  asy.Cov.par = asy.Cov.par)
  return(results)
}
