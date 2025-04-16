##### No shift, normal regression model

AD_EY_normal_Lagrange <- function(X, alpha, beta, sigma, phi, eta)
{
  n_n <- dim(X)[1]

  Psi <- alpha + X %*% beta - phi
  denominator = 1 + eta * Psi
  value <- sum(log(denominator))
  gradient <- sum(Psi / denominator)
  hessian <- -sum((Psi / denominator) ^ 2)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- 0
    hessian <- 0
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

AD_EY_normal_SolveLagrange <- function(X, alpha, beta, sigma, phi,
                                       eta.initial, iter.max,
                                       step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EY_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, eta = eta
  )

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EY_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EY_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EY_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma, phi,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EY_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, eta = eta
  )

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EY_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EY_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EY_normal <- function(X, alpha, beta, sigma, phi)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  Psi_i <- alpha + X %*% beta - phi
  Psi <- sum(Psi_i) / n_n
  Psi_square <- sum(Psi_i ^ 2) / n_n
  Psi_gradient <- c(1, colSums(X) / n_n, 0, -1)

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- 0
    Psi_square <- 0
    Psi_gradient <- rep(0, n_p + 3)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADvar_EY_normal <- function(X, alpha, beta, sigma, phi, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  Psi_i <- as.vector(alpha + X %*% beta - phi)
  p_i <- (1 / (1 + Psi_i * eta)) / n_n

  var.A <- -1
  var.B <- sum((Psi_i ^ 2) * p_i) + sigma ^ 2

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

AD_EXsubY_normal_Lagrange <- function(X, alpha, beta, sigma, phi, y.pts, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  SI <- as.vector(alpha + X %*% beta)
  dist.cdf <- pnorm(outer_minus_rcpp(y.pts[, 2], SI) / sigma, mean = 0, sd = 1) -
    pnorm(outer_minus_rcpp(y.pts[, 1], SI) / sigma, mean = 0, sd = 1)
  Psi_i <- aperm(array(
    X[rep(1:n_n, times = n_k), ] - phi[rep(1:n_k, each = n_n), ],
    c(n_n, n_k, n_p)), c(2, 1, 3)) * as.vector(dist.cdf)

  denominator <- as.vector(1 + matrix(aperm(Psi_i, c(2, 3, 1)),
                                      nrow = n_n, ncol = n_m) %*% eta)
  value <- sum(log(denominator))
  Psi_i_scaled <- matrix(aperm(Psi_i, c(2, 3, 1)) / denominator,
                         nrow = n_n, ncol = n_m)
  gradient <- as.matrix(colSums(Psi_i_scaled))
  hessian <- -matrix(colSums(Psi_i_scaled[, rep(1:n_m, times = n_m)] *
    Psi_i_scaled[, rep(1:n_m, each = n_m)]), nrow = n_m, ncol = n_m)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_m, 1)
    hessian <- matrix(0, n_m, n_m)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

AD_EXsubY_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                           phi, y.pts,
                                           eta.initial, iter.max,
                                           step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EXsubY_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, y.pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EXsubY_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, y.pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EXsubY_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, y.pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EXsubY_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                              phi, y.pts,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EXsubY_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EXsubY_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EXsubY_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, y_pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EXsubY_normal <- function(X, alpha, beta, sigma, phi, y.pts)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  Psi_gradient <- matrix(0, nrow = n_m, ncol = n_p + n_m + 2)

  SI <- as.vector(alpha + X %*% beta)
  a.pts.scaled <- outer_minus_rcpp(y.pts[, 1], SI) / sigma
  b.pts.scaled <- outer_minus_rcpp(y.pts[, 2], SI) / sigma
  dist.cdf <- pnorm(b.pts.scaled, mean = 0, sd = 1) -
    pnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf <- dnorm(b.pts.scaled, mean = 0, sd = 1) -
    dnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf.w <- dnorm(b.pts.scaled, mean = 0, sd = 1) * b.pts.scaled -
    dnorm(a.pts.scaled, mean = 0, sd = 1) * a.pts.scaled
  dist.x.phi <- array(
    X[rep(1:n_n, times = n_k), ] - phi[rep(1:n_k, each = n_n), ],
    c(n_n, n_k, n_p))
  Psi_i <- matrix(aperm(aperm(dist.x.phi, c(2, 1, 3)) * as.vector(dist.cdf), c(3, 1, 2)),
    nrow = n_m, ncol = n_n)

  Psi <- as.matrix(rowSums(Psi_i) / n_n)
  Psi_square <- matrix(rowSums(Psi_i[rep(1:n_m, times = n_m), ] *
                                 Psi_i[rep(1:n_m, each = n_m), ]) / n_n,
                       nrow = n_m, ncol = n_m)
  Psi_gradient[, 1:(n_p + 1)] <- matrix(apply(-aperm(array(matrix(aperm(
    dist.x.phi, c(3, 2, 1)),
    nrow = n_m, ncol = n_n)[rep(1:n_m, times = n_p + 1), ] *
      t(cbind(1, X))[rep(1:(n_p + 1), each = n_m), ],
    c(n_p, n_k, n_p + 1, n_n)),
    c(2, 4, 1, 3)) * as.vector(dist.pdf) / sigma, c(3, 1, 4), sum) / n_n,
    nrow = n_m, ncol = n_p + 1)
  Psi_gradient[, n_p + 2] <- rowSums(matrix(aperm(-aperm(
    dist.x.phi, c(2, 1, 3)) * as.vector(dist.pdf.w) / sigma,
    c(3, 1, 2)), nrow = n_m, ncol = n_n)) / n_n
  Psi_gradient[, (n_p + 3):(n_p + n_m + 2)] <- diag(
    rep(-rowSums(dist.cdf) / n_n, each = n_p))

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- matrix(0, n_m, 1)
    Psi_square <- matrix(0, n_m, n_m)
    Psi_gradient <- matrix(0, nrow = n_m, ncol = n_p + n_m + 2)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADvar_EXsubY_normal <- function(X, alpha, beta, sigma, phi, y.pts, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  SI <- as.vector(alpha + X %*% beta)
  a.pts.scaled <- outer_minus_rcpp(y.pts[, 1], SI) / sigma
  b.pts.scaled <- outer_minus_rcpp(y.pts[, 2], SI) / sigma
  dist.cdf <- pnorm(b.pts.scaled, mean = 0, sd = 1) -
    pnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf <- dnorm(b.pts.scaled, mean = 0, sd = 1) -
    dnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf.w <- dnorm(b.pts.scaled, mean = 0, sd = 1) * b.pts.scaled -
    dnorm(a.pts.scaled, mean = 0, sd = 1) * a.pts.scaled
  dist.x.phi <- array(
    X[rep(1:n_n, times = n_k), ] - phi[rep(1:n_k, each = n_n), ],
    c(n_n, n_k, n_p))
  Psi_i <- matrix(aperm(aperm(dist.x.phi, c(2, 1, 3)) * as.vector(dist.cdf),
                        c(3, 1, 2)), nrow = n_m, ncol = n_n)
  p_i <- (1 / (1 + as.vector(t(Psi_i) %*% eta))) / n_n

  var.A <- -diag(rep(colSums(t(dist.cdf) * p_i), each = n_k),
                 nrow = n_m, ncol = n_m)
  var.B.array <- array(colSums(
    dist.x.phi[, , rep(1:n_p, times = n_p)] *
      dist.x.phi[, , rep(1:n_p, each = n_p)] * as.vector(
        t(dist.cdf)) * p_i), c(n_k, n_p, n_p))
  var.B <- matrix(0, n_m, n_m)
  for (k in 1:n_k)
  {
    var.B[(1 + (k - 1) * n_p):(k * n_p),
          (1 + (k - 1) * n_p):(k * n_p)] <- var.B.array[k, , ]
  }

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

AD_EYsubX_normal_Lagrange <- function(X, alpha, beta, sigma, phi, inclusion, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  Psi <- outer_minus_rcpp(as.vector(alpha + X %*% beta), phi) * inclusion
  denominator = as.vector(1 + Psi %*% eta)
  Psi_scaled <- Psi / denominator
  value <- sum(log(denominator))
  gradient <- as.matrix(colSums(Psi_scaled))
  hessian <- matrix(-colSums(Psi_scaled[, rep(1:n_m, times = n_m)] *
                               Psi_scaled[, rep(1:n_m, each = n_m)]),
                    nrow = n_m, ncol = n_m)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_m, 1)
    hessian <- matrix(0, n_m, n_m)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

AD_EYsubX_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                           phi, inclusion,
                                           eta.initial, iter.max,
                                           step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EYsubX_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EYsubX_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EYsubX_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EYsubX_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                              phi, inclusion,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EYsubX_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EYsubX_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EYsubX_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EYsubX_normal <- function(X, alpha, beta, sigma, phi, inclusion)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  Psi_i <- outer_minus_rcpp(as.vector(alpha + X %*% beta), phi) * inclusion
  Psi <- as.matrix(colSums(Psi_i) / n_n)
  Psi_square <- matrix(colSums(Psi_i[, rep(1:n_m, times = n_m)] *
                                 Psi_i[, rep(1:n_m, each = n_m)]) / n_n,
                       nrow = n_m, ncol = n_m)
  Psi_gradient <- t(matrix(colSums(
    cbind(1, X, 0)[, rep(1:(n_p + 2), times = n_m)] *
      inclusion[, rep(1:n_m, each = n_p + 2)]) / n_n,
    nrow = n_p + 2, ncol = n_m))
  Psi_gradient <- cbind(Psi_gradient,
                        diag(-colSums(inclusion) / n_n, n_m, n_m))

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- matrix(0, n_m, 1)
    Psi_square <- matrix(0, n_m, n_m)
    Psi_gradient <- matrix(0, n_m, n_p + n_m + 2)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADvar_EYsubX_normal <- function(X, alpha, beta, sigma, phi, inclusion, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  Psi_i <- outer_minus_rcpp(as.vector(alpha + X %*% beta), phi) * inclusion
  p_i <- (1 / (1 + as.vector(Psi_i %*% eta))) / n_n

  var.A <- diag(-colSums(inclusion * p_i), nrow = n_m, ncol = n_m)
  var.B <- colSums(Psi_i[, rep(1:n_m, times = n_m)] *
                     Psi_i[, rep(1:n_m, each = n_m)] * p_i) +
    (sigma ^ 2) * matrix(colSums(
      inclusion[, rep(1:n_m, times = n_m)] *
        inclusion[, rep(1:n_m, each = n_m)] * p_i), nrow = n_m, ncol = n_m)

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

##### Covariate shift, normal regression model

ADCS_EY_normal_Lagrange <- function(X, alpha, beta, sigma, phi, CS.beta, eta)
{
  n_n <- dim(X)[1]

  Psi <- (alpha + X %*% beta - phi) * as.vector(exp(X %*% CS.beta))
  denominator = 1 + eta * Psi
  value <- sum(log(denominator))
  gradient <- sum(Psi / denominator)
  hessian <- -sum((Psi / denominator) ^ 2)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- 0
    hessian <- 0
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

ADCS_EY_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                         phi, CS.beta,
                                         eta.initial, iter.max,
                                         step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADCS_EY_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS.beta = CS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADCS_EY_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS.beta = CS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EY_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS.beta = CS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EY_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                            phi, CS.beta,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADCS_EY_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS_beta = CS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADCS_EY_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS_beta = CS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EY_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS_beta = CS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EY_normal <- function(X, alpha, beta, sigma, phi, CS.beta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  eSI_CS <- as.vector(exp(X %*% CS.beta))

  Psi_i <- (alpha + X %*% beta - phi) * eSI_CS
  Psi <- sum(Psi_i) / n_n
  Psi_square <- sum(Psi_i ^ 2) / n_n
  Psi_gradient <- colSums(cbind(
    eSI_CS, X * eSI_CS, 0, -eSI_CS, X * as.vector(Psi_i))) / n_n

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- 0
    Psi_square <- 0
    Psi_gradient <- rep(0, n_p * 2 + 3)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADCSvar_EY_normal <- function(X, alpha, beta, sigma, phi, CS.beta, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  eSI_CS <- as.vector(exp(X %*% CS.beta))

  Psi_i <- (alpha + X %*% beta - phi) * eSI_CS
  p_i <- (1 / (1 + Psi_i * eta)) / n_n

  var.A <- -1
  var.B <- sum((alpha + X %*% beta - phi) ^ 2 * eSI_CS * p_i) /
    sum(eSI_CS * p_i) + sigma ^ 2

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

ADCS_EXsubY_normal_Lagrange <- function(X, alpha, beta, sigma, phi,
                                        CS.beta, y.pts, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  SI <- as.vector(alpha + X %*% beta)
  eSI_CS <- as.vector(exp(X %*% CS.beta))
  dist.cdf <- pnorm(outer_minus_rcpp(y.pts[, 2], SI) / sigma, mean = 0, sd = 1) -
    pnorm(outer_minus_rcpp(y.pts[, 1], SI) / sigma, mean = 0, sd = 1)
  Psi_i <- aperm(array(
    X[rep(1:n_n, times = n_k), ] - phi[rep(1:n_k, each = n_n), ],
    c(n_n, n_k, n_p)) * eSI_CS, c(2, 1, 3)) * as.vector(dist.cdf)

  denominator <- as.vector(1 + matrix(aperm(Psi_i, c(2, 3, 1)),
                                      nrow = n_n, ncol = n_m) %*% eta)
  value <- sum(log(denominator))
  Psi_i_scaled <- matrix(aperm(Psi_i, c(2, 3, 1)) / denominator,
                         nrow = n_n, ncol = n_m)
  gradient <- as.matrix(colSums(Psi_i_scaled))
  hessian <- -matrix(colSums(Psi_i_scaled[, rep(1:n_m, times = n_m)] *
                               Psi_i_scaled[, rep(1:n_m, each = n_m)]), nrow = n_m, ncol = n_m)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_m, 1)
    hessian <- matrix(0, n_m, n_m)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

ADCS_EXsubY_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                             phi, CS.beta, y.pts,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADCS_EXsubY_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS.beta = CS.beta, y.pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADCS_EXsubY_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS.beta = CS.beta, y.pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EXsubY_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS.beta = CS.beta, y.pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EXsubY_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                                phi, CS.beta, y.pts,
                                                eta.initial, iter.max,
                                                step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADCS_EXsubY_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADCS_EXsubY_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EXsubY_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EXsubY_normal <- function(X, alpha, beta, sigma, phi, CS.beta, y.pts)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  Psi_gradient <- matrix(0, nrow = n_m, ncol = n_p * 2 + n_m + 2)

  SI <- as.vector(alpha + X %*% beta)
  eSI_CS <- as.vector(exp(X %*% CS.beta))
  a.pts.scaled <- outer_minus_rcpp(y.pts[, 1], SI) / sigma
  b.pts.scaled <- outer_minus_rcpp(y.pts[, 2], SI) / sigma
  dist.cdf <- pnorm(b.pts.scaled, mean = 0, sd = 1) -
    pnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf <- dnorm(b.pts.scaled, mean = 0, sd = 1) -
    dnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf.w <- dnorm(b.pts.scaled, mean = 0, sd = 1) * b.pts.scaled -
    dnorm(a.pts.scaled, mean = 0, sd = 1) * a.pts.scaled
  dist.x.phi <- array(
    X[rep(1:n_n, times = n_k), ] - phi[rep(1:n_k, each = n_n), ],
    c(n_n, n_k, n_p))
  Psi_i <- t(t(matrix(aperm(aperm(
    dist.x.phi, c(2, 1, 3)) * as.vector(dist.cdf), c(3, 1, 2)),
    nrow = n_m, ncol = n_n)) * eSI_CS)

  Psi <- as.matrix(rowSums(Psi_i) / n_n)
  Psi_square <- matrix(rowSums(Psi_i[rep(1:n_m, times = n_m), ] *
                                 Psi_i[rep(1:n_m, each = n_m), ]) / n_n,
                       nrow = n_m, ncol = n_m)
  Psi_gradient[, 1:(n_p + 1)] <- matrix(apply(-aperm(array(t(t(matrix(aperm(
    dist.x.phi, c(3, 2, 1)),
    nrow = n_m, ncol = n_n)) * eSI_CS)[rep(1:n_m, times = n_p + 1), ] *
      t(cbind(1, X))[rep(1:(n_p + 1), each = n_m), ],
    c(n_p, n_k, n_p + 1, n_n)),
    c(2, 4, 1, 3)) * as.vector(dist.pdf) / sigma, c(3, 1, 4), sum) / n_n,
    nrow = n_m, ncol = n_p + 1)
  Psi_gradient[, n_p + 2] <- colSums(t(matrix(aperm(-aperm(
    dist.x.phi, c(2, 1, 3)) * as.vector(dist.pdf.w) / sigma,
    c(3, 1, 2)), nrow = n_m, ncol = n_n)) * eSI_CS) / n_n
  Psi_gradient[, (n_p + 3):(n_p + n_m + 2)] <- diag(
    rep(-colSums(t(dist.cdf) * eSI_CS) / n_n, each = n_p))
  Psi_gradient[, (n_p + n_m + 3):(n_p * 2 + n_m + 2)] <- matrix(
    rowSums(Psi_i[rep(1:n_m, times = n_p), ] *
              t(X)[rep(1:n_p, each = n_m), ]) / n_n, nrow = n_m, ncol = n_p)

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- matrix(0, n_m, 1)
    Psi_square <- matrix(0, n_m, n_m)
    Psi_gradient <- matrix(0, nrow = n_m, ncol = n_p + n_m + 2)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADCSvar_EXsubY_normal <- function(X, alpha, beta, sigma, phi, CS.beta, y.pts, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  SI <- as.vector(alpha + X %*% beta)
  eSI_CS <- as.vector(exp(X %*% CS.beta))
  a.pts.scaled <- outer_minus_rcpp(y.pts[, 1], SI) / sigma
  b.pts.scaled <- outer_minus_rcpp(y.pts[, 2], SI) / sigma
  dist.cdf <- pnorm(b.pts.scaled, mean = 0, sd = 1) -
    pnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf <- dnorm(b.pts.scaled, mean = 0, sd = 1) -
    dnorm(a.pts.scaled, mean = 0, sd = 1)
  dist.pdf.w <- dnorm(b.pts.scaled, mean = 0, sd = 1) * b.pts.scaled -
    dnorm(a.pts.scaled, mean = 0, sd = 1) * a.pts.scaled
  dist.x.phi <- array(
    X[rep(1:n_n, times = n_k), ] - phi[rep(1:n_k, each = n_n), ],
    c(n_n, n_k, n_p))
  Psi_i <- t(t(matrix(aperm(aperm(
    dist.x.phi, c(2, 1, 3)) * as.vector(dist.cdf), c(3, 1, 2)),
    nrow = n_m, ncol = n_n)) * eSI_CS)
  p_i <- 1 / (1 + as.vector(t(Psi_i) %*% eta)) / n_n

  var.A <- -diag(rep(colSums(t(dist.cdf) * eSI_CS * p_i) / sum(eSI_CS * p_i),
                     each = n_p), nrow = n_m, ncol = n_m)
  var.B.array <- apply(dist.x.phi[, , rep(1:n_p, times = n_p)] *
                         dist.x.phi[, , rep(1:n_p, each = n_p)] *
                         as.vector(t(dist.cdf)) * eSI_CS * p_i, c(3, 2),
                       sum) / sum(eSI_CS * p_i)
  var.B <- matrix(0, n_m, n_m)
  for (k in 1:n_k)
  {
    var.B[(1 + (k - 1) * n_p):(k * n_p),
          (1 + (k - 1) * n_p):(k * n_p)] <- var.B.array[, k]
  }

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

ADCS_EYsubX_normal_Lagrange <- function(X, alpha, beta, sigma,
                                        phi, CS.beta, inclusion, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  eSI_CS <- as.vector(exp(X %*% CS.beta))
  Psi <- outer_minus_rcpp(as.vector(alpha + X %*% beta), phi) * inclusion * eSI_CS
  denominator = as.vector(1 + Psi %*% eta)
  Psi_scaled <- Psi / denominator
  value <- sum(log(denominator))
  gradient <- as.matrix(colSums(Psi_scaled))
  hessian <- matrix(-colSums(Psi_scaled[, rep(1:n_m, times = n_m)] *
                               Psi_scaled[, rep(1:n_m, each = n_m)]),
                    nrow = n_m, ncol = n_m)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_m, 1)
    hessian <- matrix(0, n_m, n_m)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

ADCS_EYsubX_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                             phi, CS.beta, inclusion,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADCS_EYsubX_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS.beta = CS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADCS_EYsubX_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS.beta = CS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EYsubX_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS.beta = CS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EYsubX_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                                phi, CS.beta, inclusion,
                                                eta.initial, iter.max,
                                                step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADCS_EYsubX_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADCS_EYsubX_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EYsubX_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EYsubX_normal <- function(X, alpha, beta, sigma, phi, CS.beta, inclusion)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  eSI_CS <- as.vector(exp(X %*% CS.beta))
  Psi_i <- outer_minus_rcpp(as.vector(alpha + X %*% beta), phi) * inclusion * eSI_CS
  Psi <- as.matrix(colSums(Psi_i) / n_n)
  Psi_square <- matrix(colSums(Psi_i[, rep(1:n_m, times = n_m)] *
                                 Psi_i[, rep(1:n_m, each = n_m)]) / n_n,
                       nrow = n_m, ncol = n_m)
  Psi_gradient <- t(matrix(colSums(
    cbind(1, X, 0)[, rep(1:(n_p + 2), times = n_m)] *
      inclusion[, rep(1:n_m, each = n_p + 2)] * eSI_CS) / n_n,
    nrow = n_p + 2, ncol = n_m))
  Psi_gradient <- cbind(Psi_gradient,
                        diag(-colSums(inclusion * eSI_CS) / n_n, n_m, n_m))
  Psi_gradient <- cbind(Psi_gradient,
                        matrix(colSums(
                          Psi_i[, rep(1:n_m, times = n_p)] *
                            X[, rep(1:n_p, each = n_m)]) / n_n,
                          nrow = n_m, ncol = n_p))

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- matrix(0, n_m, 1)
    Psi_square <- matrix(0, n_m, n_m)
    Psi_gradient <- matrix(0, n_m, n_p + n_m + 2)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADCSvar_EYsubX_normal <- function(X, alpha, beta, sigma, phi, CS.beta, inclusion, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  eSI_CS <- as.vector(exp(X %*% CS.beta))
  Psi_i <- outer_minus_rcpp(as.vector(alpha + X %*% beta), phi) * inclusion * eSI_CS
  p_i <- (1 / (1 + as.vector(Psi_i %*% eta))) / n_n

  var.A <- diag(-colSums(inclusion * eSI_CS * p_i) / sum(eSI_CS * p_i),
                nrow = n_m, ncol = n_m)
  var.B <- colSums(Psi_i[, rep(1:n_m, times = n_m)] *
                     Psi_i[, rep(1:n_m, each = n_m)] / eSI_CS * p_i) /
    sum(eSI_CS * p_i) +
    (sigma ^ 2) * matrix(colSums(
      inclusion[, rep(1:n_m, times = n_m)] *
        inclusion[, rep(1:n_m, each = n_m)] * eSI_CS * p_i) / sum(eSI_CS * p_i),
      nrow = n_m, ncol = n_m)

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

##### Prior probability shift, normal regression model

ADPPS_EX_normal_Lagrange <- function(X, alpha, beta, sigma,
                                     phi, PPS.beta, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  e_f <- as.vector(exp(
    (2 * (alpha + X %*% beta) + (sigma ^ 2) * PPS.beta) * PPS.beta / 2))
  Psi_i <- t(t(X) - as.vector(phi)) * e_f

  denominator <- as.vector(1 + Psi_i %*% eta)
  value <- sum(log(denominator))
  Psi_i_scaled <- Psi_i / denominator
  gradient <- as.matrix(colSums(Psi_i_scaled))
  hessian <- -matrix(colSums(
    Psi_i_scaled[, rep(1:n_p, times = n_p)] *
      Psi_i_scaled[, rep(1:n_p, each = n_p)]),
    nrow = n_p, ncol = n_p)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_p, 1)
    hessian <- matrix(0, n_p, n_p)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

ADPPS_EX_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                          phi, PPS.beta,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EX_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS.beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EX_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS.beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EX_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS.beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EX_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                             phi, PPS.beta,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EX_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS_beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EX_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS_beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EX_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS_beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EX_normal <- function(X, alpha, beta, sigma, phi, PPS.beta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  Psi_gradient <- matrix(0, nrow = n_p, ncol = n_p * 2 + 3)

  e_f <- as.vector(exp(
    (2 * (alpha + X %*% beta) + (sigma ^ 2) * PPS.beta) * PPS.beta / 2))
  Psi_i <- t(t(X) - as.vector(phi)) * e_f

  Psi <- as.matrix(colSums(Psi_i) / n_n)
  Psi_square <- matrix(colSums(
    Psi_i[, rep(1:n_p, times = n_p)] *
      Psi_i[, rep(1:n_p, each = n_p)]) / n_n, nrow = n_p, ncol = n_p)
  Psi_gradient[, 1:(n_p + 1)] <- matrix(colSums(
    Psi_i[, rep(1:n_p, times = n_p + 1)] *
      cbind(1, X)[, rep(1:(n_p + 1), each = n_p)]) / n_n,
    nrow = n_p, ncol = n_p + 1)
  Psi_gradient[, n_p + 2] <- Psi * sigma * (PPS.beta ^ 2)
  Psi_gradient[, (n_p + 3):(n_p * 2 + 2)] <- -diag(sum(e_f) / n_n,
                                                   nrow = n_p, ncol = n_p)
  Psi_gradient[, n_p * 2 + 3] <- Psi * (sigma ^ 2) * PPS.beta

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- matrix(0, n_p, 1)
    Psi_square <- matrix(0, n_p, n_p)
    Psi_gradient <- matrix(0, nrow = n_p, ncol = n_p * 2 + 3)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADPPSvar_EX_normal <- function(X, alpha, beta, sigma, phi, PPS.beta, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  e_f <- as.vector(exp(
    (2 * (alpha + X %*% beta) + (sigma ^ 2) * PPS.beta) * PPS.beta / 2))
  dist.x.phi <- t(t(X) - as.vector(phi))
  Psi_i <- dist.x.phi * e_f
  p_i <- 1 / (1 + as.vector(Psi_i %*% eta)) / n_n

  var.A <- -diag(1, nrow = n_p, ncol = n_p)
  var.B <- matrix(colSums(dist.x.phi[, rep(1:n_p, times = n_p)] *
                            dist.x.phi[, rep(1:n_p, each = n_p)] * e_f * p_i) /
                    sum(e_f * p_i), nrow = n_p, ncol = n_p)

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

ADPPS_EY_normal_Lagrange <- function(X, alpha, beta, sigma,
                                     phi, PPS.beta, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  Psi_i <- (SI + (sigma ^ 2) * PPS.beta - phi) * e_f

  denominator <- 1 + Psi_i * eta
  value <- sum(log(denominator))
  Psi_i_scaled <- Psi_i / denominator
  gradient <- sum(Psi_i_scaled)
  hessian <- -sum(Psi_i_scaled ^ 2)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- 0
    hessian <- 0
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

ADPPS_EY_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                          phi, PPS.beta,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EY_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS.beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EY_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS.beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EY_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS.beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EY_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                             phi, PPS.beta,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EY_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS_beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EY_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS_beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EY_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS_beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EY_normal <- function(X, alpha, beta, sigma, phi, PPS.beta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  c_f <- SI + (sigma ^ 2) * PPS.beta - phi
  Psi_i <- c_f * e_f

  Psi <- sum(Psi_i) / n_n
  Psi_square <- sum(Psi_i ^ 2) / n_n
  Psi_gradient <- c(
    colSums(cbind(1, X) * e_f * (1 + c_f * PPS.beta)) / n_n,
    sum(e_f * (2 + c_f * PPS.beta) * sigma * PPS.beta) / n_n,
    -sum(e_f) / n_n,
    sum(e_f * (1 + c_f * PPS.beta) * (sigma ^ 2)) / n_n
  )

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- 0
    Psi_square <- 0
    Psi_gradient <- rep(0, n_p + 4)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADPPSvar_EY_normal <- function(X, alpha, beta, sigma, phi, PPS.beta, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  c_f <- SI + (sigma ^ 2) * PPS.beta - phi
  Psi_i <- c_f * e_f
  p_i <- 1 / (1 + Psi_i * eta) / n_n

  var.A <- -1
  var.B <- sum(c_f ^ 2 * e_f * p_i) / sum(e_f * p_i) + sigma ^ 2

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

ADPPS_EXsubY_normal_Lagrange <- function(X, alpha, beta, sigma,
                                         phi, PPS.beta, y.pts, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  mu_shifted <- SI + (sigma ^ 2) * PPS.beta
  dist_x_phi <- array(X[rep(1:n_n, times = n_k), ] -
    phi[rep(1:n_k, each = n_n), ], c(n_n, n_k, n_p))
  cdf_dist <- t(pnorm(outer_minus_rcpp(y.pts[, 2], mu_shifted) / sigma, 0, 1) -
    pnorm(outer_minus_rcpp(y.pts[, 1], mu_shifted) / sigma, 0, 1))
  Psi_i <- matrix(aperm(dist_x_phi * e_f * as.vector(cdf_dist), c(1, 3, 2)),
                  nrow = n_n, ncol = n_m)

  Psi <- as.matrix(colSums(Psi_i))
  denominator <- as.vector(1 + Psi_i %*% eta)
  value <- sum(log(denominator))
  Psi_i_scaled <- Psi_i / denominator
  gradient <- as.matrix(colSums(Psi_i_scaled))
  hessian <- -matrix(colSums(Psi_i_scaled[, rep(1:n_m, times = n_m)] *
                               Psi_i_scaled[, rep(1:n_m, each = n_m)]),
                     nrow = n_m, ncol = n_m)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_m, 1)
    hessian <- matrix(0, n_m, n_m)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

ADPPS_EXsubY_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                              phi, PPS.beta, y.pts,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EXsubY_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS.beta = PPS.beta, y.pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EXsubY_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS.beta = PPS.beta, y.pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EXsubY_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS.beta = PPS.beta, y.pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EXsubY_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                                 phi, PPS.beta, y.pts,
                                                 eta.initial, iter.max,
                                                 step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EXsubY_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS_beta = PPS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EXsubY_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS_beta = PPS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EXsubY_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS_beta = PPS.beta, y_pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EXsubY_normal <- function(X, alpha, beta, sigma,
                                phi, PPS.beta, y.pts)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  Psi_gradient <- matrix(0, nrow = n_m, ncol = n_p + n_m + 3)

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  mu_shifted <- SI + (sigma ^ 2) * PPS.beta
  dist_x_phi <- array(X[rep(1:n_n, times = n_k), ] -
                        phi[rep(1:n_k, each = n_n), ], c(n_n, n_k, n_p))
  y_pts_a <- t(outer_minus_rcpp(y.pts[, 1], mu_shifted) / sigma)
  y_pts_b <- t(outer_minus_rcpp(y.pts[, 2], mu_shifted) / sigma)
  cdf_dist <- pnorm(y_pts_b, 0, 1) - pnorm(y_pts_a, 0, 1)
  pdf_a <- dnorm(y_pts_a, 0, 1)
  pdf_b <- dnorm(y_pts_b, 0, 1)
  pdf_dist <- pdf_b - pdf_a
  Psi_i <- matrix(aperm(dist_x_phi * e_f * as.vector(cdf_dist), c(1, 3, 2)),
                  nrow = n_n, ncol = n_m)

  Psi <- as.matrix(colSums(Psi_i) / n_n)
  Psi_square <- matrix(colSums(
    Psi_i[, rep(1:n_m, times = n_m)] *
      Psi_i[, rep(1:n_m, each = n_m)]) / n_n, nrow = n_m, ncol = n_m)
  Psi_gradient[, 1:(n_p + 1)] <- matrix(colSums(matrix(aperm(
    dist_x_phi * e_f * as.vector(cdf_dist * PPS.beta - pdf_dist / sigma),
    c(1, 3, 2)), nrow = n_n, ncol = n_m)[, rep(1:n_m, times = n_p + 1)] *
      cbind(1, X)[, rep(1:(n_p + 1), each = n_m)]) / n_n,
    nrow = n_m, ncol = n_p + 1)
  Psi_gradient[, n_p + 2] <- as.vector(apply(
    dist_x_phi * e_f * as.vector(t(
      cdf_dist * sigma * (PPS.beta ^ 2) -
        (pdf_b * (y_pts_b / sigma + PPS.beta * 2) -
           pdf_a * (y_pts_a / sigma + PPS.beta * 2)))),
    c(3, 2), sum) / n_n)
  Psi_gradient[, (n_p + 3):(n_p + n_m + 2)] <- diag(rep(
    -colSums(cdf_dist * e_f) / n_n, each = n_p), nrow = n_m, ncol = n_m)
  Psi_gradient[, n_p + n_m + 3] <- as.vector(apply(
    dist_x_phi * e_f * as.vector(
      cdf_dist * (sigma ^ 2) * PPS.beta - pdf_dist * sigma),
    c(3, 2), sum) / n_n)

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- matrix(0, n_m, 1)
    Psi_square <- matrix(0, n_m, n_m)
    Psi_gradient <- matrix(0, nrow = n_m, ncol = n_p + n_m + 3)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADPPSvar_EXsubY_normal <- function(X, alpha, beta, sigma,
                                   phi, PPS.beta, y.pts, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_k <- dim(phi)[1]
  n_m <- n_k * n_p

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  mu_shifted <- SI + (sigma ^ 2) * PPS.beta
  dist_x_phi <- array(X[rep(1:n_n, times = n_k), ] -
                        phi[rep(1:n_k, each = n_n), ], c(n_n, n_k, n_p))
  y_pts_a <- t(outer_minus_rcpp(y.pts[, 1], mu_shifted) / sigma)
  y_pts_b <- t(outer_minus_rcpp(y.pts[, 2], mu_shifted) / sigma)
  cdf_dist <- pnorm(y_pts_b, 0, 1) - pnorm(y_pts_a, 0, 1)
  pdf_a <- dnorm(y_pts_a, 0, 1)
  pdf_b <- dnorm(y_pts_b, 0, 1)
  pdf_dist <- pdf_b - pdf_a
  Psi_i <- matrix(aperm(dist_x_phi * e_f * as.vector(cdf_dist), c(1, 3, 2)),
                  nrow = n_n, ncol = n_m)
  p_i <- 1 / (1 + as.vector(Psi_i %*% eta)) / n_n

  var.A <- diag(rep(-colSums(cdf_dist * e_f * p_i) /
                      sum(e_f * p_i), each = n_p),
                nrow = n_m, ncol = n_m)
  var.B.array <- apply(dist_x_phi[, , rep(1:n_p, times = n_p)] *
    dist_x_phi[, , rep(1:n_p, each = n_p)] * as.vector(cdf_dist) *
    e_f * p_i, c(2, 3), sum) / sum(e_f * p_i)
  var.B <- matrix(0, n_m, n_m)
  for (k in 1:n_k)
  {
    var.B[(1 + (k - 1) * n_p):(k * n_p),
          (1 + (k - 1) * n_p):(k * n_p)] <- var.B.array[k, ]
  }

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

ADPPS_EYsubX_normal_Lagrange <- function(X, alpha, beta, sigma,
                                         phi, PPS.beta, inclusion, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  Psi_i <- outer_minus_rcpp(SI + (sigma ^ 2) * PPS.beta, phi) * inclusion * e_f

  denominator <- as.vector(1 + Psi_i %*% eta)
  value <- sum(log(denominator))
  Psi_i_scaled <- Psi_i / denominator
  gradient <- as.matrix(colSums(Psi_i_scaled))
  hessian <- matrix(-colSums(Psi_i_scaled[, rep(1:n_m, times = n_m)] *
                               Psi_i_scaled[, rep(1:n_m, each = n_m)]),
                    nrow = n_m, ncol = n_m)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_m, 1)
    hessian <- matrix(0, n_m, n_m)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

ADPPS_EYsubX_normal_SolveLagrange <- function(X, alpha, beta, sigma,
                                              phi, PPS.beta, inclusion,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EYsubX_normal_Lagrange(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS.beta = PPS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EYsubX_normal_Lagrange(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS.beta = PPS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EYsubX_normal_Lagrange(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS.beta = PPS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EYsubX_normal_SolveLagrange_v1 <- function(X, alpha, beta, sigma,
                                                 phi, PPS.beta, inclusion,
                                                 eta.initial, iter.max,
                                                 step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- ADPPS_EYsubX_normal_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    direction.step <- solve_rcpp(
      step$hessian - eps_inv * diag(dim(step$hessian)[1]),
      step$gradient)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- ADPPS_EYsubX_normal_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EYsubX_normal_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EYsubX_normal <- function(X, alpha, beta, sigma, phi, PPS.beta, inclusion)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  c_f <- outer_minus_rcpp(SI + (sigma ^ 2) * PPS.beta, phi)
  Psi_i <- c_f * e_f * inclusion

  Psi <- as.matrix(colSums(Psi_i) / n_n)
  Psi_square <- matrix(colSums(
    Psi_i[, rep(1:n_m, times = n_m)] *
      Psi_i[, rep(1:n_m, each = n_m)]) / n_n, nrow = n_m, ncol = n_m)
  Psi_gradient <- cbind(
    matrix(colSums(
      (inclusion * (1 + c_f * PPS.beta))[, rep(1:n_m, times = n_p + 1)] *
        cbind(1, X)[, rep(1:(n_p + 1), each = n_m)] * e_f) / n_n,
      nrow = n_m, ncol = n_p + 1),
    colSums(inclusion * (2 + c_f * PPS.beta) * e_f * sigma * PPS.beta) / n_n,
    diag(-colSums(inclusion * e_f) / n_n, nrow = n_m, ncol = n_m),
    colSums(inclusion * (1 + c_f * PPS.beta) * e_f * (sigma ^ 2)) / n_n
  )

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- 0
    Psi_square <- 0
    Psi_gradient <- rep(0, n_p + 4)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADPPSvar_EYsubX_normal <- function(X, alpha, beta, sigma,
                                   phi, PPS.beta, inclusion, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- length(phi)

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * PPS.beta) * PPS.beta / 2)
  c_f <- outer_minus_rcpp(SI + (sigma ^ 2) * PPS.beta, phi)
  Psi_i <- c_f * e_f * inclusion
  p_i <- 1 / (1 + as.vector(Psi_i %*% eta)) / n_n

  var.A <- diag(-colSums(inclusion * e_f * p_i) / sum(e_f * p_i),
                nrow = n_m, ncol = n_m)
  var.B <- matrix(colSums(
    Psi_i[, rep(1:n_m, times = n_m)] *
      Psi_i[, rep(1:n_m, each = n_m)] / e_f * p_i) / sum(e_f * p_i),
    nrow = n_m, ncol = n_m) + (sigma ^ 2) * matrix(colSums(
      inclusion[, rep(1:n_m, times = n_m)] *
        inclusion[, rep(1:n_m, each = n_m)] * e_f * p_i) / sum(e_f * p_i),
      nrow = n_m, ncol = n_m)

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

##### No shift, logistic regression model

AD_EY_logistic_Lagrange <- function(X, alpha, beta, phi, eta)
{
  n_n <- dim(X)[1]

  eSI <-exp(alpha + X %*% beta)
  Psi <- ((1 - phi) * eSI - phi) / (1 + eSI)
  denominator = 1 + eta * Psi
  value <- sum(log(denominator))
  gradient <- sum(Psi / denominator)
  hessian <- -sum((Psi / denominator) ^ 2)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- 0
    hessian <- 0
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

AD_EY_logistic_SolveLagrange <- function(X, alpha, beta, phi,
                                         eta.initial, iter.max,
                                         step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EY_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, eta = eta
  )

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EY_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EY_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EY_logistic_SolveLagrange_v1 <- function(X, alpha, beta, phi,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol, eps_inv)
{
  eta <- eta.initial
  step <- AD_EY_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, eta = eta
  )

  for (k in 1:iter.max)
  {
    direction.step <- step$gradient / (step$hessian - eps_inv)

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EY_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EY_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EY_logistic <- function(X, alpha, beta, phi)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  eSI_i <- as.vector(exp(alpha + X %*% beta))
  Psi_i <- ((1 - phi) * eSI_i - phi) / (1 + eSI_i)
  Psi <- sum(Psi_i) / n_n
  Psi_square <- sum(Psi_i ^ 2) / n_n

  Psi_gradient <- colSums(
    cbind(cbind(1, X) * eSI_i / (1 + eSI_i), -1) / (1 + eSI_i)
    ) / n_n

  if (any(!is.finite(Psi))|
      any(!is.finite(Psi_square))|
      any(!is.finite(Psi_gradient)))
  {
    Psi <- 0
    Psi_square <- 0
    Psi_gradient <- rep(0, n_p + 3)
  }

  results <- list(score = Psi,
                  score_square = Psi_square,
                  score_gradient = Psi_gradient)
  return(results)
}

ADvar_EY_logistic <- function(X, alpha, beta, phi, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]

  eSI_i <- as.vector(exp(alpha + X %*% beta))
  Psi_i <- ((1 - phi) * eSI_i - phi) / (1 + eSI_i)
  p_i <- (1 / (1 + Psi_i * eta)) / n_n

  var.A <- -1
  var.B <- sum(((1 - phi) ^ 2 * eSI_i + phi ^ 2) / (1 + eSI_i) * p_i)

  results <- list(var.A = var.A,
                  var.B = var.B)
  return(results)
}

AD_EXsubY_logistic_Lagrange <- function(X, alpha, beta, phi, eta)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- 2 * n_p

  eSI_i <- as.vector(exp(alpha + X %*% beta))
  Psi_i <- cbind(t(t(X) - phi[1, ]) * eSI_i / (1 + eSI_i),
                 t(t(X) - phi[2, ]) / (1 + eSI_i))
  denominator <- as.vector(1 + Psi_i %*% eta)
  value <- sum(log(denominator))
  Psi_i_scaled <- Psi_i / denominator
  gradient <- as.matrix(colSums(Psi_i_scaled))
  hessian <- -matrix(colSums(Psi_i_scaled[, rep(1:n_m, times = n_m)] *
                               Psi_i_scaled[, rep(1:n_m, each = n_m)]),
                     nrow = n_m, ncol = n_m)

  if (any(!is.finite(value))|
      any(!is.finite(gradient))|
      any(!is.finite(hessian)))
  {
    value <- -10000 * n_n
    gradient <- matrix(0, n_m, 1)
    hessian <- matrix(0, n_m, n_m)
  }

  results <- list(value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}












##### No shift, Gamma regression model

AD_EY_Gamma_SolveLagrange <- function(X, alpha, beta, nu, phi,
                                      eta.initial, iter.max,
                                      step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- AD_EY_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, eta = eta
  )

  for (k in 1:iter.max)
  {
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EY_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EY_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EXsubY_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                          phi, y.pts,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- AD_EXsubY_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EXsubY_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EXsubY_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, y_pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EYsubX_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                          phi, inclusion,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- AD_EYsubX_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- AD_EYsubX_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EYsubX_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

##### Prior probability shift, Gamma regression model

ADPPS_EX_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                         phi, PPS.beta,
                                         eta.initial, iter.max,
                                         step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADPPS_EX_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, PPS_beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADPPS_EX_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, PPS_beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EX_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, PPS_beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EY_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                         phi, PPS.beta,
                                         eta.initial, iter.max,
                                         step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADPPS_EY_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, PPS_beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADPPS_EY_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, PPS_beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EY_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, PPS_beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EXsubY_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                             phi, PPS.beta, y.pts,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADPPS_EXsubY_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, PPS_beta = PPS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADPPS_EXsubY_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, PPS_beta = PPS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EXsubY_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, PPS_beta = PPS.beta, y_pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EYsubX_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                             phi, PPS.beta, inclusion,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADPPS_EYsubX_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADPPS_EYsubX_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EYsubX_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

##### Covariate shift, Gamma regression model

ADCS_EY_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                        phi, CS.beta,
                                        eta.initial, iter.max,
                                        step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADCS_EY_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, CS_beta = CS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADCS_EY_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, CS_beta = CS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EY_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, CS_beta = CS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EXsubY_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                            phi, CS.beta, y.pts,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADCS_EXsubY_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADCS_EXsubY_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EXsubY_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EYsubX_Gamma_SolveLagrange <- function(X, alpha, beta, nu,
                                            phi, CS.beta, inclusion,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADCS_EYsubX_Gamma_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADCS_EYsubX_Gamma_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EYsubX_Gamma_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}













AD_EXsubY_logistic_SolveLagrange <- function(X, alpha, beta, phi,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- AD_EXsubY_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, eta = eta)

  for (k in 1:iter.max)
  {
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }

    step.size <- 1
    eta.new <- eta - direction.step
    step.new <- AD_EXsubY_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EXsubY_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

AD_EYsubX_logistic_SolveLagrange <- function(X, alpha, beta,
                                             phi, inclusion,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- AD_EYsubX_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- AD_EYsubX_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- AD_EYsubX_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

##### Prior probability shift, logistic regression model

ADPPS_EX_logistic_SolveLagrange <- function(X, alpha, beta,
                                            phi, PPS.beta,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADPPS_EX_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, PPS_beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADPPS_EX_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, PPS_beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EX_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, PPS_beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EY_logistic_SolveLagrange <- function(X, alpha, beta,
                                            phi, PPS.beta,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADPPS_EY_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, PPS_beta = PPS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADPPS_EY_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, PPS_beta = PPS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EY_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, PPS_beta = PPS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADPPS_EYsubX_logistic_SolveLagrange <- function(X, alpha, beta,
                                                phi, PPS.beta, inclusion,
                                                eta.initial, iter.max,
                                                step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADPPS_EYsubX_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADPPS_EYsubX_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADPPS_EYsubX_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, PPS_beta = PPS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

##### Covariate shift, logistic regression model

ADCS_EY_logistic_SolveLagrange <- function(X, alpha, beta,
                                           phi, CS.beta,
                                           eta.initial, iter.max,
                                           step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADCS_EY_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, CS_beta = CS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADCS_EY_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, CS_beta = CS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EY_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, CS_beta = CS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EXsubY_logistic_SolveLagrange <- function(X, alpha, beta,
                                               phi, CS.beta,
                                               eta.initial, iter.max,
                                               step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADCS_EXsubY_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, CS_beta = CS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADCS_EXsubY_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, CS_beta = CS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EXsubY_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, CS_beta = CS.beta, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

ADCS_EYsubX_logistic_SolveLagrange <- function(X, alpha, beta,
                                               phi, CS.beta, inclusion,
                                               eta.initial, iter.max,
                                               step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- ADCS_EYsubX_logistic_Lagrange_rcpp(
    X = X, alpha = alpha, beta = beta,
    phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- ADCS_EYsubX_logistic_Lagrange_rcpp(
      X = X, alpha = alpha, beta = beta,
      phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value + tol)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- ADCS_EYsubX_logistic_Lagrange_rcpp(
          X = X, alpha = alpha, beta = beta,
          phi = phi, CS_beta = CS.beta, inclusion = inclusion, eta = eta.new)
      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new
    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}
