######################
### linear algebra ###
######################

A <- matrix(rnorm(100 ^ 2), 100, 100)

test1 <- pracma::pinv(A)
test2 <- pinv_rcpp(A)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "pracma" = pracma::pinv(A),
    "Rcpp" = pinv_rcpp(A),
    "Rcpp_v1" = pinv_rcpp_v1(A),
    "Rcpp_v2" = pinv_rcpp_v2(A),
    "Rcpp_v3" = pinv_rcpp_v3(A),
    "Rcpp_v4" = pinv_rcpp_v4(A)
  )
)

###

A <- matrix(rnorm(50 ^ 2), 50, 50)
B <- matrix(rnorm(50 ^ 2), 50, 10)

test1 <- solve(A, B)
test2 <- solve_rcpp(A, B)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "solve" = solve(A, B),
    "Rcpp" = solve_rcpp(A, B)
  )
)

###

A <- matrix(rnorm(50 ^ 2), 50, 50)
B <- t(A) %*% A

test1 <- solve(B)
test2 <- inv_sympd_rcpp(B)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "solve" = solve(B),
    "Rcpp" = inv_sympd_rcpp(B)
  )
)

###

A <- matrix(rnorm(50 ^ 2), 50, 50)
B <- t(A) %*% A

test1 <- eigen(B)
sum(abs(test1$vectors %*% diag(test1$values) %*% t(test1$vectors) - B))

test2 <- eigen_rcpp(B)
sum(abs(test2$vector %*% diag(as.vector(test2$value)) %*% t(test2$vector) - B))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "eigen" = eigen(B),
    "Rcpp" = eigen_rcpp(B),
    "Rcpp_v1" = eigen_rcpp_v1(B)
  )
)

###

A <- matrix(rnorm(50 ^ 2), 50, 50)
B <- t(A) %*% A

test1 <- chol(B)
test2 <- chol_rcpp(B)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "chol" = chol(B),
    "Rcpp" = chol_rcpp(B)
  )
)

###

A <- matrix(rnorm(50 ^ 2), 50, 50)

test1 <- rcond(A)
test2 <- rcond_rcpp(A)
abs(test1 - test2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "rcond" = rcond(A),
    "Rcpp" = rcond_rcpp(A)
  )
)

#########################
### statistal algebra ###
#########################

n <- 1000
p <- 10
a <- sample(1:100, size = n, replace = TRUE)
A <- matrix(rnorm(n * p), n, p)
Aa <- data.frame(A = A, a = a)

test1 <- aggregate(A, by = list(a), FUN = sum)[, -1]
test2 <- GroupSum_rcpp(A, a)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "aggregate" = aggregate(A, by = list(a), FUN = sum)[, -1],
    "Rcpp" = GroupSum_rcpp(A, a)
  )
)

###

n <- 1000
p <- 10
w <- rexp(n)

A <- matrix(rnorm(n * p), n, p)

test1 <- eXsq_rcpp(A)
test2 <- eXsq_R(A)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = eXsq_R(A),
    "Rcpp" = eXsq_rcpp(A)
  )
)

test1 <- eXsq_w_rcpp(A, w)
test2 <- eXsq_w_R(A, w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = eXsq_w_R(A, w),
    "Rcpp" = eXsq_w_rcpp(A, w)
  )
)

###

n <- 1000
p <- 10

A <- matrix(rnorm(n * p), n, p)

test1 <- Xsq_lowtri_rcpp(A)
test2 <- Xsq_lowtri_R(A)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = Xsq_lowtri_R(A),
    "Rcpp" = Xsq_lowtri_rcpp(A)
  )
)

###

n <- 1000
p <- 10

A <- matrix(rnorm(n * p), n, p)
B <- matrix(rnorm(n * p), n, p)

test1 <- twoXYsym_lowtri_rcpp(A, B)
test2 <- twoXYsym_lowtri_R(A, B)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = twoXYsym_lowtri_R(A, B),
    "Rcpp" = twoXYsym_lowtri_rcpp(A, B)
  )
)

###

n <- 1000
p <- 1
k <- 50

Y <- matrix(rnorm(n * p), n, p)
y <- matrix(rnorm(k * p), k, p)

test1 <- ctingP_R(Y, y)
test2 <- ctingP_uni_R(as.vector(Y), as.vector(y))
test3 <- ctingP_uni_R_outer(as.vector(Y), as.vector(y))
test4 <- ctingP_rcpp(Y, y)
test5 <- ctingP_uni_rcpp(as.vector(Y), as.vector(y))
sum(abs(test1 - test2))
sum(abs(test1 - test3))
sum(abs(test1 - test4))
sum(abs(test1 - test5))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = ctingP_R(Y, y),
    "Rcpp" = ctingP_rcpp(Y, y),
    "R_uni" = ctingP_uni_R(as.vector(Y), as.vector(y)),
    "R_uni_outer" = ctingP_uni_R_outer(as.vector(Y), as.vector(y)),
    "Rcpp_uni" = ctingP_uni_rcpp(as.vector(Y), as.vector(y))
  )
)

######################
### count and rank ###
######################

n <- 1000
a <- sample(1:100, size = n, replace = TRUE)
w <- rexp(n)

test1 <- countAinB_outer(1:100, a)
test2 <- countAinB_rcpp(1:100, a)
test3 <- countAinB_for(1:100, a)
sum(abs(test1 - test2))
sum(abs(test1 - test3))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "outer" = countAinB_outer(1:100, a),
    "Rcpp" = countAinB_rcpp(1:100, a),
    "for" = countAinB_for(1:100, a)
  )
)

test1 <- countAinB_W_outer(1:100, a, w)
test2 <- countAinB_W_rcpp(1:100, a, w)
test3 <- countAinB_W_for(1:100, a, w)
sum(abs(test1 - test2))
sum(abs(test1 - test3))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "outer" = countAinB_W_outer(1:100, a, w),
    "Rcpp" = countAinB_W_rcpp(1:100, a, w),
    "for" = countAinB_W_for(1:100, a, w)
  )
)

###

a <- runif(1000)
b <- seq(0, 1, 0.1)

test1 <- rankAinB_for(a, b)
test2 <- rankAinB_outer(a, b)
test3 <- rankAinB_rcpp(a, b)
sum(abs(test1 - test2))
sum(abs(test1 - test3))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "outer" = rankAinB_outer(a, b),
    "Rcpp" = rankAinB_rcpp(a, b),
    "for" = rankAinB_for(a, b)
  )
)

#############
### outer ###
#############

a <- 1:1000
b <- 1000:1

test1 <- outer(a, b, FUN = "*")
test2 <- outer_times_rcpp(a, b)
sum(abs(test1 - test2))

test1 <- outer(a, b, FUN = "-")
test2 <- outer_minus_rcpp(a, b)
sum(abs(test1 - test2))

test1 <- outer(a, b, FUN = "+")
test2 <- outer_plus_rcpp(a, b)
sum(abs(test1 - test2))

test1 <- outer(a, b, FUN = "<=")
test2 <- outer_leq_rcpp(a, b)
sum(abs(test1 - test2))

test1 <- outer(a, b, FUN = ">=")
test2 <- outer_geq_rcpp(a, b)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "times_outer" = outer(a, b, FUN = "*"),
    "times_Rcpp" = outer_times_rcpp(a, b),
    "minus_outer" = outer(a, b, FUN = "-"),
    "minus_Rcpp" = outer_minus_rcpp(a, b),
    "plus_outer" = outer(a, b, FUN = "+"),
    "plus_Rcpp" = outer_plus_rcpp(a, b),
    "leq_outer" = outer(a, b, FUN = "<="),
    "leq_Rcpp" = outer_leq_rcpp(a, b),
    "geq_outer" = outer(a, b, FUN = ">="),
    "geq_Rcpp" = outer_geq_rcpp(a, b)
  )
)





