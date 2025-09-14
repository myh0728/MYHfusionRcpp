n <- 100
p <- 5

Xi <- matrix(rnorm(n * p), n, p)
Ci <- 1 + rexp(n)

test.data.1 <- simSurv.AFT.normal(Xi = Xi, Ci = Ci,
                                  alpha0 = 1, beta0 = rep(1, p),
                                  sigma.error = 1)

test.data.2 <- simSurv.PH(Xi = Xi, Ci = Ci,
                          beta0 = rep(1, p),
                          cumhazard.inv = function(t){t})

w <- rexp(n)
w <- w / sum(w)

###

test1 <- KME_R_outer(time.last = test.data.1$time.last,
                     is.event = test.data.1$is.event,
                     time.event = sort(unique(
                       test.data.1$time.last[test.data.1$is.event == 1])))
test2 <- KME_rcpp(time_last = test.data.1$time.last,
                  is_event = test.data.1$is.event,
                  time_event = sort(unique(
                    test.data.1$time.last[test.data.1$is.event == 1])))
sum(abs(test1$hazard - test2))

time.jumps <- sort(unique(c(seq(0, 3, 0.1), test.data.1$time.last[test.data.1$is.event == 1])))

test1 <- KME_R_outer(time.last = test.data.1$time.last,
                     is.event = test.data.1$is.event,
                     time.event = time.jumps)
test2 <- KME_rcpp(time_last = test.data.1$time.last,
                  is_event = test.data.1$is.event,
                  time_event = time.jumps)
sum(abs(test1$hazard - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = KME_R_outer(time.last = test.data.1$time.last,
                      is.event = test.data.1$is.event,
                      time.event = sort(unique(
                        test.data.1$time.last[test.data.1$is.event == 1]))),
    "Rcpp" = KME_rcpp(time_last = test.data.1$time.last,
                      is_event = test.data.1$is.event,
                      time_event = sort(unique(
                        test.data.1$time.last[test.data.1$is.event == 1])))
  )
)

###

test1 <- SKME_R_outer(time.last = test.data.1$time.last,
                      is.event = test.data.1$is.event,
                      time.event = sort(unique(
                        test.data.1$time.last[test.data.1$is.event == 1])),
                      X = Xi, x = Xi, K = K2_Bw, h = 1.5)


