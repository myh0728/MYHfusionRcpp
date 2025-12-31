n <- 200
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

ord <- order(test.data.1$time.last)
time_last <- test.data.1$time.last[ord]
is_event <- test.data.1$is.event[ord]
time_event <- sort(unique(
  test.data.1$time.last[test.data.1$is.event == 1]))

test1 <- KME_R(time.last = time_last,
               is.event = is_event,
               time.event = time_event)
test2 <- KME_exact_rcpp(time_last = time_last,
                        is_event = is_event,
                        time_event = time_event)

all.equal(as.vector(test1$survival), as.vector(test2$survival))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = KME_R(time.last = time_last,
                is.event = is_event,
                time.event = time_event),
    "Rcpp" = KME_exact_rcpp(time_last = time_last,
                            is_event = is_event,
                            time_event = time_event)
  )
)

test <- KME_at_times_rcpp(time_event = time_event,
                          S_event = test2$survival,
                          H_event = test2$cum_hazard,
                          time_points = seq(0, 5, 0.1))

test3 <- SKME_R(time.last = time_last,
                is.event = is_event,
                time.event = time_event,
                X = Xi, x = Xi, K = K2_Ep_R, h = 1.5)

test4 <- SKME_exact_rcpp(time_last = time_last,
                         is_event = is_event,
                         time_event = time_event,
                         X = Xi, x = Xi, h = rep(1.5, length = p))

all.equal(as.vector(test3$survival), as.vector(test4$survival))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = SKME_R(time.last = time_last,
                 is.event = is_event,
                 time.event = time_event,
                 X = Xi, x = Xi, K = K2_Ep_R, h = 1.5),
    "Rcpp" = SKME_exact_rcpp(time_last = time_last,
                             is_event = is_event,
                             time_event = time_event,
                             X = Xi, x = Xi, h = rep(1.5, length = p))
  )
)












