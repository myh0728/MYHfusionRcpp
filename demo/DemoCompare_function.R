square_fun <- function(a){

  return(a ^ 2)
}

testfunction_R <- function(a, b){

  return(b(a))
}

testfunction_R(5, square_fun)
testfunction_rcpp(5, square_fun)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    "R" = testfunction_R(5, square_fun),
    "Rcpp" = testfunction_rcpp(5, square_fun),
    "direct" = 5 ^ 2
  )
)




