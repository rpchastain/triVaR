#' Inverse CDF of the truncated exponential distribution intended to be used 
#' for sampling
#'
#' @param x Input value, expected to be between 0 and 1
#' @param rate Single value specifying the rate.
#' @param a Lower truncation bound, if not provided it assumes the default value
#' of 0.
#' @param b Upper truncation bound.
#'
#' @return Numeric vector
#' @export 
#'
#' @examples
#' inv_texp(0.25, b = 6)
#' inv_texp(0.5, rate = -1, a = 2, b = 6)
inv_texp <- function(x, rate = 1, a = 0, b) {
  log(exp(-rate * a) - x * (exp(-rate * a ) - exp(-rate * b)))/(-rate)
}


#' Truncated exponential random number generator
#'
#' @param n Number of random samples to draw
#' @param rate Single value specifying the rate.
#' @param a Lower truncation bound, if not provided it assumes the default value
#' of 0.
#' @param b Upper truncation bound.
#'
#' @return Numeric vector
#' @export
#'
#' @examples
#' rtexp(5, b = 6)
#' rtexp(2, rate = -1, a = 2, b = 6)
rtexp <- function(n, rate = 1, a = 0, b) {
  u <- stats::runif(n)
  inv_texp(u, rate, a, b)
}
