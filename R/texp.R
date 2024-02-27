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


#' MSTE-N aggregation function intended for package use only.
#'
#' @param ... Parameters required for the `rtexp` function.
#'
#' @return Dataframe consisting of `sum` and `max` of `n` `rtexp` random 
#' variables
#' @export
#'
#' @examples get_msten(n = 4, rate = -1, b = 6)
get_msten <- function(...) {
  rand_sample <- triVaR::rtexp(...)
  data.frame(x = sum(rand_sample), y = max(rand_sample))
}

get_tetlg <- function(...) {
  rand_sample <- stats::rexp(...)
  data.frame(x = sum(rand_sample), y = max(rand_sample))
}

#' MSTE-N random vector generator
#'
#' @param n_dist The counting distribution which N from the MSTE-N vector will
#' be drawn. The number of random variables drawn from this distribution 
#' determines the number of MSTE-N random vectors drawn.
#' @param rate Single value specifying the rate.
#' @param a Lower truncation bound, if not provided it assumes the default value
#' of 0.
#' @param b Upper truncation bound.
#' @param ... The named parameters required for distribution supplied to 
#' `n_dist`.
#'
#' @return Numeric vector
#' @export
#'
#' @examples
#' 
#' rmsten(n_dist = rgeom, rate = -1, a = 2, b = 6, p = 0.25, n = 6)
rmsten <- function(n_dist, rate = 1, a = 0, b, ...) {
  event_len <- n_dist(...)
  rand_sample <- do.call(
    rbind.data.frame,
    lapply(event_len, get_msten, rate, a, b)
    )
  rand_sample[['N']] <- event_len
  rand_sample
}

#' TETLG random vector generator
#'
#' @param n_dist The counting distribution which N from the TETLG vector will
#' be drawn. The number of random variables drawn from this distribution 
#' determines the number of TETLG random vectors drawn.
#' @param rate Single value specifying the rate.
#' @param ... The named parameters required for distribution supplied to 
#' `n_dist`.
#'
#' @return Numeric vector
#' @export
#'
#' @examples
#' 
#' rtetlg(n_dist = rgeom, rate = 2, p = 0.25, n = 6)

rtetlg <- function(n_dist, rate = 1, ...) {
  event_len <- n_dist(...)
  rand_sample <- do.call(
    rbind.data.frame,
    lapply(event_len, get_tetlg, rate)
  )
  rand_sample[['N']] <- event_len
  rand_sample
}

#' Simulate Geometric Distribution with Modified Probability of Success
#' 
#' @param n number of trials to simulate
#' @param p probability of success in each trial (0 <= p < 1)
#' @param success logical indicating whether or not to use a modified version of
#' the geometric distribution that returns ceiling(log(runif())/log(1-p)) 
#' instead of rgeom() when TRUE. Default is TRUE.
#' @return A vector of length n with random numbers drawn from the geometric 
#' distribution with probability p if success = FALSE, or modified geometric
#' distribution as described above otherwise. 
#' @export
#' 
#' @examples 
#' 
#' sim_geom <- r_geom(n = 100, p = .5, success=TRUE)


r_geom <- function(n,p, success = TRUE) {
  if (success) {
    ceiling(log(stats::runif(n))/log(1-p))
  } else {
    stats::rgeom(n, prob = p)
  }
}






