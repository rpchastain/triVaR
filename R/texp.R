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
  if (rate == 0) {
    runif(n, max = delta)
  }
  
  u <- runif(n)
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
  rand_sample <- rtexp(...)
  data.frame(x = sum(rand_sample), y = max(rand_sample))
}

get_tetlg <- function(...) {
  rand_sample <- rexp(...)
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

#' Simulate Geometric Distribution defined in n = 1, 2,...
#' 
#' @param n Number of trials to simulate
#' match the duration of episodes with a minimum duration
#' @param ... Optional parameters which are the same as `rgeom`
#' @return A vector of length n with random numbers drawn from the geometric 
#' distribution with probability p if success = FALSE, or modified geometric
#' distribution as described above otherwise. 
#' @export
#' 
#' @examples 
#' 
#' sim_geom <- r_geom(n = 100, p = .5)

r_geom <- function(n, ...){ rgeom(n, ...) + 1}

#' Returns the distribution value at `x`
#' 
#' @param x number of trials to simulate
#' @param ... Optional parameters which are the same as `dgeom`
#' @return A vector of length n with random numbers drawn from the geometric 
#' distribution with probability p if success = FALSE, or modified geometric
#' distribution as described above otherwise. 
#' @export
#' 
#' @examples 
#' 
#' dist_geom <- d_geom(4, prob = 0.25)

d_geom <- function(x, ...){ dgeom(x - 1, ...) }

#' Simulate Geometric Distribution with Modified Probability of Success
#' 
#' @param q probability of the value `q`
#' @param ... Optional parameters which are the same as `pgeom`
#' @return A vector of length n with random numbers drawn from the geometric 
#' distribution with probability p if success = FALSE, or modified geometric
#' distribution as described above otherwise. 
#' @export
#' 
#' @examples 
#' 
#' sim_geom <- p_geom(q = 0.75, p = .5)

p_geom <- function(q, ...){ pgeom(q - 1, ...) }

#' Simulate Geometric Distribution with Modified Probability of Success
#' 
#' @param p the 
#' @param ... Optional parameters which are the same as `pgeom`
#' @return A vector of length n with random numbers drawn from the geometric 
#' distribution with probability p if success = FALSE, or modified geometric
#' distribution as described above otherwise. 
#' @export
#' 
#' @examples 
#' 
#' sim_geom <- q_geom(p = 0.25, prob = 0.25)

q_geom <- function(p, ...){ qgeom(p, ...) + 1 }


#' @title Rate Maximum Likelihood Estimation (MLE) Function This function 
#' calculates the MLE for the rate parameter of the MSTE-N 
#' distribution. This function has no analytical solution so numerical methods
#' are required to find an approximate solution.
#' 
#' @param rate The rate parameter of a Truncated Exponential distribution.
#' @param b The upper truncation parameter.
#' @param samp A data frame containing two columns: "N" (typically the duration
#' of the Trivariate Event) and "x" (the sum of the values measured during the
#' Trivariate Event).
#'
#' @return The estimate of the rate parameter.
#' @export
#'
#' @examples \dontrun{
#' rate_mle(
#'   rate = 0.5,
#'   b = 100, 
#'   samp = data.frame(
#'     N = c(rep(20, 10), rep(30, 10)),
#'     x = rnorm(20)
#'     )
#'   )
#'}


rate_mle <- function(rate, b, samp) {
  first <- mean(samp$N)/rate
  second <- mean(samp$x)
  third <- mean(samp$N) * b /(exp(rate * b) - 1)
  abs(first - second - third)
}

#' @title Function to optimize the `rate_mle()` function. This is a convenience
#' function  which use the stats `R` function `optimize()`
#' to numerically find the MLE of the MSTE-N distribution.
#' 
#' @param rate_interval Initial interval for estimating the rate parameter.
#' Defaults to c(-5, 5).
#' @param optim_tries Maximum number of optimization attempts. Defaults to 25.
#' @param b The upper truncation parameter to pass to `rate_mle()`.
#' @param samp The data sample to pass to `rate_mle()`.
#'
#' @return The estimated rate parameter as a single value or a warning message 
#' if no solution is found within the specified interval and number of 
#' optimization tries.
#' 
#' @export
#'
#' @examples \dontrun{
#' rate_estimator(rate_interval = c(-3, 4), optim_tries = 50)
#' }

rate_estimator <- function(
    rate_interval = c(-5, 5), optim_tries = 25,
    b, samp
) {
  
  while (optim_tries > 0) {
    beta_hat <- optimize(
      rate_mle, b = b, samp = samp,
      interval = rate_interval, maximum = FALSE
    )$minimum
    
    beta_edges <- abs(beta_hat - rate_interval)
    if (any(abs(beta_hat - rate_interval) < 1e-4)) {
      rate_interval <- c(beta_hat - 2, beta_hat + 2)
      optim_tries <- optim_tries - 1
      no_solution <- TRUE
    } else {
      optim_tries <- 0
      no_solution <- FALSE
    }
  }
  if (no_solution) {
      warning_msg <- sprintf(
        "Convergence for estimate of rate not found in %s rounds of 
         optimization. Please change the rate_interval or increase the number
         of optimization rounds. Final estimate for the rate: %s",
        optim_tries, beta_hat
      )
      warning(strwrap(warning_msg))
  } else {
    beta_hat
  }
}






