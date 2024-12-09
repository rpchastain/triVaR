% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/texp.R
\name{q_geom}
\alias{q_geom}
\title{Simulate Geometric Distribution with Modified Probability of Success}
\usage{
q_geom(p, ...)
}
\arguments{
\item{p}{the}

\item{...}{Optional parameters which are the same as `pgeom`}
}
\value{
A vector of length n with random numbers drawn from the geometric 
distribution with probability p if success = FALSE, or modified geometric
distribution as described above otherwise.
}
\description{
Simulate Geometric Distribution with Modified Probability of Success
}
\examples{

sim_geom <- q_geom(p = 0.25, prob = 0.25)
}
