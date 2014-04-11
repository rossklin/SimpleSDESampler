\name{synthetic.dataset}
\alias{synthetic.dataset}
\title{Generate Synthetic Dataset}
\usage{
synthetic.dataset(num.entities = 10, tmax = 10, steps = 400 * tmax,
  process.noise.sd = 0, observation.noise.sd = 0.01, do.standardise = F,
  initial.generator = function(i) {     rnorm(3) }, det.deriv = det.lorenz,
  stoch.deriv = function(x, t) diag(rep(1, 3)), jacobian = jacob.lorenz,
  at.times = seq(0, tmax, length.out = 1001), save.to = NULL)
}
\arguments{
  \item{num.entities}{Number of sample trajectories to
  generate.}

  \item{tmax}{End of the time interval to integrate the SDE
  over.}

  \item{steps}{Number of time steps to use in the
  integration.}

  \item{process.noise.sd}{Standard deviation of the
  brownian motion component.}

  \item{observation.noise.sd}{Standard deviation of
  synthetic observation noise.}

  \item{do.standardise}{Standardise the output?}

  \item{initial.generator}{Function that takes an index and
  generates a starting point for a sample SDE trajectory.}

  \item{det.deriv}{Function computing the deterministic
  component of the time derivative given the state.}

  \item{stoch.deriv}{Function computing the stochastic
  component of the time derivative given the state.}

  \item{jacobian}{Function computing the jacobian of
  det.deriv with respect to the state variables.}

  \item{at.times}{Sequence of times to include in the
  output.

  Produces sample trajectories for an SDE on Ito form:
  dx(t) = f(x(t)) dt + g(x(t), t) e(t) sqrt(dt) where
  det.deriv is f and stoch.deriv is (at the moment) g * e
  (probably will change so stoch.deriv is only g in the
  future)}
}
\description{
Simulate sample trajectories of an SDE and store the
results in a time.table
}

