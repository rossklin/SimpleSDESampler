\name{synthetic.dataset}
\alias{synthetic.dataset}
\title{Generate Synthetic Dataset}
\usage{
synthetic.dataset(num.entities = 10, tmax = 10, steps = 400 * tmax,
  process.noise.sd = 0, observation.noise.sd = 0.01, do.standardise = F,
  initial.generator = function(i) {     rnorm(3) },
  det.deriv = examples.gensys.det.lorenz, stoch.deriv = function(x, t)
  diag(rep(1, 3)), jacobian = examples.gensys.jacob.lorenz,
  at.times = seq(0, tmax, length.out = 1001), save.to = NULL,
  retries = 10)
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

  \item{det.deriv}{Function f: (m x n matrix of m states,
  scalar t time) -> m x n matrix of m state derivs
  computing the deterministic component of the dynamic
  given the state.}

  \item{stoch.deriv}{Function g: (1 x n matrix state,
  scalar t time) -> n x n matrix of noise weights,
  computing the stochastic coefficient of the time dynamic
  the state.}

  \item{jacobian}{Function J: (1 x n matrix state, scalar t
  time) -> n x n matrix d fi / d xj computing the jacobian
  of det.deriv with respect to the state variables.}

  \item{at.times}{Array of times to include in the output.

  Produces sample trajectories for an SDE on Ito form:
  dx(t) = f(x(t), t) dt + g(x(t), t) e(t) sqrt(dt) where
  det.deriv is f and stoch.deriv is g}
}
\description{
Simulate sample trajectories of an SDE and store the
results in a time.table
}

