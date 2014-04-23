\name{synthetic.dataset.quick}
\alias{synthetic.dataset.quick}
\title{Generate Synthetic Dataset: optimised for Laurent polynomial systems}
\usage{
synthetic.dataset.quick(num.entities = 10, tmax = 10, steps = 400 * tmax,
  process.noise.sd = 0, observation.noise.sd = 0.01, do.standardise = F,
  initial.generator = function(i) {     rnorm(3) },
  sys = examples.lpsys.lorenz(), at.times = seq(0, tmax, length.out = 1001),
  save.to = NULL)
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

  \item{sys}{lpoly system constructed with
  lpoly_make_system, see examples.lpsys.lorenz()}

  \item{at.times}{Sequence of times to include in the
  output.

  Produces sample trajectories for an SDE on Ito form:
  dx(t) = f(x(t)) dt + g(x(t), t) e(t) sqrt(dt) where
  det.deriv is f and stoch.deriv is g}
}
\description{
Simulate sample trajectories of an SDE and store the
results in a time.table
}

