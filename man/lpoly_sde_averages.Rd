\name{lpoly_sde_averages}
\alias{lpoly_sde_averages}
\title{LPoly System Implicit SDE Simulator: time-wise averages}
\usage{
lpoly_sde_averages(nrep, sys, sigma, start, from, to, steps)
}
\arguments{
  \item{nrep}{Number of repetitions to average over:
  integer}

  \item{sys}{lpoly_system_type XPtr object created with
  lpoly_make_system}

  \item{sigma}{Amplitude of noise: scalar}

  \item{start}{Initial position: n vector}

  \item{from}{Initial time: scalar}

  \item{to}{Final time: scalar}

  \item{steps}{Number of points to take, s.t. dt = (from -
  to) / (steps + 1): integer}
}
\description{
Samples multiple times from lpoly_sde and returns the
time-wise average
}

