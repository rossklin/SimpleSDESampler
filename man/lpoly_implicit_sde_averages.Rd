\name{lpoly_implicit_sde_averages}
\alias{lpoly_implicit_sde_averages}
\usage{
lpoly_implicit_sde_averages(nrep, sys, sigma, start, from, to, steps)
}
\arguments{
  \item{nrep}{Number of repetitions to average over:
  integer}

  \item{d_det}{Deterministic component: an R function: m x
  n matrix of m states -> m x n matrix of m time
  derivatives}

  \item{d_stoch}{Stochastic component: an R function: (1 x
  n matrix state, scalar time) -> 1 x n matrix state}

  \item{jacobian}{Jacobian of deterministic component: an R
  function: (n vector state, scalar time) -> n x n matrix
  df_i / du_j}

  \item{sigma}{Amplitude of noise: scalar}

  \item{start}{Initial position: n vector}

  \item{from}{Initial time: scalar}

  \item{to}{Final time: scalar}

  \item{steps}{Number of points to take, s.t. dt = (from -
  to) / (steps + 1): integer}
}

