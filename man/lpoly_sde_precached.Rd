\name{lpoly_sde_precached}
\alias{lpoly_sde_precached}
\title{LPoly System Implicit SDE Simulator using NLOPT, with precached noise data}
\usage{
lpoly_sde_precached(sys, noise, start, from, to, steps, x_tol = 0,
  algorithm = "TNEWTON")
}
\arguments{
  \item{sys}{lpoly_system_type XPtr object created with
  lpoly_make_system}

  \item{noise}{Noise buffer (brownian motion normalised
  steps ~N(0,sigma^2))}

  \item{start}{Initial position: n vector}

  \item{from}{Initial time: scalar}

  \item{to}{Final time: scalar}

  \item{steps}{Number of points to take, s.t. dt = (from -
  to) / (steps + 1): integer}
}
\description{
Simulates a trajectory to the SDE specified by *sys* and
the noise level sigma, starting at the point start and
integrating over times [*from*,*to*] on *steps + 1* time
points.
}

