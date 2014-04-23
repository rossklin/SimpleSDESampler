\name{lpoly_implicit_sde}
\alias{lpoly_implicit_sde}
\title{LPoly System Implicit SDE Simulator}
\usage{
lpoly_implicit_sde(sys, sigma, start, from, to, steps)
}
\arguments{
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
Simulates a trajectory to the SDE specified by *sys* and
the noise level sigma, starting at the point start and
integrating over times [*from*,*to*] on *steps + 1* time
points.
}

