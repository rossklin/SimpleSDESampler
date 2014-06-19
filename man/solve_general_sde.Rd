\name{solve_general_sde}
\alias{solve_general_sde}
\title{Simulates an SDE explicitly (derivative free)}
\usage{
solve_general_sde(d_det, d_stoch, start, from, to, steps)
}
\arguments{
  \item{d_det}{R function representing f(x)}

  \item{d_stoch}{R function representing g(x) E}

  \item{start}{numeric vector with initial state}

  \item{from}{scalar with initial time}

  \item{to}{scalar with final time}

  \item{steps}{number of steps to take}
}
\description{
For an Ito form SDE dx = f(x) dt + g(x) E sqrt(dt)
}

