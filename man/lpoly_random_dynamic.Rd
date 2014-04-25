\name{lpoly_random_dynamic}
\alias{lpoly_random_dynamic}
\title{LPoly random dynamics generator}
\usage{
lpoly_random_dynamic(d, ord = 1:2, scale.factor = 1, damp.factor = 0.2,
  pzero = 0, porder.decay = 0.5)
}
\arguments{
  \item{d}{Dimension of the system}

  \item{ord}{Polynomial orders to include: array}

  \item{scale.factor}{Scaling for normally distributed
  polynomial coefficients}

  \item{damp.factor}{Scaling for the damping terms}

  \item{pzero}{Probability of removing terms from the set
  of possible terms}

  \item{porder.decay}{Exponential reduction of coefficients
  w.r.t. order}
}
\description{
Generates an lpoly_system_type object representing a
polynomial dynamic with randomly generated coefficients.
}

