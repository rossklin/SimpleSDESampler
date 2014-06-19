\name{lpoly_compute_jacobian}
\alias{lpoly_compute_jacobian}
\title{LPoly Jacobian}
\usage{
lpoly_compute_jacobian(lps, state)
}
\arguments{
  \item{lps}{XPtr to an lpoly_system_type object created by
  a call to lpoly_make_system}

  \item{state}{Single row matrix}
}
\description{
Builds the model matrix from the data using the
lpoly_system_type object
}

