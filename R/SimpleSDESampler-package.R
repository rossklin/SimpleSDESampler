#' @title SimpleSDESampler.
#'
#' @name SimpleSDESampler
#' @docType package
#' 
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp
#' @import plyr 
#' @import timetablr
#' @import data.table
#' @description
#' Runs simulations on SDE systems
#' @section Included functions:
#' \itemize{
#' \item \code{\link{synthetic.dataset}}
#' \item \code{\link{synthetic.dataset.quick}}
#' \item \code{\link{lpoly_sde}}
#' \item \code{\link{lpoly_sde_precached}}
#' \item \code{\link{lpoly_sde_averages}}
#' \item \code{\link{solve_general_sde}}
#' \item \code{\link{solve_implicit_sde}}
#' \item \code{\link{solve_implicit_sde_averages}}
#' \item \code{\link{lpoly_make_system}}
#' \item \code{\link{lpoly_system_modelfun}}
#' \item \code{\link{lpoly_system_evalfun}}
#' \item \code{\link{lpoly_system_jacobian}}
#' \item \code{\link{lpoly_system_specs}}
#' \item \code{\link{lpoly_model_spec}}
#' \item \code{\link{lpoly_random_dynamic}}
#' \item \code{\link{lpoly_examples_lorenz}}
#' }

.onLoad <- function(libloc, pkgname){
    dyn.load(file.path(libloc, pkgname, "libs", paste0("libnlopt",.Platform$dynlib.ext)))
    library.dynam("SimpleSDESampler", pkgname, libloc)
}
