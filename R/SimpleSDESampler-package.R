#' @title SimpleSDESampler.
#'
#' @name SimpleSDESampler
#' @docType package
#' 
#' @useDynLib SimpleSDESampler
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
#' \item \code{\link{lpoly_implicit_sde}}
#' \item \code{\link{lpoly_implicit_sde_averages}}
#' }