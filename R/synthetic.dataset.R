## Software License Agreement (BSD License)
##
## Copyright (c) 2014, Ross Linscott (rossklin@gmail.com)
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
##     The names of its contributors may not be used to endorse or promote products
##     derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' Generate Synthetic Dataset
#'
#' Simulate sample trajectories of an SDE and store the results in a time.table
#'
#' @param num.entities Number of sample trajectories to generate.
#' @param tmax End of the time interval to integrate the SDE over.
#' @param steps Number of time steps to use in the integration.
#' @param process.noise.sd Standard deviation of the brownian motion component.
#' @param observation.noise.sd Standard deviation of synthetic observation noise.
#' @param do.standardise Standardise the output?
#' @param initial.generator Function that takes an index and generates a starting point for a sample SDE trajectory.
#' @param det.deriv Function f: (m x n matrix of m states, scalar t time) -> m x n matrix of m state derivs computing the deterministic component of the dynamic given the state.
#' @param stoch.deriv Function g: (1 x n matrix state, scalar t time) -> n x n matrix of noise weights, computing the stochastic coefficient of the time dynamic the state.
#' @param jacobian Function J: (1 x n matrix state, scalar t time) -> n x n matrix d fi / d xj computing the jacobian of det.deriv with respect to the state variables.
#' @param at.times Array of times to include in the output.
#'
#' Produces sample trajectories for an SDE on Ito form: 
#' dx(t) = f(x(t), t) dt + g(x(t), t) e(t) sqrt(dt)
#' where det.deriv is f and stoch.deriv is g
#' @export

synthetic.dataset <- function( num.entities = 10
                              , tmax = 10
                              , steps = 400 * tmax
                              , process.noise.sd = 0
                              , observation.noise.sd = .01
                              , do.standardise = F
                              , initial.generator = function(i){
                                  rnorm(3)
                              }
                              , det.deriv = examples.gensys.det.lorenz
                              , stoch.deriv = function(x, t) diag(rep(1,3))
			      , jacobian = examples.gensys.jacob.lorenz
                              , at.times = seq(0, tmax, length.out = 1001)
			      , save.to = NULL
                              , retries = 10){

    if (!isTRUE(all.equal(at.times * steps / tmax, as.integer(at.times * steps / tmax + 0.5)))) stop("Times don't match!")

    dimension = length(initial.generator(0))

    df <- adply(seq_len(num.entities), 1, function(i){
        rt <- 0

        while (rt < retries){
            u = tryCatch( solve_implicit_sde( d_det = det.deriv
                , d_stoch = stoch.deriv
                , jacobian = jacobian
                , sigma = process.noise.sd
                , start = initial.generator(i)
                , from = 0
                , to = tmax
                , steps = steps)
                , error = function(e) NULL
                )
            rt <- if (is.null(u)) rt + 1 else retries
        }

        if (is.null(u)) u <- matrix(NaN, steps + 1, length(initial.generator(i)))
        
        data.frame( time = seq(0, tmax, length.out = steps + 1), u = u)
    }, .progress = "text")
    
    colnames(df)[[1]] <- "entity"
    df[,c(-1,-2)] <- df[,c(-1,-2)] + matrix(rnorm(dimension * num.entities * (steps+1), 0, observation.noise.sd), num.entities * (steps + 1), dimension)
    if (do.standardise) df[,c(-1,-2)] <- sapply(df[,c(-1,-2)], standardise)
    tt <- as.time.table(df, "entity", "time")
    if (!is.null(at.times)) tt <- subset(tt, times = at.times, index = unique(index(tt)))
    
    if (!is.null(save.to)){
      write.csv(tt, save.to, row.names = F)
    }
    tt
}

#' Generate Synthetic Dataset: optimised for Laurent polynomial systems
#'
#' Simulate sample trajectories of an SDE and store the results in a time.table
#'
#' @param num.entities Number of sample trajectories to generate.
#' @param tmax End of the time interval to integrate the SDE over.
#' @param steps Number of time steps to use in the integration.
#' @param process.noise.sd Standard deviation of the brownian motion component.
#' @param observation.noise.sd Standard deviation of synthetic observation noise.
#' @param do.standardise Standardise the output?
#' @param initial.generator Function that takes an index and generates a starting point for a sample SDE trajectory.
#' @param sys lpoly system constructed with lpoly_make_system, see examples.lpsys.lorenz()
#' @param at.times Sequence of times to include in the output.
#'
#' Produces sample trajectories for an SDE on Ito form: 
#' dx(t) = f(x(t)) dt + g(x(t), t) e(t) sqrt(dt)
#' where det.deriv is f and stoch.deriv is g
#' @export

synthetic.dataset.quick <- function( num.entities = 10
                                    , tmax = 10
                                    , steps = 400 * tmax
                                    , process.noise.sd = 0
                                    , observation.noise.sd = .01
                                    , do.standardise = F
                                    , initial.generator = function(i){
                                        rnorm(3)
                                    }
                                    , sys = examples.lpsys.lorenz()
                                    , at.times = seq(0, tmax, length.out = 1001)
                                    , save.to = NULL
                                    , retries = 10){

    if (!isTRUE(all.equal(at.times * steps / tmax, as.integer(at.times * steps / tmax + 0.5)))) stop("Times don't match!")

    dimension = length(initial.generator(0))

    df <- adply(seq_len(num.entities), 1, function(i){
        rt <- 0

        while (rt < retries){
            u = tryCatch( lpoly_implicit_sde(sys = sys
                , sigma = process.noise.sd
                , start = initial.generator(i)
                , from = 0
                , to = tmax
                , steps = steps)
                , error = function(e) NULL
                )
            rt <- if (is.null(u)) rt + 1 else retries
        }

        if (is.null(u)) u <- matrix(NaN, steps + 1, length(initial.generator(i)))
        
        data.frame( time = seq(0, tmax, length.out = steps + 1), u = u)
    }, .progress = "text")
    
    colnames(df)[[1]] <- "entity"
    df[,c(-1,-2)] <- df[,c(-1,-2)] + matrix(rnorm(dimension * num.entities * (steps+1), 0, observation.noise.sd), num.entities * (steps + 1), dimension)
    if (do.standardise) df[,c(-1,-2)] <- sapply(df[,c(-1,-2)], standardise)
    tt <- as.time.table(df, "entity", "time")
    if (!is.null(at.times)) tt <- subset(tt, times = at.times, index = unique(index(tt)))

    if (nrow(na.omit(tt)) != nrow(tt)) warning("Non-numeic values introduced, may be caused by too long time steps or by a mismatch between generation and output times")
    
    if (!is.null(save.to)){
      write.csv(tt, save.to, row.names = F)
    }
    tt
}

