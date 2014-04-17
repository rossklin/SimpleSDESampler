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
#' @param det.deriv Function computing the deterministic component of the time derivative given the state.
#' @param stoch.deriv Function computing the stochastic component of the time derivative given the state.
#' @param jacobian Function computing the jacobian of det.deriv with respect to the state variables.
#' @param at.times Sequence of times to include in the output.
#'
#' Produces sample trajectories for an SDE on Ito form: 
#' dx(t) = f(x(t)) dt + g(x(t), t) e(t) sqrt(dt)
#' where det.deriv is f and stoch.deriv is (at the moment) g * e
#' (probably will change so stoch.deriv is only g in the future)
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
                              , det.deriv = det.lorenz
                              , stoch.deriv = function(x, t) diag(rep(1,3))
			      , jacobian = jacob.lorenz
                              , at.times = seq(0, tmax, length.out = 1001)
			      , save.to = NULL){

    dimension = length(initial.generator(0))

    df <- adply(seq_len(num.entities), 1, function(i) data.frame( time = seq(0, tmax, length.out = steps + 1)
                                                                 , u = solve_implicit_sde( det.deriv
                                                                       , stoch.deriv
								       , jacobian
								       , process.noise.sd
                                                                       , initial.generator(i)
                                                                       , 0
                                                                       , tmax
                                                                       , steps)), .progress = "text")
    colnames(df)[[1]] <- "entity"
    df[,c(-1,-2)] <- df[,c(-1,-2)] + matrix(rnorm(dimension * num.entities * (steps+1), 0, observation.noise.sd), num.entities * (steps + 1), dimension)
    if (do.standardise) df[,c(-1,-2)] <- sapply(df[,c(-1,-2)], standardise)
    ## tt <- as.time.table(df, "entity", "time", paste0("u.", seq_len(dimension)))
    ## if (!is.null(at.times)) tt <- subset(tt, times = data.table(at.times), index = data.table(unique(index(tt))))
    ## tt
    ## for now, don't use time.table
    
    dt <- data.table(df, key = c("entity", "time"))

    if (!is.null(at.times)) dt <- dt[time %in% at.times]

    if (!is.null(save.to)){
      write.csv(dt, save.to, row.names = F)
    }

    dt
}

det.lorenz <- function(u,t){
    ## lorenz system proposed by Chi-Chung Chen and Kung Yao
    ## in STOCHASTIC CALCULUS NUMERICAL EVALUATION OF CHAOTIC COMMUNICATION SYSTEM PERFORMANCE
    ## note: behaviour seems stable with respect to each parameter value
    sigma <- 16
    r <- 45.6
    b <- 5
    
    x <- u[,1]
    y <- u[,2]
    z <- u[,3]

    du <- cbind(sigma * (y - x), r * x - y - 20 * x * z, 5 * x * y - b * z)
}

attr(det.lorenz, "description") <- list("dx/dt = 16 y - 16 x", "dy/dt = 45.6 x - y - 20 x z", "dz/dt = 5 x y - 5 z")

jacob.lorenz <- function(u,t){
    sigma <- 16
    r <- 45.6
    b <- 5
    
    x <- u[1]
    y <- u[2]
    z <- u[3]

    matrix(c( -sigma, sigma, 0 
    	      , r - 20 * z, -1, -20 * x
	      , 5 * y, 5 * x, -b), nrow = 3, ncol = 3)
}

det.linear2d <- function(u,t){
    du <- c(-u[2], u[1])
}

#concistent.lindep.noise <- function(ntimes, d, sd){
#  n <- ntimes
#  stoch.data <- matrix(rnorm(n * d, 0, sd), n, d)
#  function(idx){
#    stoch.data[idx %% n,]
#  }
#}