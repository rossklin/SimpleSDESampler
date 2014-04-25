#' LPoly model spec builder
#' @param d number of measurement variables
#' @param ord array of polynomial orders to include in the model
#' @export

lpoly_model_spec <- function(d, ord){
    x <- data.frame(0:max(ord))
    s <- expand.grid(rep(x, d))
    s <- s[rowSums(abs(s)) %in% ord,]
    rownames(s) <- paste0("T", seq_len(nrow(s)))
    colnames(s) <- paste0("V", seq_len(d))
    as.matrix(s)
}

#' LPoly random dynamics generator
#' 
#' Generates an lpoly_system_type object representing a polynomial dynamic 
#' with randomly generated coefficients.
#' @param d Dimension of the system
#' @param ord Polynomial orders to include: array
#' @param scale.factor Scaling for normally distributed polynomial coefficients
#' @param damp.factor Scaling for the damping terms
#' @param pzero Probability of removing terms from the set of possible terms
#' @param porder.decay Exponential reduction of coefficients w.r.t. order
#' @export

lpoly_random_dynamic <- function(d, ord = 1:2, scale.factor = 1, damp.factor = 0.2, pzero = 0, porder.decay = 0.5){
    next_odd <- function(x) 2 * as.integer(x / 2) + 1

    ms <- lpoly_model_spec(d, ord)
    ms <- ms[runif(nrow(ms), 0, 1) > pzero,]
    mspec <- rbind(ms, diag(rep(next_odd(max(ord) + 1), d)))
    A <- matrix(scale.factor * rnorm(nrow(ms) * d), d, nrow(ms))
    A <- t(aaply(seq_len(ncol(A)), 1, function(cidx) porder.decay^(sum(ms[cidx,]) - 1) * A[,cidx]))
    B <- diag(rep(-damp.factor, d))
    mcm <- cbind(A, B)
    lsys <- lpoly_make_system(mcm, mspec)
}

#' LPoly example dynamic: lorenz system
#' @export

examples.lpsys.lorenz <- function(s = 16, r = 45.6, b = 5){
  cm = t(matrix(c( -s, s, 0, 0, 0, 0,
       	  	      r, -1, 0, 0, -20, 0,
		      0, 0, -b, 5, 0, 0), ncol = 3))
  trm = t(matrix(c( 1,0,0,
    		    0,1,0,
		    0,0,1,
		    1,1,0,
		    1,0,1,
		    0,1,1), nrow = 3))
  
  lpoly_make_system(cm, trm)
}

examples.gensys.det.lorenz <- function(u,t){
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

attr(examples.gensys.det.lorenz, "description") <- list("dx/dt = 16 y - 16 x", "dy/dt = 45.6 x - y - 20 x z", "dz/dt = 5 x y - 5 z")

examples.gensys.jacob.lorenz <- function(u_,t){
    u <- as.numeric(u_)
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
