#' @param ord array of orders to include in the model
lpoly_model_spec <- function(d, ord){
    x <- data.frame(ord)
    s <- expand.grid(rep(x, d))
    s <- s[rowSums(s) %in% ord,]
    rownames(s) <- paste0("v", seq_len(nrow(s)))
    as.matrix(s)
}

examples.lorenz.sys <- function(s = 16, r = 45.6, b = 5){
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
