
synthetic.dataset <- function( num.entities = 10
                              , tmax = 10
                              , steps = 1000 * tmax
                              , process.noise.sd = 0
                              , observation.noise.sd = .01
                              , do.standardise = F
                              , initial.generator = function(i){
                                  rnorm(3)
                              }
                              , det.deriv = det.lorenz
                              , stoch.deriv = function(t) concistent.lindep.noise(t, seq(0, tmax, length.out = steps + 1), 3)
			      , jacobian = jacob.lorenz
                              , at.times = seq(0, tmax, length.out = 1001)){

    dimension = length(initial.generator(0))

    df <- adply(seq_len(num.entities), 1, function(i) data.frame( time = seq(0, tmax, length.out = steps + 1)
                                                                 , u = solve_implicit_sde( det.deriv
                                                                       , stoch.deriv
								       , jacobian
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

    dt
}

det.lorenz <- function(u){
    ## lorenz system proposed by Chi-Chung Chen and Kung Yao
    ## in STOCHASTIC CALCULUS NUMERICAL EVALUATION OF CHAOTIC COMMUNICATION SYSTEM PERFORMANCE
    ## note: behaviour seems stable with respect to each parameter value
    sigma <- 16
    r <- 45.6
    b <- 5
    
    x <- u[1]
    y <- u[2]
    z <- u[3]

    du <- c(sigma * (y - x), r * x - y - 20 * x * z, 5 * x * y - b * z)
}

jacob.lorenz <- function(u){
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

det.linear2d <- function(u){
    du <- c(-u[2], u[1])
}

concistent.lindep.noise <- function(t, tlist, d){
  set.seed(which.min(abs(tlist - t)))
  rnorm(d)
}