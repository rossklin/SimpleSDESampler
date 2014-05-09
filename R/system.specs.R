#' LPoly system builder
#' @param mcm model coefficient matrix
#' @param mspec model specification, as created by lpoly_model_spec

lpoly_make_system <- function(mcm, mspec){
    lsys <- lpoly_make_system_xptr(mcm, mspec)
    attr(lsys, "mcm") <- mcm
    attr(lsys, "mspec") <- mspec
    lsys
}

#' LPoly system modelfun
#' @param sys lpoly_system_type object created with lpoly_make_system

lpoly_system_modelfun <- function(sys){
    function(x){
        m <- lpoly_model_matrix(sys, x)
        colnames(m) <- rownames(attr(sys, "mspec"))
        m
    }
}

#' LPoly system evalfun
#' @param sys lpoly_system_type object created with lpoly_make_system

lpoly_system_evalfun <- function(sys){
    function(x) {
        m <- lpoly_model_matrix(sys, x) %*% t(attr(sys, "mcm"))
        colnames(m) <- rownames(attr(sys, "mcm"))
        m
    }
}

#' LPoly system jacobian
#' @param sys lpoly_system_type object created with lpoly_make_system

lpoly_system_jacobian <- function(sys){
    function(x) lpoly_compute_jacobian(sys, x)
}

#' LPoly system specs
#' @param sys lpoly_system_type object created with lpoly_make_system

lpoly_system_specs <- function(sys){
    res <- NULL
    res$mcm <- attr(sys, "mcm")
    res$mspec <- attr(sys, "mspec")
    res
}

#' LPoly model spec builder
#' @param d number of measurement variables
#' @param ord array of polynomial orders to include in the model
#' @export

lpoly_model_spec <- function(d, ord){
    x <- data.frame(unique(c(0:max(ord), ord)))
    s <- expand.grid(rep(x, d))
    s <- s[rowSums(abs(s)) %in% abs(ord),,drop=F]
    rownames(s) <- alply(s, 1, function(v) paste0(v, collapse = "."))
    colnames(s) <- paste0("V", seq_len(d))
    as.matrix(s)
}

#' LPoly steplength guesser: sampling method
#' @export

lpoly_guess_steplength_sample <- function(sys, goal = 0.1, freq = 400, h0 = 1){
    h <- 2 * h0
    ml2 <- 1
    steps <- 10
    f <- lpoly_system_evalfun(sys)

    while(ml2 > goal){
        h <- h / 2
        print.noquote(sprintf("Testing h = %e...", h))
        tt <- synthetic.dataset.quick(sys = sys
                                      , num.entities = 10
                                      , steps = steps * as.integer(freq * h + 1)
                                      , tmax = steps * h
                                      , at.times = seq(0, steps*h, h)
                                      , observation.noise.sd = 0
                                      , process.noise.sd = 0)
        
        z <- as.matrix(measurement(step_derivatives(tt)))
        x <- as.matrix(measurement(tt))
        y <- f(x)
        r <- y - z
        l2e <- sqrt(rowSums(r^2) / rowSums(y^2))
        l2e[is.na(l2e)] <- 1
        ml2 <- median(l2e)
        print.noquote(sprintf("ml2 = %f", ml2))
    }

    h
}

#' LPoly steplength guesser: linear method
#' @export

lpoly_guess_steplength_linear <- function(sys){
    f <- lpoly_system_evalfun(sys)
    j <- lpoly_system_jacobian(sys)
    d <- nrow(lpoly_system_specs(sys)$mcm)
    n <- 100
    x <- matrix(rnorm(n * d), n, d)
    js <- alply(1:n, 1, function(x) j(matrix(x, 1, d)))
    ls <- laply(js, function(m) eigen(m)$values)
    ys <- aaply(1:n, 1, function(i) log(2 * f(matrix(x[i,], 1, d)) / min(abs(ls[i,]))^2) / min(abs(ls[i,])))
    h <- 0.1 * quantile(ys, 0.1)
    h
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

lpoly_random_dynamic <- function(d, ord = 1:2, scale.factor = 1, damp.factor = 0.02, pzero = 0.8, model.size = NULL, porder.decay = 0.5){
    next_odd <- function(x) 2 * as.integer(x / 2) + 1

    ms <- lpoly_model_spec(d, ord)

    if (!(nrow(ms) > 0)) stop("Can not generate dynamic without model spec!")

    if (is.null(model.size)){
        use.rows <- 0
        while (sum(use.rows) < 1) use.rows <- which(runif(nrow(ms), 0, 1) > pzero)
    }else {
        use.rows <- sample(seq_len(nrow(ms)), model.size)
    }
    
    ms <- ms[use.rows,,drop=F]
    msB <- diag(rep(next_odd(max(ord) + 1), d))
    msC <- diag(rep(-next_odd(max(-ord) + 1), d))
    mspec <- ms
    msnames <- rownames(ms)
    if (max(ord) > 0) {
        mspec <- rbind(mspec, msB)
        msnames <- c(msnames, paste0("neg.feed.", seq_len(d)))
    }
    if (min(ord) < 0) {
        mspec <- rbind(mspec, msC)
        msnames <- c(msnames, paste0("pos.feed.", seq_len(d)))
    }
    rownames(mspec) <- msnames
    colnames(mspec) <- paste0("u.", seq_len(d))
    
    A <- matrix(scale.factor * rnorm(nrow(ms) * d), d, nrow(ms))

    ## Yes, R is wonderfully concistent, that's why we have to write
    ## special case clauses for different sizes of arrays...

    if (ncol(A) > 1){
        A <- t(aaply(seq_len(ncol(A)), 1, function(cidx) porder.decay^(sum(abs(ms[cidx,])) - 1) * A[,cidx]))
    }else{
        A <- porder.decay^sum(abs(ms)) * A
    }

    B <- diag(rep(-damp.factor, d))
    C <- diag(rep(damp.factor, d))
    mcm <- A
    if (max(ord) > 0) mcm <- cbind(mcm, B)
    if (min(ord) < 0) mcm <- cbind(mcm, C)
    rownames(mcm) <- colnames(mspec)
    colnames(mcm) <- rownames(mspec)
    
    lsys <- lpoly_make_system(mcm, mspec)
}

#' LPoly example dynamic: lorenz system
#' @export

lpoly_examples_lorenz <- function(s = 16, r = 45.6, b = 5){
  cm = t(matrix(c( -s, s, 0, 0, 0, 0,
       	  	      r, -1, 0, 0, -20, 0,
		      0, 0, -b, 5, 0, 0), ncol = 3))
  trm = t(matrix(c( 1,0,0,
    		    0,1,0,
		    0,0,1,
		    1,1,0,
		    1,0,1,
		    0,1,1), nrow = 3))

  rownames(trm) <- c("x", "y", "z", "xy", "xz", "yz")
  colnames(trm) <- c("dx", "dy", "dz")
  rownames(cm) <- colnames(trm)
  colnames(cm) <- rownames(trm)
  
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

    t(matrix(c( -sigma, sigma, 0 
               , r - 20 * z, -1, -20 * x
               , 5 * y, 5 * x, -b), nrow = 3, ncol = 3))
}
