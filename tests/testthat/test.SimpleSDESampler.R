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

test_that( "synthetic.dataset produces expected structure", {
    print.noquote("Testing data structure...")
    with(data = list(tt = synthetic.dataset(num.entities = 2, tmax = 1, at.times = seq(0,1,0.01))),{
        expect_that(nrow(tt), equals(202))
        expect_that(colnames(tt), equals(c("entity", "time", "u.1", "u.2", "u.3")))
    })
})

test_that( "solve_implicit_sde produces correct solution to deterministic test equation", {
    print.noquote("Testing scalar test equation...")
    with(data = list(tt = solve_implicit_sde( d_det = function(u, t) -u
    	      			    , d_stoch = function(u, t) matrix(0, 1, 1)
				    , jacobian = function(u, t) matrix(-1)
				    , sigma = 0
				    , start = 1
				    , from = 0
				    , to = 1
				    , steps = 400)), {
        expect_less_than( abs(as.numeric(tail(tt, 1)) - exp(-1)), 0.001)
    })
})

test_that( "solve_implicit_sde_averages produces concistent estimates for a linear system", {
    # solution to du/dt = v, dv/dt = -u:
    # u = sin(t) + c1
    # v = cos(t) + c2
    # so with [u,v](0) = [0,1]
    # we expect E([u,v])(pi/2) = [1,0]
    
    print.noquote("Testing concistency on linear system...")
    
    A <- matrix(c(0, -1, 1, 0), 2, 2)
    with(data = list(res = tail( solve_implicit_sde_averages(
				 nrep = 100
    	      		     	 , d_det = function(u, t) t(A %*% t(u))
    	      			 , d_stoch = function(u, t) diag(rep(1,2))
				 , jacobian = function(u, t) A
				 , sigma = 0.1
				 , start = c(0,1)
				 , from = 0
				 , to = pi/2
				 , steps = 400), 1)), {
        expect_less_than(sqrt(sum((res - c(1,0))^2)), 0.05)
    })
})

test_that( "solve_implicit_sde uses different random seeds", {
    print.noquote("Testing random seeds...")
    with( data = list( x1 = solve_implicit_sde( d_det = function(u, t) -u
    	      			    , d_stoch = function(u, t) matrix(1)
				    , jacobian = function(u, t) matrix(-1)
				    , sigma = 0.1
				    , start = 1
				    , from = 0
				    , to = 1
				    , steps = 400), 
                       x2 = solve_implicit_sde( d_det = function(u, t) -u
    	      			    , d_stoch = function(u, t) matrix(1)
				    , jacobian = function(u, t) matrix(-1)
				    , sigma = 0.1
				    , start = 1
				    , from = 0
				    , to = 1
				    , steps = 400)), {
        expect_that(sum(as.numeric(x1 == x2)), equals(1))
    })
})
