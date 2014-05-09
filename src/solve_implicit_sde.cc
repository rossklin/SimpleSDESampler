/* Software License Agreement (BSD License)
 *
 * Copyright (c) 2014, Ross Linscott (rossklin@gmail.com)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *
 *     The names of its contributors may not be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if defined(NDEBUG)
#undef NDEBUG
#endif

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstring>
#include <ctime>

#include <vector>
#include <utility>
#include <exception>

#include <Rcpp.h>

#include "../inst/include/SimpleSDESampler.h"


//' General implicit SDE simulator
//'
//' @param d_det Deterministic component: an R function: (m x n matrix of m states, scalar time) -> m x n matrix of m time derivatives
//' @param jacobian Jacobian of deterministic component: an R function: (1 x n matrix state, scalar time) -> n x n matrix df_i / du_j
//' @param sigma Amplitude of noise: scalar
//' @param start Initial position: n vector
//' @param from Initial time: scalar
//' @param to Final time: scalar
//' @param steps Number of points to take, s.t. dt = (from - to) / (steps + 1): integer
//' @export
// [[Rcpp::export]]
NumericMatrix solve_implicit_sde(Rcpp::Function d_det
				 , Rcpp::Function jacobian
				 , double sigma
				 , NumericVector start
				 , double from, double to, int steps
				 , double x_tol = 0
				 , const char* algorithm = "TNEWTON") {

  const double dt = (to - from)/steps;
  std::vector<double> state = as<std::vector<double> >(start);
  NumericMatrix result(steps+1, start.size());

  g_system s(&d_det, &jacobian);
  nlopt_stepper stepper(&s, dt, start.size(), sigma, steps, x_tol, algorithm);

  for(int j = 0; j < start.size(); ++j){
    result(0, j) = state[j];
  }

  for(int i = 1; i <= steps; ++i) {
    stepper.do_step(state, i*dt);
    for(int j = 0; j < start.size(); ++j){
      result(i, j) = state[j];
    }
  }

  return result;
}

//' General implicit SDE simulator: time-wise averages
//'
//' @param nrep Number of repetitions to average over: integer
//' @param d_det Deterministic component: an R function: m x n matrix of m states -> m x n matrix of m time derivatives
//' @param d_stoch Stochastic component: an R function: (1 x n matrix state, scalar time) -> 1 x n matrix state
//' @param jacobian Jacobian of deterministic component: an R function: (n vector state, scalar time) -> n x n matrix df_i / du_j
//' @param sigma Amplitude of noise: scalar
//' @param start Initial position: n vector
//' @param from Initial time: scalar
//' @param to Final time: scalar
//' @param steps Number of points to take, s.t. dt = (from - to) / (steps + 1): integer
//' @export
// [[Rcpp::export]]
NumericMatrix solve_implicit_sde_averages( int nrep
					   , Rcpp::Function d_det
					   , Rcpp::Function jacobian
					   , double sigma
					   , NumericVector start
					   , double from, double to, int steps
					   , double x_tol = 0
					   , const char* algorithm = "TNEWTON"){
  
  NumericMatrix result(steps+1, start.size());
  NumericMatrix buf;
  int i,j;
  double *p, *q;
  int smax = result.nrow() * result.ncol();

  p = &result(0,0);

  memset(p, 0, smax * sizeof(double));
  
  for (i = 0; i < nrep; i++){
    buf = solve_implicit_sde( d_det
			      , jacobian
			      , sigma
			      , start
			      , from
			      , to
			      , steps
			      , x_tol
			      , algorithm);

    q = &buf(0,0);
    for (j = 0; j < smax; j++){
      p[j] += q[j];
    }
  }

  for (j = 0; j < smax; j++){
    p[j] /= nrep;
  }

  return result;
} 
