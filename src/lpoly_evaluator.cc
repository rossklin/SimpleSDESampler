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
#include <cmath>
#include <cstring>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#ifdef STANDALONE
#include <RInside.h>
#endif

#include "../inst/include/SimpleSDESampler.h"
  
lpoly_evaluator::lpoly_evaluator(NumericMatrix cm, NumericMatrix trm) : terms(trm), coef_matrix(as_ublas_matrix(cm)){}

void lpoly_evaluator::set_h(double v){
  h = v;
}

NumericMatrix lpoly_evaluator::build(NumericMatrix data){
  double value;
  NumericMatrix result(data.nrow(), terms.nrow());
  int i,j,k;

  for (i = 0; i < data.nrow(); i++){
    for (j = 0; j < terms.nrow(); j++){
      value = 1;
      for (k = 0; k < terms.ncol(); k++){
	value *= pow(data(i, k), terms(j,k));
      }
      result(i,j) = value;
    }

  }
  return result;
}

// compute evaluation function
void lpoly_evaluator::operator()(uvector &query, uvector &out, double t) {
  if (query.size() != terms.ncol()){
    Rcpp::Rcout << "lpoly_evaluator::(): dimension mismatch: " << query.size() << " != " << terms.ncol() << endl;
    exit(-1);
  }

  NumericMatrix state(1, query.size(), query.begin());
  NumericMatrix state_model = build(state);
  uvector state_vector = as_ublas_vector(NumericVector(state_model.begin(), state_model.end()));

  uvector dpart = prod(coef_matrix, state_vector);
  out = dpart + (1 / sqrt(h)) * noise(t);
}

uvector lpoly_evaluator::noise(double t){
  int idx = t / h;
  int dim = noise_precache.size2();
  uvector x(dim);

  if (idx >= noise_precache.size1()){
    Rcpp::Rcout << "Invalid time index: " << idx << endl;
    exit(1);
  }

  memcpy(&x[0], &noise_precache(idx,0), dim * sizeof(double));
  return x;
}

void lpoly_evaluator::set_noise(umatrix m){
  noise_precache = m;
}

void lpoly_evaluator::build_noise(int n, int d, double sigma){
  int i;
  int size = n * d;
  double *p;

  boost::mt19937 gener;
  boost::normal_distribution<> normal(0,sigma*sigma);
  boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal);

  boost::posix_time::time_duration diff_time = boost::posix_time::microsec_clock::local_time() - boost::posix_time::second_clock::local_time();
  unsigned int the_seed = diff_time.total_microseconds();

  rng.engine().seed(the_seed);
  rng.distribution().reset();

  noise_precache.resize(n, d);
  p = &noise_precache(0,0);
    
  for (i = 0; i < size; i++){
    p[i] = rng();
  }
}

uvector lpoly_evaluator::evalfun(uvector query){
  NumericMatrix state(1, query.size(), query.begin());
  NumericMatrix state_model = build(state);
  uvector state_vector = as_ublas_vector(NumericVector(state_model.begin(), state_model.end()));
  return prod(coef_matrix, state_vector);
}
 
lpoly_jacobian::lpoly_jacobian(NumericMatrix cm, NumericMatrix trm) : terms(trm), coef_matrix(as_ublas_matrix(cm)){}

// compute jacobian
void lpoly_jacobian::operator()(uvector &q, umatrix &out, double t){
  int i,j,k;
  umatrix term_derivs = umatrix(terms.nrow(), terms.ncol());
  double x;

  for (i = 0; i < terms.nrow(); i++){
    for (j = 0; j < terms.ncol(); j++){
      // compute derivative of term i with respect to variable j
      term_derivs(i,j) = 1;
      for (k = 0; k < terms.ncol(); k++){
	term_derivs(i,j) *= pow(terms(i, k), j == k) * pow(q(k), terms(i,k) - (j == k));
      }
    }
  }

  out = prod(coef_matrix, term_derivs);
}
