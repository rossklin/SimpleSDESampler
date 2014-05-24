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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "../inst/include/SimpleSDESampler.h"

g_system::g_system(Rcpp::Function *fd, Rcpp::Function *fj){
  d = fd;
  j = fj;
}

uvector g_system::evalfun(uvector state){
  return as_ublas_vector((*d)(as_r_matrix(state), -1));
}

umatrix g_system::jacobian(uvector state){
  return as_ublas_matrix((*j)(as_r_matrix(state), -1));
}

uvector g_system::noisefun(double t){
  int idx = t / h + 0.5;
  int dim = noise_precache.size2();
  uvector x(dim);

  if (idx >= noise_precache.size1()){
    Rcpp::Rcout << "Invalid time index: " << idx << endl;
    exit(1);
  }

  memcpy(&x[0], &noise_precache(idx,0), dim * sizeof(double));
  return x;
}

void g_system::build_noise(int n, int dim, double sigma){
  int i;
  int size = n * dim;
  double *p;

  boost::mt19937 gener;
  boost::normal_distribution<> normal(0,sigma*sigma);
  boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal);

  boost::posix_time::time_duration diff_time = boost::posix_time::microsec_clock::local_time() - boost::posix_time::second_clock::local_time();
  unsigned int the_seed = diff_time.total_microseconds();
  
  rng.engine().seed(the_seed);
  rng.distribution().reset();

  noise_precache.resize(n, dim);
  p = &noise_precache(0,0);
    
  for (i = 0; i < size; i++){
    p[i] = rng();
  }

}

void g_system::set_noise(umatrix noise){
  noise_precache = noise;
}

void g_system::set_h(double hh){
  h = hh;
}

double g_system::get_h(){
  return h;
}
