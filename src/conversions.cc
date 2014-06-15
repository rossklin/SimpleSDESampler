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

#include "../inst/include/SimpleSDESampler.h"

NumericMatrix as_r_matrix(umatrix m_){
  // as odeint::implicit_euler insists on using row-major ublas::matrix, we must transpose here...
  umatrix m = boost::numeric::ublas::trans(m_);
  NumericMatrix x(m.size2(), m.size1());
  memcpy(&x[0], &m(0,0), m.size1() * m.size2() * sizeof(double));
  return x;
}

NumericMatrix as_r_matrix(uvector v){
  NumericMatrix x(1, v.size());
  memcpy(&x[0], &v[0], v.size() * sizeof(double));
  return x;
}

umatrix as_ublas_matrix(NumericMatrix m){
  // as odeint::implicit_euler insists on using row-major ublas::matrix, we must transpose here...
  umatrix x_(m.ncol(), m.nrow());
  memcpy(&x_(0,0), &m[0], m.nrow() * m.ncol() * sizeof(double));
  return boost::numeric::ublas::trans(x_);
}

NumericVector as_r_vector(uvector v){
  NumericVector x(v.size());
  memcpy(&x[0], &v[0], v.size() * sizeof(double));
  return x;
}

uvector as_ublas_vector(NumericVector v){
  uvector x(v.size());
  memcpy(&x[0], &v[0], v.size() * sizeof(double));
  return x;
}

uvector as_ublas_vector(std::vector<double> v){
  uvector x(v.size());
  memcpy(&x[0], &v[0], v.size() * sizeof(double));
  return x;
}
