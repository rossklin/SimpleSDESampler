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

#include <vector>
#include <utility>
#include <exception>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "nlopt.hpp"

#include <Rcpp.h>

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::CharacterVector;
using Rcpp::XPtr;
using Rcpp::wrap;
using Rcpp::as;

using namespace boost::numeric::odeint;
using namespace boost::math;
using namespace std;

typedef boost::numeric::ublas::vector<double> uvector;
typedef boost::numeric::ublas::matrix<double> umatrix;

class lpoly_evaluator{
public:
  NumericMatrix terms;
  umatrix coef_matrix;
  umatrix noise_precache;
  double h;
  
  lpoly_evaluator(NumericMatrix cm, NumericMatrix trm);
  void set_h(double v);
  NumericMatrix build(NumericMatrix data);
  void operator()(uvector &query, uvector &out, double t);
  uvector noise(double t);
  void build_noise(int n, int d, double sigma);
  void set_noise(umatrix m);
  uvector evalfun(uvector q);
};

class lpoly_jacobian{
public:
  NumericMatrix terms;
  umatrix coef_matrix;
  
  lpoly_jacobian(NumericMatrix cm, NumericMatrix trm);
  void operator()(uvector &q, umatrix &out, double t);
};

typedef std::pair<lpoly_evaluator, lpoly_jacobian> lpoly_system_type;

class general_system{
 public:
  int dim;
  uvector u_old;
  double t;
  
  virtual uvector evalfun(uvector state) = 0;
  virtual umatrix jacobian(uvector state) = 0;
  virtual uvector noisefun(double t) = 0;
  virtual void build_noise(int n, int dim, double sigma) = 0;
  virtual void set_noise(umatrix noise) = 0;
  virtual void set_h(double h) = 0;
  virtual double get_h() = 0;
};

class lpoly_system : public general_system{
 public:
  lpoly_system_type *sys;
  
  lpoly_system(lpoly_system_type *s);
  uvector evalfun(uvector state);
  umatrix jacobian(uvector state);
  uvector noisefun(double t);
  void build_noise(int n, int dim, double sigma);
  void set_noise(umatrix noise);
  void set_h(double h);
  double get_h();
};

class g_system : public general_system{
 public:
  double h;
  Rcpp::Function *d;
  Rcpp::Function *j;
  umatrix noise_precache;

  g_system(Rcpp::Function *fd, Rcpp::Function *fj);
  uvector evalfun(uvector state);
  umatrix jacobian(uvector state);
  uvector noisefun(double t);
  void build_noise(int n, int dim, double sigma);
  void set_noise(umatrix noise);
  void set_h(double h);
  double get_h();
};

struct nlopt_stepper{
  general_system *sys;
  nlopt::opt *optimizer;

  nlopt_stepper(general_system *sys, double h, int dim, double sigma, int steps, double xtol = 0, const char* algorithm = "TNEWTON");
  nlopt_stepper(general_system *sys, double h, int dim, umatrix noise, double xtol = 0, const char* algorithm = "TNEWTON");
  ~nlopt_stepper();
  void setup(general_system *sys, double h, int dim, double xtol, const char* algorithm);
  void do_step(std::vector<double> &x, double t);
};

NumericMatrix lpoly_sde(nlopt_stepper &stepper
			, NumericVector start
			, double from, double to, int steps 
			, double x_tol = 0
			, const char* algorithm = "TNEWTON");


NumericMatrix as_r_matrix(umatrix m_);
NumericMatrix as_r_matrix(uvector v);
umatrix as_ublas_matrix(NumericMatrix m);
NumericVector as_r_vector(uvector v);
uvector as_ublas_vector(NumericVector v);
uvector as_ublas_vector(std::vector<double> v);
