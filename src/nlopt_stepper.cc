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

double nlopt_objective(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
  
nlopt_stepper::nlopt_stepper(general_system *s, double h, int dim, double sigma, int steps, double xtol, const char* algorithm){
  setup(s, h, dim, xtol, algorithm);
  s -> build_noise(steps + 2, dim, sigma);
}

nlopt_stepper::nlopt_stepper(general_system *s, double h, int dim, umatrix noise, double xtol, const char* algorithm){
  setup(s, h, dim, xtol, algorithm);
  s -> set_noise(noise);
}

nlopt_stepper::~nlopt_stepper(){
  delete optimizer;
}

void nlopt_stepper::setup(general_system *s, double h, int dim, double xtol, const char* algorithm){

  sys = s;
  sys -> dim = dim;
  sys -> set_h(h);

  if (!strcmp(algorithm, "TNEWTON")){
    optimizer = new nlopt::opt(nlopt::LD_TNEWTON, sys -> dim);
  }else if (!strcmp(algorithm, "TNEWTON_RESTART")){
    optimizer = new nlopt::opt(nlopt::LD_TNEWTON_RESTART, sys -> dim);
  }else if (!strcmp(algorithm, "LBFGS")){
    optimizer = new nlopt::opt(nlopt::LD_LBFGS, sys -> dim);
  }else{
    Rcpp::Rcout << "nlopt_stepper: invalid algorithm: " << algorithm << endl;
    exit(-1);
  }

  optimizer -> set_min_objective(nlopt_objective, sys);
  //optimizer -> set_ftol_abs(xtol);
  //optimizer -> set_ftol_rel(0.01);
  optimizer -> set_stopval(xtol);
  //optimizer -> set_xtol_abs(xtol);
  
}

void nlopt_stepper::do_step(std::vector<double> &x, double t){
  double f_opt;

#ifdef VERBOSE
  Rcpp::Rcout << "STARTING STEP FROM " << x[0] << " at t = " << t << endl;
#endif
  
  sys -> t = t;
  sys -> u_old = as_ublas_vector(x);

  optimizer -> optimize(x, f_opt);

#ifdef VERBOSE
  Rcpp::Rcout << "RESULTING POSTION: " << x[0] << ", obj = " << f_opt << endl;
#endif
}

double nlopt_objective(const std::vector<double> &x, std::vector<double> &grad, void* f_data){
  // X = x + f(X) h + sqrt(h) E
  // minimize sum_i(-X_i + f_i(X) h + x_i + sqrt(h) E_i)^2 w.r.t. X
  // g_j = sum_i 2 (-X_i + f_i(X) h + x_i + sqrt(h) E_i) (-d_ij + J_ij) 
  // g = 2 (f(X) h - X + x + sqrt(h) E) (J - I) 

  general_system *s = (general_system*)f_data;
  uvector u = as_ublas_vector(x);
  uvector f;
  uvector e;
  umatrix J;
  umatrix I = boost::numeric::ublas::identity_matrix<double> (s -> dim);
  uvector g;
  uvector vobj;
  double obj;
  double h = s -> get_h();
  uvector u_old = s -> u_old;
  double t = s -> t;

  int i;

  f = s -> evalfun(u);
  J = s -> jacobian(u);
  e = s -> noisefun(t);

  vobj = h * f + sqrt(h) * e + u_old - u;
  obj = 0;
  for(i = 0; i < s -> dim; i++){
    obj += pow(vobj(i), 2);
  }

#ifdef VERBOSE
  Rcpp::Rcout << "t " << t << ": from " << u_old(0) << " to " << u(0) << ": sqrt(h) e = " << sqrt(h) * e(0) << "(du = " << u(0) - u_old(0) <<"), evalfun = " << f(0) << ", J = " << J(0,0) << endl;
#endif

  g = prod(2 * vobj, h * J - I);

  memcpy(&grad[0], &g(0), grad.size() * sizeof(double));

  return obj;
}
