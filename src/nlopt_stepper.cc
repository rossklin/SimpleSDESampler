
#include "../inst/include/SimpleSDESampler.h"

double nlopt_objective(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
  
uvector nlopt_stepper_data::evalfun(uvector state){
  return sys -> first.evalfun(state);
}

umatrix nlopt_stepper_data::jacobian(uvector state){
  umatrix out;
  sys -> second(state, out, 0);
  return out;
}

uvector nlopt_stepper_data::noisefun(double t){
  return sys -> first.noise(t);
}

nlopt_stepper::nlopt_stepper(lpoly_system_type *sys, double h, int dim, double sigma, int steps){
  data.sys = sys;
  data.dim = dim;
  data.h = h;

  data.sys -> first.build_noise(steps + 2, data.dim, sigma);
  data.sys -> first.set_h(h);

  optimizer = new nlopt::opt(nlopt::LD_TNEWTON, data.dim);
  optimizer -> set_min_objective(nlopt_objective, &data);
  //optimizer -> set_ftol_rel(0.01);
  //optimizer -> set_stopval(0.0005);
  optimizer -> set_xtol_abs(0.001);
}

nlopt_stepper::~nlopt_stepper(){
  delete optimizer;
}

void nlopt_stepper::do_step(std::vector<double> &x, double t){
  double f_opt;
  
  data.t = t;
  data.u_old = as_ublas_vector(x);

  optimizer -> optimize(x, f_opt);
}

double nlopt_objective(const std::vector<double> &x, std::vector<double> &grad, void* f_data){
  // X = x + f(X) h + sqrt(h) E
  // minimize sum_i(-X_i + f_i(X) h + x_i + sqrt(h) E_i)^2 w.r.t. X
  // g_j = sum_i 2 (-X_i + f_i(X) h + x_i + sqrt(h) E_i) (-d_ij + J_ij) 
  // g = 2 (f(X) h - X + x + sqrt(h) E) (J - I) 

  nlopt_stepper_data *s = (nlopt_stepper_data *)f_data;
  uvector u = as_ublas_vector(x);
  uvector f;
  uvector e;
  umatrix J;
  umatrix I = boost::numeric::ublas::identity_matrix<double> (s -> dim);
  uvector g;
  uvector vobj;
  double obj;
  double h = s -> h;
  uvector u_old = s -> u_old;
  double t = s -> t + h;

  int i;

  f = s -> evalfun(u);
  J = s -> jacobian(u);
  e = s -> noisefun(t);
  
  vobj = h * f + sqrt(h) * e + u_old - u;
  obj = 0;
  for(i = 0; i < s -> dim; i++){
    obj += pow(vobj(i), 2);
  }

  g = prod(2 * vobj, h * J - I);

  memcpy(&grad[0], &g(0), grad.size() * sizeof(double));

  return obj;
}
