
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
