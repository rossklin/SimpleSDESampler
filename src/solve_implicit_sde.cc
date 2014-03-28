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

#include <Rcpp.h>

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

typedef boost::numeric::ublas::vector<double> uvector;
typedef boost::numeric::ublas::matrix<double> umatrix;

using namespace boost::numeric::odeint;
using namespace boost::math;

using namespace std;

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

struct r_deriv{
  double h;
  Rcpp::Function det;
  Rcpp::Function stoch;

  r_deriv(Rcpp::Function d, Rcpp::Function s, double delta_t) : det(d), stoch(s), h(delta_t){}

  void operator()(uvector &q, uvector &out, double t){
    uvector dpart = as_ublas_vector(det(as_r_matrix(q)));
    uvector spart = as_ublas_vector(stoch(t));
    out = dpart + (1 / sqrt(h)) * spart;
  }
};

struct r_jacobian{
  Rcpp::Function J;

  r_jacobian(Rcpp::Function j) : J(j){}

  void operator()(uvector &q, umatrix &out, double t){
    out = as_ublas_matrix(J(as_r_vector(q)));
  }
};

//' @export
// [[Rcpp::export]]
NumericMatrix solve_implicit_sde(Rcpp::Function d_det
				 , Rcpp::Function d_stoch
				 , Rcpp::Function jacobian
				 , NumericVector start
				 , double from, double to, int steps ) {

  implicit_euler<double> stepper;
  const double dt = (to - from)/steps;
  r_deriv sd = r_deriv(d_det, d_stoch, dt);
  r_jacobian sj = r_jacobian(jacobian);

  uvector state = as_ublas_vector(start);

  NumericMatrix result(steps+1, start.size());

  for(int j = 0; j < start.size(); ++j)
    result(0, j) = state[j];

  for(int i = 1; i <= steps; ++i) {
    stepper.do_step(std::make_pair(sd, sj), state, i*dt, dt);
    for(int j = 0; j < start.size(); ++j)
      result(i, j) = state[j];
  }

  return result;
}
