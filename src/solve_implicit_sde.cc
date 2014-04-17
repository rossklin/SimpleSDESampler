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

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

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
  umatrix noise_precache;

  r_deriv(Rcpp::Function d, Rcpp::Function s, double delta_t) : det(d), stoch(s), h(delta_t){}

  void operator()(uvector &q, uvector &out, double t){
    NumericMatrix rq = as_r_matrix(q);
    uvector dpart = as_ublas_vector(det(rq, t));
    umatrix spart = as_ublas_matrix(stoch(rq, t));
    uvector buf(spart.size1());
    
    axpy_prod(spart, noise(t), buf);

    out = dpart + (1 / sqrt(h)) * buf;
  }

  uvector noise(double t){
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

  void build_noise(int n, int d, double sigma){
    int i;
    int size = n * d;
    double *p;

    boost::mt19937 gener;
    boost::normal_distribution<> normal(0,sigma);
    boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal);

    rng.engine().seed(clock());
    rng.distribution().reset();

    noise_precache.resize(n, d);
    p = &noise_precache(0,0);
    
    for (i = 0; i < size; i++){
      p[i] = rng();
    }
  }
};

struct r_jacobian{
  Rcpp::Function J;

  r_jacobian(Rcpp::Function j) : J(j){}

  void operator()(uvector &q, umatrix &out, double t){
    out = as_ublas_matrix(J(as_r_vector(q), t));
  }
};

//' @param d_det Deterministic component: an R function: (m x n matrix of m states, scalar time) -> m x n matrix of m time derivatives
//' @param d_stoch Stochastic component: an R function: (1 x n matrix state, scalar time) -> 1 x n matrix state
//' @param jacobian Jacobian of deterministic component: an R function: (n vector state, scalar time) -> n x n matrix df_i / du_j
//' @param sigma Amplitude of noise: scalar
//' @param start Initial position: n vector
//' @param from Initial time: scalar
//' @param to Final time: scalar
//' @param steps Number of points to take, s.t. dt = (from - to) / (steps + 1): integer
//' @export
// [[Rcpp::export]]
NumericMatrix solve_implicit_sde(Rcpp::Function d_det
				 , Rcpp::Function d_stoch
				 , Rcpp::Function jacobian
				 , double sigma
				 , NumericVector start
				 , double from, double to, int steps ) {

  implicit_euler<double> stepper;
  const double dt = (to - from)/steps;
  r_deriv sd = r_deriv(d_det, d_stoch, dt);
  r_jacobian sj = r_jacobian(jacobian);
  uvector state = as_ublas_vector(start);
  NumericMatrix result(steps+1, start.size());

  sd.build_noise(steps + 2, start.size(), sigma);

  for(int j = 0; j < start.size(); ++j)
    result(0, j) = state[j];

  for(int i = 1; i <= steps; ++i) {
    stepper.do_step(std::make_pair(sd, sj), state, i*dt, dt);
    for(int j = 0; j < start.size(); ++j)
      result(i, j) = state[j];
  }

  return result;
}

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
					   , Rcpp::Function d_stoch
					   , Rcpp::Function jacobian
					   , double sigma
					   , NumericVector start
					   , double from, double to, int steps){
  
  NumericMatrix result(steps+1, start.size());
  NumericMatrix buf;
  int i,j;
  double *p, *q;
  int smax = result.nrow() * result.ncol();

  p = &result(0,0);

  memset(p, 0, smax * sizeof(double));
  
  for (i = 0; i < nrep; i++){
    buf = solve_implicit_sde( d_det
			      , d_stoch
			      , jacobian
			      , sigma
			      , start
			      , from
			      , to
			      , steps);

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
