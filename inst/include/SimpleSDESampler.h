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
};

class lpoly_jacobian{
public:
  NumericMatrix terms;
  umatrix coef_matrix;
  
  lpoly_jacobian(NumericMatrix cm, NumericMatrix trm);
  void operator()(uvector &q, umatrix &out, double t);
};

typedef std::pair<lpoly_evaluator, lpoly_jacobian> lpoly_system_type;

NumericMatrix as_r_matrix(umatrix m_);
NumericMatrix as_r_matrix(uvector v);
umatrix as_ublas_matrix(NumericMatrix m);
NumericVector as_r_vector(uvector v);
uvector as_ublas_vector(NumericVector v);
