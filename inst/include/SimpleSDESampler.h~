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

#include <boost/optional/optional.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <Rcpp.h>

using std::vector;

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::CharacterVector;
using Rcpp::XPtr;
using Rcpp::wrap;
using Rcpp::as;

using namespace boost::numeric::odeint;
using namespace boost::math;
using namespace std;
