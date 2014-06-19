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

/*
  Thanks 
  http://www.boost.org/doc/libs/1_55_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html#boost_numeric_odeint.odeint_in_detail.steppers.implicit_solvers 
  for example code on an sde
*/

template< class T > class stochastic_euler
{
public:

    typedef T state_type;
    typedef T deriv_type;
    typedef double value_type;
    typedef double time_type;
    typedef unsigned short order_type;
    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static order_type order( void ) { return 1; }

    template< class System >
    void do_step( System system , state_type &x , time_type t , time_type dt ) const
    {
        deriv_type det , stoch ;
        system.first( x , det );
        system.second( x , stoch );
        for( size_t i=0 ; i<x.size() ; ++i )
            x[i] += dt * det[i] + sqrt( dt ) * stoch[i];
    }
};

struct r_compute{
  Rcpp::Function f;

  r_compute(Rcpp::Function g) : f(g){}

  void operator()(vector<double> &q, vector<double> &out){
    out = as<vector<double> >(f(q));
  }
  
};

//' Simulates an SDE explicitly (derivative free)
//' 
//' For an Ito form SDE dx = f(x) dt + g(x) E sqrt(dt)
//' 
//' @param d_det R function representing f(x)
//' @param d_stoch R function representing g(x) E
//' @param start numeric vector with initial state
//' @param from scalar with initial time
//' @param to scalar with final time
//' @param steps number of steps to take
//' @export
// [[Rcpp::export]]
NumericMatrix solve_general_sde(Rcpp::Function d_det
			   , Rcpp::Function d_stoch
			   , vector<double> start
			   , double from, double to, int steps ) {
  stochastic_euler< vector<double> > se;
  const double dt = (to - from)/steps;
  r_compute sd = r_compute(d_det);
  r_compute ss = r_compute(d_stoch);

  vector<double> state(start);

  NumericMatrix result(steps+1, start.size());

  for(int j = 0; j < start.size(); ++j)
    result(0, j) = state[j];

  for(int i = 1; i <= steps; ++i) {
    se.do_step(std::make_pair(sd, ss), state, i*dt, dt);
    for(int j = 0; j < start.size(); ++j)
      result(i, j) = state[j];
  }

  return result;
}
