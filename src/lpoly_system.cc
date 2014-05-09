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

lpoly_system::lpoly_system(lpoly_system_type *s){
  sys = s;
}

uvector lpoly_system::evalfun(uvector state){
  return sys -> first.evalfun(state);
}

umatrix lpoly_system::jacobian(uvector state){
  umatrix out;
  sys -> second(state, out, 0);
  return out;
}

uvector lpoly_system::noisefun(double t){
  return sys -> first.noise(t);
}

void lpoly_system::build_noise(int n, int dim, double sigma){
  sys -> first.build_noise(n, dim, sigma);
}

void lpoly_system::set_noise(umatrix noise){
  sys -> first.set_noise(noise);
}

void lpoly_system::set_h(double hh){
  sys -> first.set_h(hh);
}

double lpoly_system::get_h(){
  return sys -> first.h;
}
