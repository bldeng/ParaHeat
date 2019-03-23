// BSD 3-Clause License
//
// Copyright (c) 2019, Jiong Tao, Bailin Deng, Yue Peng
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>

struct Parameters {
  Parameters()
      : heat_solver_max_iter(500),
        heat_solver_eps(1e-8),
        heat_solver_convergence_check_frequency(10),
        grad_solver_max_iter(500),
        grad_solver_eps(1e-4),
        penalty(50),
        grad_solver_output_frequency(50),
        grad_solver_convergence_check_frequency(10),
        solver_type(0) {
    source_vertices.push_back(0);
  }

  // Parameters for heat solver
  int heat_solver_max_iter;
  double heat_solver_eps;
  int heat_solver_convergence_check_frequency;

  // Parameters for the ADMM solver for gradient correction
  int grad_solver_max_iter;
  double grad_solver_eps;
  double penalty;
  int grad_solver_output_frequency;
  int grad_solver_convergence_check_frequency;

  // Indices for source vertices
  std::vector<int> source_vertices;

  // SolverType. 0 for face based algorithm; 1 for edge based algorithm
  int solver_type;

  // Load options from file
  bool load(const char* filename);

  // Check whether the parameter values are valid
  bool valid_parameters() const;

  // Print options
  void output_options();
};

#endif /* PARAMETERS_H_ */
