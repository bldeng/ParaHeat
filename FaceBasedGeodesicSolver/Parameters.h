#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>

struct Parameters
{
  Parameters()
  :heat_solver_max_iter(500), heat_solver_eps(1e-8),
   heat_solver_convergence_check_frequency(10),
   grad_solver_max_iter(500), grad_solver_eps(1e-4), penalty(50),
   grad_solver_output_frequency(50), grad_solver_convergence_check_frequency(10)
  {
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

  // Load options from file
  bool load(const char* filename);

  // Check whether the parameter values are valid
  bool valid_parameters() const;

  // Print options
  void output_options();
};



#endif /* PARAMETERS_H_ */
