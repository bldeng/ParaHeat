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



#ifndef FACEBASEDGEODESICSOLVER_H_
#define fACEBASEDGEODESICSOLVER_H_

#include "EigenTypes.h"
#include "surface_mesh/Surface_mesh.h"
#include "Parameters.h"
#include <fstream>

class FaceBasedGeodesicSolver {
 public:
  FaceBasedGeodesicSolver();

  bool solve(const char* mesh_file, const Parameters &para);

  const DenseVector& get_distance_values();

  const DenseVector& get_heat_solution();

 private:

  typedef surface_mesh::Surface_mesh MeshType;
  MeshType mesh;
  double model_scaling_factor;

  Parameters param;

  IndexVector bfs_vertex_list;  // Non-source vertex indices stored according to their BFS order
  IndexVector bfs_segment_addr;  // Addresses within bfs_vertex_list for the first vertex in each BFS layer

  std::pair<int, double>* bfs_laplacian_coef;  // Vertices and their weights for evaluating cotan Laplacian
  IndexVector bfs_laplacian_coef_addr;  // Starting addresses for the segment of each Laplacian within bfs_laplacian_vtx

  IndexVector transition_halfedge_idx;
  IndexVector transition_from_vtx;
  Matrix3X transition_edge_vector;  // For each vertex in bfs_vertex_list, store the vector of a halfedge pointing to that vertex and representing the transition direction for recovering distance values from gradients
  Matrix2Xi transition_edge_neighbor_faces;  // Neighboring face indices for each transition edge

  Matrix2Xi S;  // Paper : S  selection matrix, each column storing the two face indices associated with an internal edge
  Matrix3X edge_vector;   // paper : e   edges unit vector
  Matrix3X init_grad;   // initial gradients computed from heat flow

  Matrix3X G;   // Paper : G   gradients for each face
  Matrix3X Y;  // paper : Y   auxiliary variable for the compatibility condition (Y = S * G)
  Matrix3X D;  // Paper : scaled dual variables lambda / (mu * sqrt(area));

  Matrix3X e;  // Unit vectors of internal edges
  DenseVector face_area;
  DenseVector Y_area;   // Face area associated with each column of Y
  DenseVector Y_area_squared;  // squared values of Y_area, used for computing dual residual squared norm
  bool need_compute_residual_norms;

  Matrix3X SG1, SG2;  // Storage for S * G
  Matrix3X *prev_SG, *current_SG;  // Pointer to current and previous S*G buffers

  Matrix3Xi faces_Y_index;  // the set of columns in matrix Y that corresponding to each face

  DenseVector geod_dist_values;

  int n_vertices;         // number of vertices
  int n_faces;            // number of faces
  int n_edges;            // number of edges
  int n_interior_edges;            // number of interior edges

  int iter_num;

  // Variables for primal and dual residuals
  double primal_residual_sqr_norm, dual_residual_sqr_norm;
  double primal_residual_sqr_norm_threshold, dual_residual_sqr_norm_threshold;

  // Variables for the progress of the solver
  bool output_progress;
  bool optimization_converge, optimization_end;

  bool load_input(const char* mesh_file);
  void normalize_mesh();

  void init_bfs_paths();
  void prepare_integrate_geodesic_distance();
  void gauss_seidel_init_gradients();
  void compute_integrable_gradients();
  void integrate_geodesic_distance();

  void update_Y();                          // update Y
  void update_G();                          // update G
  void update_dual_variables();  // update dual variables and check if the solver converges

  typedef long double HeatScalar;
  typedef Eigen::Matrix<HeatScalar, Eigen::Dynamic, 1> VectorHS;
  typedef Eigen::Matrix<HeatScalar, 3, 1> Vector3HS;
  typedef Eigen::Matrix<HeatScalar, 3, 3> Matrix3HS;
  void compute_heatflow_residual(const VectorHS &heat_values,
                                 HeatScalar init_source_val,
                                 VectorHS &residuals);

};

#endif /* GEODESICDISTANCESOLVER_H_ */
