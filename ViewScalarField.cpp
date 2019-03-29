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

#include "DistanceFile.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/isolines.h>
#include <iostream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: ViewScalarField MESH_FILE DATA_FILE" << std::endl;
    return 1;
  }

  // Load a mesh
  if (!igl::read_triangle_mesh(argv[1], V, F)) {
    std::cerr << "Error: unable to read mesh from file " << argv[1]
              << std::endl;
    return 1;
  }

  // Load scalar values
  Eigen::VectorXd scalar_vals;
  if (!DistanceFile::load(argv[2], scalar_vals)) {
    std::cerr << "Error: unable to read scalar values from file " << argv[2]
              << std::endl;
  }

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  // Compute per-vertex colors
  igl::jet(scalar_vals, true, C);

  // Add per-vertex colors
  viewer.data().set_colors(C);

  Eigen::MatrixXd isoV;
  Eigen::MatrixXi isoE;
  igl::isolines(V, F, scalar_vals, 60, isoV, isoE);
  viewer.data().set_edges(isoV, isoE, Eigen::RowVector3d(0, 0, 0));
  viewer.data().show_lines = false;
  viewer.data().line_width = 1.5;

  // Launch the viewer
  viewer.launch();
}
