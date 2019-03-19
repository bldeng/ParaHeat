#include "DistanceFile.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/isolines.h>
#include <iostream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

int main(int argc, char *argv[])
{
  if(argc != 3){
    std::cerr << "Usage: ViewScalarField MESH_FILE DATA_FILE" << std::endl;
    return 1;
  }

  // Load a mesh
  if(!igl::read_triangle_mesh(argv[1], V, F)){
    std::cerr << "Error: unable to read mesh from file " << argv[1] << std::endl;
    return 1;
  }

  // Load scalar values
  Eigen::VectorXd scalar_vals;
  if(DistanceFile::load(argv[2], scalar_vals)){
    std::cerr << "Error: unable to read scalar values from file " << argv[2] << std::endl;
  }

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  // Compute per-vertex colors
  igl::jet(scalar_vals,true,C);

  // Add per-vertex colors
  viewer.data().set_colors(C);

  Eigen::MatrixXd isoV;
  Eigen::MatrixXi isoE;
  igl::isolines(V, F, scalar_vals, 60, isoV, isoE);
  viewer.data().set_edges(isoV,isoE,Eigen::RowVector3d(0,0,0));
  viewer.data().show_lines = false;
  viewer.data().line_width = 1.5;

  // Launch the viewer
  viewer.launch();
}
