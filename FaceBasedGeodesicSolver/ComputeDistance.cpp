#include "GeodesicDistanceSolver.h"
#include "DistanceFile.h"
#include "GetRSS.h"
#include <iostream>


int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    std::cerr << "Usage: GeodDistSolver PARAMETERS_FILE MESH_FILE DISTANCE_FILE" << std::endl;
    return 1;
  }

  Parameters param;
  if(!param.load(argv[1])){
    std::cerr << "Error: unable to load parameter file" << std::endl;
    return 1;
  }
  param.output_options();

	GeodesicDistanceSolver solver;
	if(!solver.solve(argv[2],param)){
		std::cerr << "Error in solving geodesic distance" << std::endl;
		return 1;
	}

	if(!DistanceFile::save(argv[3], solver.get_distance_values())){
		std::cerr << "Error in saving distance file" << std::endl;
		return 1;
	}

	size_t peak_mem = getPeakRSS();
	std::cout << "Peak memory usage in bytes: " << peak_mem << std::endl;

	return 0;
}




