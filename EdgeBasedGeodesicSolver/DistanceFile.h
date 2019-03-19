#ifndef DISTANCEFILE_H
#define DISTANCEFILE_H

#include "EigenTypes.h"
#include <fstream>
#include <iostream>

class DistanceFile
{

 public:

    static bool save(const char *file_name, const DenseVector &dist_values)
    {
      int n_values = dist_values.size();
      if(n_values <= 0){
        std::cerr << "Error: empty distance value vector" << std::endl;
        return false;
      }

      std::ofstream ofile(file_name);
      if (!ofile.is_open()){
        std::cerr << "Unable to open file " << file_name << std::endl;
        return false;
      }

      ofile << n_values << std::endl;
      if(!ofile){
        std::cerr << "Error writing to file " << file_name << std::endl;
        return false;
      }

      for(int i = 0; i < n_values; ++ i){
        ofile << dist_values(i) << std::endl;
        if(!ofile){
          std::cerr << "Error writing to file " << file_name << std::endl;
          return false;
        }
      }

        return true;
    }

    static bool load(const char *file_name, DenseVector &dist_values)
    {
    std::ifstream ifile(file_name);
    if (!ifile.is_open()){
      std::cerr << "Unable to open file " << file_name << std::endl;
      return false;
    }

    int n_values = 0;
    if(!(ifile >> n_values)){
      std::cerr << "Error parsing the number of values" << std::endl;
      return false;
    }

    if(n_values <= 0){
      std::cerr << "Error: invalid number or values" << std::endl;
      return false;
    }

    dist_values.setZero(n_values);

    double val = 0;
    for(int i = 0; i < n_values; ++ i){
      if(ifile >> val){
        dist_values(i) = val;
      }
      else{
        std::cerr << "Error parsing distance value at position " << i << std::endl;
        return false;
      }
    }

        return true;
    }
};

#endif // DISTANCEFILE_H
