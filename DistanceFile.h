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


#ifndef DISTANCEFILE_H_
#define DISTANCEFILE_H_

#include "EigenTypes.h"
#include <fstream>
#include <iostream>

class DistanceFile {

 public:

  static bool save(const char *file_name, const DenseVector &dist_values) {
    int n_values = dist_values.size();
    if (n_values <= 0) {
      std::cerr << "Error: empty distance value vector" << std::endl;
      return false;
    }

    std::ofstream ofile(file_name);
    if (!ofile.is_open()) {
      std::cerr << "Unable to open file " << file_name << std::endl;
      return false;
    }

    ofile << n_values << std::endl;
    if (!ofile) {
      std::cerr << "Error writing to file " << file_name << std::endl;
      return false;
    }

    for (int i = 0; i < n_values; ++i) {
      ofile << dist_values(i) << std::endl;
      if (!ofile) {
        std::cerr << "Error writing to file " << file_name << std::endl;
        return false;
      }
    }

    return true;
  }

  static bool load(const char *file_name, DenseVector &dist_values) {
    std::ifstream ifile(file_name);
    if (!ifile.is_open()) {
      std::cerr << "Unable to open file " << file_name << std::endl;
      return false;
    }

    int n_values = 0;
    if (!(ifile >> n_values)) {
      std::cerr << "Error parsing the number of values" << std::endl;
      return false;
    }

    if (n_values <= 0) {
      std::cerr << "Error: invalid number or values" << std::endl;
      return false;
    }

    dist_values.setZero(n_values);

    double val = 0;
    for (int i = 0; i < n_values; ++i) {
      if (ifile >> val) {
        dist_values(i) = val;
      } else {
        std::cerr << "Error parsing distance value at position " << i
                  << std::endl;
        return false;
      }
    }

    return true;
  }
};

#endif /* DISTANCEFILE_H_ */
