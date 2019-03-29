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


#include "Parameters.h"
#include "DistanceFile.h"
#include <iostream>

void shift_and_normalize_distance(const std::vector<int> &source_vtx,
                                  DenseVector &dist_values) {
  double mean_source_dist = 0;
  for (int i = 0; i < static_cast<int>(source_vtx.size()); ++i) {
    mean_source_dist += dist_values(source_vtx[i]);
  }
  mean_source_dist /= double(source_vtx.size());

  DenseVector source_dist;
  source_dist.setConstant(dist_values.size(), mean_source_dist);
  dist_values -= source_dist;

  dist_values /= dist_values.maxCoeff();
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr
        << "Usage: CompareDistance PARAMETERS_FILE DISTANCE_FILE REFERENCE_DISTANCE_FILE"
        << std::endl;
    return 1;
  }

  Parameters param;
  if (!param.load(argv[1])) {
    std::cerr << "Error: unable to load parameter file" << std::endl;
    return 1;
  }

  DenseVector distance, ref_distance;
  if (!DistanceFile::load(argv[2], distance)) {
    std::cerr << "Error: unable to read distance values from file " << argv[2]
              << std::endl;
    return 1;
  }

  if (!DistanceFile::load(argv[3], ref_distance)) {
    std::cerr << "Error: unable to read distance values from file " << argv[3]
              << std::endl;
    return 1;
  }

  if (distance.size() != ref_distance.size()) {
    std::cerr << "Error: different length of input distance arrays"
              << std::endl;
  }

  int n_vtx = distance.size();
  for (int i = 0; i < static_cast<int>(param.source_vertices.size()); ++i) {
    if (param.source_vertices[i] >= n_vtx) {
      std::cerr << "Error: source vertex index " << param.source_vertices[i]
                << " is out of bound" << std::endl;
      return 1;
    }
  }

  // Shift distance array such that the value at sources become zero
  shift_and_normalize_distance(param.source_vertices, distance);
  shift_and_normalize_distance(param.source_vertices, ref_distance);

  // Compute mean error
  double err = 0;
  double eps = 1e-14;
  std::vector<bool> source_flag(n_vtx, false);
  for (int i = 0; i < static_cast<int>(param.source_vertices.size()); ++i) {
    source_flag[param.source_vertices[i]] = true;
  }

  for (int i = 0; i < n_vtx; ++i) {
    if ((!source_flag[i]) && (std::fabs(ref_distance(i)) > eps)) {
      err += std::fabs(distance(i) - ref_distance(i))
          / std::fabs(ref_distance(i));
    }
  }
  err /= double(n_vtx - param.source_vertices.size());

  std::cout << "Mean relative error: " << (err * 100) << "%" << std::endl;

  return 0;
}

