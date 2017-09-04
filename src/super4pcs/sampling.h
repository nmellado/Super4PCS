// Copyright 2017 Nicolas Mellado
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// -------------------------------------------------------------------------- //
//
// Authors: Nicolas Mellado, Dror Aiger
//
// An implementation of the Super 4-points Congruent Sets (Super 4PCS)
// algorithm presented in:
//
// Super 4PCS: Fast Global Pointcloud Registration via Smart Indexing
// Nicolas Mellado, Dror Aiger, Niloy J. Mitra
// Symposium on Geometry Processing 2014.
//
// Data acquisition in large-scale scenes regularly involves accumulating
// information across multiple scans. A common approach is to locally align scan
// pairs using Iterative Closest Point (ICP) algorithm (or its variants), but
// requires static scenes and small motion between scan pairs. This prevents
// accumulating data across multiple scan sessions and/or different acquisition
// modalities (e.g., stereo, depth scans). Alternatively, one can use a global
// registration algorithm allowing scans to be in arbitrary initial poses. The
// state-of-the-art global registration algorithm, 4PCS, however has a quadratic
// time complexity in the number of data points. This vastly limits its
// applicability to acquisition of large environments. We present Super 4PCS for
// global pointcloud registration that is optimal, i.e., runs in linear time (in
// the number of data points) and is also output sensitive in the complexity of
// the alignment problem based on the (unknown) overlap across scan pairs.
// Technically, we map the algorithm as an 'instance problem' and solve it
// efficiently using a smart indexing data organization. The algorithm is
// simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in
// significant speedup over alternative approaches and allows unstructured
// efficient acquisition of scenes at scales previously not possible. Complete
// source code and datasets are available for research use at
// http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.

#ifndef _PCS_UTILS_H_
#define _PCS_UTILS_H_

#include <vector>
#include <array>
#include "super4pcs/shared4pcs.h"


namespace GlobalRegistration {
namespace Sampling {


struct UniformDistSampler{
private:
    template <typename _Scalar>
    class HashTable {
    public:
        using Scalar = _Scalar;

    private:
        const uint64_t MAGIC1 = 100000007;
        const uint64_t MAGIC2 = 161803409;
        const uint64_t MAGIC3 = 423606823;
        const uint64_t NO_DATA = 0xffffffffu;
        Scalar voxel_;
        Scalar scale_;
        using VoxelType = std::array<int,3>;
        std::vector<VoxelType> voxels_;
        std::vector<uint64_t> data_;

    public:
        HashTable(int maxpoints, Scalar voxel) : voxel_(voxel), scale_(1.0f / voxel) {
            uint64_t n = maxpoints;
            voxels_.resize(n);
            data_.resize(n, NO_DATA);
        }
        template <typename Point>
        uint64_t& operator[](const Point& p) {
            // TODO: use eigen power here.
            VoxelType c {int(floor(p.x() * scale_)),
                         int(floor(p.y() * scale_)),
                         int(floor(p.z() * scale_))};

            uint64_t key = (MAGIC1 * c[0] + MAGIC2 * c[1] + MAGIC3 * c[2]) % data_.size();
            while (1) {
                if (data_[key] == NO_DATA) {
                    voxels_[key] = c;
                    break;
                } else if (voxels_[key] == c) {
                    break;
                }
                key++;
                if (key == data_.size()) key = 0;
            }
            return data_[key];
        }
    };
public:
    template <typename Point>
    inline
    void operator() (const std::vector<Point>& inputset,
                     const Match4PCSOptions& options,
                     std::vector<Point>& output) const {
      int num_input = inputset.size();
      output.clear();
      HashTable<typename Point::Scalar> hash(num_input, options.delta);
      for (int i = 0; i < num_input; i++) {
        uint64_t& ind = hash[inputset[i]];
        if (ind >= num_input) {
          output.push_back(inputset[i]);
          ind = output.size();
        }
      }
    }
};


} // namespace Sampling
} // namespace Super4PCS


#endif
