// Copyright 2014 Nicolas Mellado
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
// Authors: Moos Hueting, Nicolas Mellado
//
// This file is part of the implementation of the Super 4-points Congruent Sets
// (Super 4PCS) algorithm presented in:
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

#ifndef _SUPER4PCS_ACCELERATORS_BBOX_H
#define _SUPER4PCS_ACCELERATORS_BBOX_H

#include "super4pcs/utils/disablewarnings.h"
#include <Eigen/Geometry>
#include <limits>

namespace GlobalRegistration{

template <typename _Scalar, int _Dim>
class AABB : public Eigen::AlignedBox<_Scalar, _Dim>
{
public:
    using Base = Eigen::AlignedBox<_Scalar, _Dim>;
    using Scalar = typename Base::Scalar;
    static constexpr int Dim = _Dim;

    using VectorType = Eigen::Matrix<Scalar, Dim, 1>;

    using Base::extend;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AABB(Scalar min =  std::numeric_limits<Scalar>::max() / 2,
         Scalar max = -std::numeric_limits<Scalar>::max() / 2)
      : Base(VectorType::Constant(min), VectorType::Constant(max)) {}

    template <class InputIt>
    inline void extend(InputIt first, InputIt last)
    { std::for_each(first, last,
        std::bind1st(std::mem_fun(&Base::extend), this)); }

}; // class AABB

template <typename _Scalar>
using AABB3D = AABB< _Scalar, 3 >;

} // namespace Super4PCS

#endif // BBOX_H




