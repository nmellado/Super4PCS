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

#ifndef BBOX_H
#define BBOX_H

#include "utils/disablewarnings.h"

#include <limits>

namespace Super4PCS{

template <typename _Scalar, int _Dim>
class AABB
{
public:
    typedef _Scalar Scalar;
    enum { Dim = _Dim };

    typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

    AABB() { _min.setConstant( std::numeric_limits<Scalar>::max() / 2);
             _max.setConstant(-std::numeric_limits<Scalar>::max() / 2); }
    AABB(const VectorType& min, const VectorType& max) : _min(min), _max(max) {}
    AABB(const AABB& bb) : _min(bb._min), _max(bb._max) {}
    template <class InputIt>
    AABB(InputIt first, InputIt last) { extendTo(first, last); }

    inline AABB<Scalar, Dim>& operator=(const AABB<Scalar, Dim>& bb)
    { _min = bb._min; _max = bb._max; return (*this); }

    template <class VectorTypeDerived>
    inline void extendTo(const VectorTypeDerived& q)
    { _min = (q.array() < _min.array()).select(q, _min);
      _max = (q.array() > _max.array()).select(q, _max); }

    template <class InputIt>
    inline void extendTo(InputIt first, InputIt last)
    { std::for_each(first, last,
        std::bind1st(std::mem_fun(&AABB::extendTo), this)); }

    inline bool contains(const VectorType& q) const
    { return ((q.array() > _min.array()) && (q.array() < _max.array())).all(); }

    inline Scalar diagonal() const
    { return (_max - _min).norm(); }

    inline VectorType center() const
    { return _min + ((_max - _min) / Scalar(2.)); }

    inline const VectorType& min() const
    { return _min; }

    inline const VectorType& max() const
    { return _max; }

    inline Scalar depth() const
    { return (_max - _min)(2); }

    inline Scalar width() const
    { return (_max - _min)(1); }

    inline Scalar height() const
    { return (_max - _min)(0); }

protected:
    VectorType _min, _max;

private:

}; // class AABB

template <typename _Scalar>
using AABB3D = AABB< _Scalar, 3 >;

} // namespace Super4PCS

#endif // BBOX_H




