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
// Authors: Nicolas Mellado
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

#ifndef CGAL_POINT_SET_SUPER4PCS_ALIGN_H
#define CGAL_POINT_SET_SUPER4PCS_ALIGN_H

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

namespace CGAL {

  /*!
   * \ingroup PkgPointSetProcessing
   *
   * \brief Align two 3d point-clouds using the Super4PCS algorithm
   *
   * @param b1 iterator over the first input point of the first point clouds.
   * @param e1 past-the-end iterator over the input points of the first point clouds.
   * @param b2 iterator over the first input point of the second point clouds.
   * @param e2 past-the-end iterator over the input points of the second point clouds.
   * @param point_pmap property map: value_type of InputIterator -> Point_3 or Point_3
   * @param np sequence of \ref namedparameters among the ones listed below
   *
   * @tparam ForwardIterator iterator over input points.
   * @tparam PointPmap is a model of `ReadablePropertyMap` with value type
   * `Point_3<Kernel>`. It can be omitted if the value type of `ForwardIterator`
   * is convertible to `Point_3<Kernel>`.
   * @tparam NamedParameters a sequence of \ref namedparameters
   *
   * \cgalNamedParamsBegin
   *    \cgalParamBegin{delta} Scalar value used to configuration the LCP computation\cgalParamEnd
   *    \cgalParamBegin{overlap} Scalar \in[0:1] defining the amount of overlap between the two clouds\cgalParamEnd
   *    \cgalParamBegin{sample_size} Number of samples used for registration\cgalParamEnd
   * \cgalNamedParamsEnd
   *
   * \return The Transformation aligning the second input point-cloud to the first
   */
  template <typename Kernel,
            class ForwardIterator,
            class PointPmap,
            typename NamedParameters>
  Aff_transformation_3<Kernel> align(
      ForwardIterator b1,
      ForwardIterator e1,
      ForwardIterator b2,
      ForwardIterator e2,
      PointPmap point_pmap,
      const NamedParameters & np){
    Aff_transformation_3<Kernel> tr;

    return tr;
  }

} // namespace CGAL

#endif
