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

#ifndef CGAL_POLYGON_MESH_PROCESSING_SUPER4PCS_ALIGN_H
#define CGAL_POLYGON_MESH_PROCESSING_SUPER4PCS_ALIGN_H

#include <CGAL/Aff_transformation_3.h>

namespace CGAL {
namespace Polygon_mesh_processing {

  /*!
   * \ingroup PMP_meshing_grp
   *
   * \brief Align two 3D polygon meshes using the Super4PCS algorithm
   *
   * @param P first polygon mesh
   * @param Q second polygon mesh
   * @param np sequence of \ref namedparameters among the ones listed below
   *
   * @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
   * that has an internal property map for `boost::vertex_point_t`
   * @tparam NamedParameters a sequence of \ref namedparameters
   *
   * \cgalNamedParamsBegin
   *    \cgalParamBegin{delta} Scalar value used to configuration the LCP computation\cgalParamEnd
   *    \cgalParamBegin{overlap} Scalar \in[0:1] defining the amount of overlap between the two clouds\cgalParamEnd
   *    \cgalParamBegin{sample_size} Number of samples used for registration\cgalParamEnd
   * \cgalNamedParamsEnd
   *
   * \return The transformation aligning the second input mesh to the first
   */
  template <class PolygonMesh, typename NamedParameters>
  static inline
  Aff_transformation_3 align(const PolygonMesh     & P,
                             const PolygonMesh     & Q,
                             const NamedParameters & np );

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif
