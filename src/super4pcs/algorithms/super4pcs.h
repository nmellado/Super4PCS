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

#ifndef _SUPER4PCS_ALGO_SUPER4PCS_H_
#define _SUPER4PCS_ALGO_SUPER4PCS_H_

#include "super4pcs/algorithms/match4pcsBase.h"
#include "super4pcs/algorithms/pairCreationFunctor.h"

namespace GlobalRegistration {

/// Class for the computation of the 4PCS algorithm.
class MatchSuper4PCS : public Match4PCSBase {
 public:
  using Base        = Match4PCSBase;
  using Scalar      = typename Base::Scalar;
  using PairsVector = typename Base::PairsVector;

  explicit MatchSuper4PCS(const Match4PCSOptions& options,
                          const Utils::Logger &logger);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ~MatchSuper4PCS();

private:
 /// Private data contains parameters and internal variables that are computed
 /// and change during the match computation. All parameters have default
 /// values.

 /// Internal data members.

 mutable PairCreationFunctor<Scalar> pcfunctor_;

protected:
  /// Constructs pairs of points in Q, corresponding to a single pair in the
  /// in basein P.
  /// @param [in] pair_distance The distance between the pairs in P that we have
  /// to match in the pairs we select from Q.
  /// @param [in] pair_normal_distance The angle between the normals of the pair
  /// in P.
  /// @param [in] pair_distance_epsilon Tolerance on the pair distance. We allow
  /// candidate pair in Q to have distance of
  /// pair_distance+-pair_distance_epsilon.
  /// @param [in] base_point1 The index of the first point in P.
  /// @param [in] base_point2 The index of the second point in P.
  /// @param [out] pairs A set of pairs in Q that match the pair in P with
  /// respect to distance and normals, up to the given tolerance.
  void
  ExtractPairs(
          Scalar pair_distance,
          Scalar pair_normals_angle,
          Scalar pair_distance_epsilon,
          int base_point1,
          int base_point2,
          PairsVector* pairs) const override;

  /// Finds congruent candidates in the set Q, given the invariants and threshold
  /// distances. Returns true if a non empty set can be found, false otherwise.
  /// @param invariant1 [in] The first invariant corresponding to the set P_pairs
  /// of pairs, previously extracted from Q.
  /// @param invariant2 [in] The second invariant corresponding to the set
  /// Q_pairs of pairs, previously extracted from Q.
  /// @param [in] distance_threshold1 The distance for verification.
  /// @param [in] distance_threshold2 The distance for matching middle points due
  /// to the invariants (See the paper for e1, e2).
  /// @param [in] P_pairs The first set of pairs.
  /// @param [in] Q_pairs The second set of pairs.
  /// @param [out] quadrilaterals The set of congruent quadrilateral. In fact,
  /// it's a super set from which we extract the real congruent set.
  bool FindCongruentQuadrilaterals(
          Scalar invariant1,
          Scalar invariant2,
          Scalar distance_threshold1,
          Scalar distance_threshold2,
          const PairsVector& P_pairs,
          const PairsVector& Q_pairs,
          std::vector<Quadrilateral>* quadrilaterals) const override;

  /// Initializes the data structures and needed values before the match
  /// computation.
  /// @param [in] point_P First input set.
  /// @param [in] point_Q Second input set.
  /// expected to be in the inliers.
  void Initialize(const std::vector<Point3D>& P,
                 const std::vector<Point3D>& Q) override;
};
} /// namespace Super4PCS

#endif  /// _4PCS_H_
