// Copyright 2012 Dror Aiger
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
// Authors: Dror Aiger, Yoni Weill
//
// An implementation of the 4-points Congruent Sets (4PCS) algorithm presented
// in:
//
// 4-points Congruent Sets for Robust Surface Registration
// Dror Aiger, Niloy J. Mitra, Daniel Cohen-Or
// ACM SIGGRAPH 2008 and ACM Transaction of Graphics.
//
// Given two sets of points in 3-space, P and Q, the algorithm applies RANSAC
// in roughly O(n^2) time instead of O(n^3) for standard RANSAC, using an
// efficient method based on invariants, to find the set of all 4-points in Q
// that can be matched by rigid transformation to a given set of 4-points in P
// called a base. This avoids the need to examine all sets of 3-points in Q
// against any base of 3-points in P as in standard RANSAC.
// The algorithm can use colors and normals to speed-up the matching
// and to improve the quality. It can be easily extended to affine/similarity
// transformation but then the speed-up is smaller because of the large number
// of congruent sets. The algorithm can also limit the range of transformations
// when the application knows something on the initial pose but this is not
// necessary in general (though can speed the runtime significantly).

// Home page of the 4PCS project (containing the paper, presentations and a
// demo): http://graphics.stanford.edu/~niloy/research/fpcs/fpcs_sig_08.html
// Use google search on "4-points congruent sets" to see many related papers
// and applications.

#ifndef _SUPER4PCS_ALGO_4PCS_H_
#define _SUPER4PCS_ALGO_4PCS_H_

#include "super4pcs/algorithms/match4pcsBase.h"

namespace GlobalRegistration {

/// Class for the computation of the 4PCS algorithm.
class Match4PCS : public Match4PCSBase {
public:
    using Base        = Match4PCSBase;
    using Scalar      = typename Base::Scalar;
    using PairsVector = typename Base::PairsVector;
    using VectorType  = typename Base::VectorType;

    explicit Match4PCS(const Match4PCSOptions& options,
                       const Utils::Logger logger);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    ~Match4PCS();

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
            Scalar pair_distance_epsilon, int base_point1,
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
} /// namespace match_4pcs

#endif  /// _4PCS_H_
