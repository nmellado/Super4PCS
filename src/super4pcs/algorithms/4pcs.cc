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


#include "super4pcs/algorithms/4pcs.h"
#include "super4pcs/algorithms/match4pcsBase.h"
#include "super4pcs/accelerators/utils.h"

#include <fstream>
#include <time.h>  //clock

namespace GlobalRegistration {

Match4PCS::Match4PCS(const Match4PCSOptions& options,
                     const Utils::Logger logger)
    : Base(options, logger) { }

Match4PCS::~Match4PCS() { }

// Finds congruent candidates in the set Q, given the invariants and threshold
// distances.
bool Match4PCS::FindCongruentQuadrilaterals(
        Scalar invariant1,
        Scalar invariant2,
        Scalar /*distance_threshold1*/,
        Scalar distance_threshold2,
        const std::vector<std::pair<int, int>>& P_pairs,
        const std::vector<std::pair<int, int>>& Q_pairs,
        std::vector<GlobalRegistration::Quadrilateral>* quadrilaterals) const {
  using RangeQuery = GlobalRegistration::KdTree<Scalar>::RangeQuery<>;

  if (quadrilaterals == nullptr) return false;

  size_t number_of_points = 2 * P_pairs.size();

  // We need a temporary kdtree to store the new points corresponding to
  // invariants in the P_pairs and then query them (for range search) for all
  // the new points corresponding to the invariants in Q_pairs.
  quadrilaterals->clear();

  GlobalRegistration::KdTree<Scalar> kdtree (number_of_points);

  // Build the kdtree tree using the invariants on P_pairs.
  for (size_t i = 0; i < P_pairs.size(); ++i) {
    const VectorType& p1 = sampled_Q_3D_[P_pairs[i].first].pos();
    const VectorType& p2 = sampled_Q_3D_[P_pairs[i].second].pos();
    kdtree.add(p1 + invariant1 * (p2-p1));
  }
  kdtree.finalize();

    //Point3D invRes;
  // Query the Kdtree for all the points corresponding to the invariants in Q_pair.
  for (size_t i = 0; i < Q_pairs.size(); ++i) {
    const VectorType& p1 = sampled_Q_3D_[Q_pairs[i].first].pos();
    const VectorType& p2 = sampled_Q_3D_[Q_pairs[i].second].pos();

    RangeQuery query;
    query.queryPoint = p1 + invariant2 * (p2 - p1);
    query.sqdist     = distance_threshold2;

    kdtree.doQueryDistProcessIndices(query,
        [quadrilaterals, i, &P_pairs, &Q_pairs](int id){
        quadrilaterals->emplace_back(P_pairs[id/2].first, P_pairs[id/2].second,
                                     Q_pairs[i].first, Q_pairs[i].second);
    });
  }

  return quadrilaterals->size() != 0;
}

// Constructs two sets of pairs in Q, each corresponds to one pair in the base
// in P, by having the same distance (up to some tolerantz) and optionally the
// same angle between normals and same color.
void
Match4PCS::ExtractPairs(Scalar pair_distance,
                        Scalar pair_normals_angle,
                        Scalar pair_distance_epsilon,
                        int base_point1,
                        int base_point2,
                        PairsVector* pairs) const {
  if (pairs == nullptr) return;

  pairs->clear();
  pairs->reserve(2 * sampled_Q_3D_.size());

  VectorType segment1 = (base_3D_[base_point2].pos() -
                         base_3D_[base_point1].pos()).normalized();


  // Go over all ordered pairs in Q.
  for (size_t j = 0; j < sampled_Q_3D_.size(); ++j) {
    const Point3D& p = sampled_Q_3D_[j];
    for (size_t i = j + 1; i < sampled_Q_3D_.size(); ++i) {
      const Point3D& q = sampled_Q_3D_[i];
      // Compute the distance and two normal angles to ensure working with
      // wrong orientation. We want to verify that the angle between the
      // normals is close to the angle between normals in the base. This can be
      // checked independent of the full rotation angles which are not yet
      // defined by segment matching alone..
      const Scalar distance = (q.pos() - p.pos()).norm();
#ifndef MULTISCALE
      if (std::abs(distance - pair_distance) > pair_distance_epsilon) continue;
#endif

      if ( options_.max_normal_difference > 0 &&
              q.normal().squaredNorm() > 0 &&
              p.normal().squaredNorm() > 0) {
        const Scalar norm_threshold =
              0.5 * options_.max_normal_difference * M_PI / 180.0;
        const double first_normal_angle = (q.normal() - p.normal()).norm();
        const double second_normal_angle = (q.normal() + p.normal()).norm();
        // Take the smaller normal distance.
        const Scalar first_norm_distance =
            std::min(std::abs(first_normal_angle - pair_normals_angle),
                std::abs(second_normal_angle - pair_normals_angle));
        // Verify appropriate angle between normals and distance.

        if (first_norm_distance > norm_threshold) continue;
      }
      // Verify restriction on the rotation angle, translation and colors.
      if (options_.max_color_distance > 0) {
          const bool use_rgb = (p.rgb()[0] >= 0 && q.rgb()[0] >= 0 &&
                  base_3D_[base_point1].rgb()[0] >= 0 &&
                  base_3D_[base_point2].rgb()[0] >= 0);
          bool color_good = (p.rgb() - base_3D_[base_point1].rgb()).norm() <
                  options_.max_color_distance &&
                  (q.rgb() - base_3D_[base_point2].rgb()).norm() <
                  options_.max_color_distance;

          if (use_rgb && ! color_good) return;
      }

      if (options_.max_translation_distance > 0) {
          const bool dist_good = (p.pos() - base_3D_[base_point1].pos()).norm() <
                  options_.max_translation_distance &&
                  (q.pos() - base_3D_[base_point2].pos()).norm() <
                  options_.max_translation_distance;
          if (! dist_good) return;
      }

      // need cleaning here
      if (options_.max_angle > 0){
          VectorType segment2 = (q.pos() - p.pos()).normalized();
          if (std::acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0) {
              pairs->emplace_back(j, i);
          }

          if (std::acos(segment1.dot(- segment2)) <= options_.max_angle * M_PI / 180.0) {
              // Add ordered pair.
              pairs->emplace_back(i, j);
          }
      }else {
          pairs->emplace_back(j, i);
          pairs->emplace_back(i, j);
      }
    }
  }
}


// Initialize all internal data structures and data members.
void Match4PCS::Initialize(const std::vector<Point3D>& /*P*/,
                               const std::vector<Point3D>& /*Q*/) {}

}
