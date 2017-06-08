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


#include "4pcs.h"
#include "match4pcsBase.h"
#include "accelerators/utils.h"

#include "ANN/ANN.h"
#include <opencv2/highgui/highgui.hpp>

#include <fstream>
#include <time.h>  //clock

namespace match_4pcs {

Match4PCS::Match4PCS(const Match4PCSOptions& options)
    : Base(options) { }

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
        std::vector<Super4PCS::Quadrilateral>* quadrilaterals) const {
  if (quadrilaterals == NULL) return false;

  int number_of_points = 2 * P_pairs.size();

  // We need a temporary ANN tree to store the new points corresponding to
  // invariants in the P_pairs and then query them (for range search) for all
  // the new points corresponding to the invariants in Q_pairs.
  ANNpointArray data_points = annAllocPts(number_of_points, 3);
  ANNpoint query_point = annAllocPt(3);
  ANNidxArray near_neighbor_index = new ANNidx[number_of_points];
  ANNdistArray distances = new ANNdist[number_of_points];

  quadrilaterals->clear();

  // Build the ANN tree using the invariants on P_pairs.
  for (int i = 0; i < P_pairs.size(); ++i) {
    const Point3D& p1 = sampled_Q_3D_[P_pairs[i].first];
    const Point3D& p2 = sampled_Q_3D_[P_pairs[i].second];
    data_points[i * 2][0] = p1.x + invariant1 * (p2.x - p1.x);
    data_points[i * 2][1] = p1.y + invariant1 * (p2.y - p1.y);
    data_points[i * 2][2] = p1.z + invariant1 * (p2.z - p1.z);
    data_points[i * 2 + 1][0] = p1.x + (1.0-invariant1) * (p2.x - p1.x);
    data_points[i * 2 + 1][1] = p1.y + (1.0-invariant1) * (p2.y - p1.y);
    data_points[i * 2 + 1][2] = p1.z + (1.0-invariant1) * (p2.z - p1.z);
  }

  ANNkd_tree* tree = new ANNkd_tree(data_points, number_of_points, 3);

    //Point3D invRes;
  // Query the ANN for all the points corresponding to the invariants in Q_pair.
  for (int i = 0; i < Q_pairs.size(); ++i) {
    const Point3D& p1 = sampled_Q_3D_[Q_pairs[i].first];
    const Point3D& p2 = sampled_Q_3D_[Q_pairs[i].second];
    query_point[0] = p1.x + invariant2 * (p2.x - p1.x);
    query_point[1] = p1.y + invariant2 * (p2.y - p1.y);
    query_point[2] = p1.z + invariant2 * (p2.z - p1.z);

    tree->annkFRSearch(query_point, distance_threshold2, number_of_points,
                       near_neighbor_index, distances, 0);

    // This is a new candidate of a quadrilateral.
    for (int j = 0; j < number_of_points; ++j) {
      if (distances[j] != ANN_DIST_INF) {
        int id = near_neighbor_index[j] / 2;

        const Point3D& pp1 = sampled_Q_3D_[P_pairs[id].first];
        const Point3D& pp2 = sampled_Q_3D_[P_pairs[id].second];

        quadrilaterals->push_back(
            Super4PCS::Quadrilateral(P_pairs[id].first, P_pairs[id].second,
                          Q_pairs[i].first, Q_pairs[i].second));
      } else
        break;
    }

    // We test the other order as our pairs are not ordered.
    query_point[0] = p1.x + (1.0-invariant2) * (p2.x - p1.x);
    query_point[1] = p1.y + (1.0-invariant2) * (p2.y - p1.y);
    query_point[2] = p1.z + (1.0-invariant2) * (p2.z - p1.z);

    tree->annkFRSearch(query_point, distance_threshold2, number_of_points,
                       near_neighbor_index, distances, 0);

    for (int j = 0; j < number_of_points; ++j) {
      if (distances[j] != ANN_DIST_INF) {
        int id = near_neighbor_index[j] / 2;

        const Point3D& pp1 = sampled_Q_3D_[P_pairs[id].first];
        const Point3D& pp2 = sampled_Q_3D_[P_pairs[id].second];

        quadrilaterals->push_back(
            Super4PCS::Quadrilateral(P_pairs[id].first, P_pairs[id].second,
                          Q_pairs[i].first, Q_pairs[i].second));
      } else
        break;
    }
  }

  annDeallocPt(query_point);
  annDeallocPts(data_points);
  delete[] near_neighbor_index;
  delete[] distances;
  delete tree;

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
  if (pairs == NULL) return;

  pairs->clear();
  pairs->reserve(2 * sampled_Q_3D_.size());

  cv::Point3d segment1 = base_3D_[base_point2] - base_3D_[base_point1];
  segment1 *= 1.0 / cv::norm(segment1);


  // Go over all ordered pairs in Q.
  for (int j = 0; j < sampled_Q_3D_.size(); ++j) {
    const Point3D& p = sampled_Q_3D_[j];
    for (int i = j + 1; i < sampled_Q_3D_.size(); ++i) {
      const Point3D& q = sampled_Q_3D_[i];
      // Compute the distance and two normal angles to ensure working with
      // wrong orientation. We want to verify that the angle between the
      // normals is close to the angle between normals in the base. This can be
      // checked independent of the full rotation angles which are not yet
      // defined by segment matching alone..
      const Scalar distance = cv::norm(q - p);
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

        if (!first_norm_distance < norm_threshold) continue;
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

          if (use_rgb && ! color_good) continue;
      }

      if (options_.max_translation_distance > 0) {
          const bool dist_good = cv::norm(p - base_3D_[base_point1]) <
                  options_.max_translation_distance &&
                  cv::norm(q - base_3D_[base_point2]) <
                  options_.max_translation_distance;
          if (! dist_good) continue;
      }

      if (options_.max_angle > 0){
          cv::Point3d segment2 = q - p;
          segment2 *= 1.0 / cv::norm(segment2);
          if (acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0) {
            pairs->push_back(std::make_pair(j, i));
          }

          segment2 = p - q;
          segment2 *= 1.0 / cv::norm(segment2);
          if (acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0) {
            // Add ordered pair.
            pairs->push_back(std::make_pair(i, j));
          }
      }else {
          pairs->push_back(std::make_pair(j, i));
          pairs->push_back(std::make_pair(i, j));
      }
    }
  }
}


// Initialize all internal data structures and data members.
void Match4PCS::Initialize(const std::vector<Point3D>& /*P*/,
                               const std::vector<Point3D>& /*Q*/) {}

}
