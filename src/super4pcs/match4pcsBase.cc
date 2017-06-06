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
// Authors: Dror Aiger, Yoni Weill, Nicolas Mellado
//
// This file is part of the implementation of the 4-points Congruent Sets (4PCS)
// algorithm presented in:
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

#include "match4pcsBase.h"

#include <vector>
#include <opencv2/highgui/highgui.hpp>

#include "shared4pcs.h"
#include "sampling.h"
#include "accelerators/kdtree.h"

namespace Super4PCS{

Match4PCSBase::Match4PCSBase(const match_4pcs::Match4PCSOptions& options)
  :number_of_trials_(0),
    max_base_diameter_(-1),
    P_mean_distance_(1.0),
    best_LCP_(0.0F),
    options_(options) {
  base_3D_.resize(4);
}

Match4PCSBase::Scalar
Match4PCSBase::MeanDistance() {
  const float kDiameterFraction = 0.2;
  Super4PCS::KdTree<Scalar>::VectorType query_point;

  int number_of_samples = 0;
  Scalar distance = 0.0;

  for (int i = 0; i < sampled_P_3D_.size(); ++i) {
    query_point << sampled_P_3D_[i].x, sampled_P_3D_[i].y, sampled_P_3D_[i].z;

    Super4PCS::KdTree<Scalar>::Index resId =
        kd_tree_.doQueryRestrictedClosestIndex(query_point, P_diameter_ * kDiameterFraction, i);

    if (resId != Super4PCS::KdTree<Scalar>::invalidIndex()) {
      distance += cv::norm(sampled_P_3D_[i] - sampled_P_3D_[resId]);
      number_of_samples++;
    }
  }

  return distance / number_of_samples;
}

void Match4PCSBase::init(const std::vector<Point3D>& P,
                         const std::vector<Point3D>& Q){
        const Scalar kSmallError = 0.00001;
        const int kMinNumberOfTrials = 4;
        const Scalar kDiameterFraction = 0.3;

        centroid_P_.x = 0;
        centroid_P_.y = 0;
        centroid_P_.z = 0;
        centroid_Q_.x = 0;
        centroid_Q_.y = 0;
        centroid_Q_.z = 0;

        sampled_P_3D_.clear();
        sampled_Q_3D_.clear();

        int sample_fraction_P = 1;  // We prefer not to sample P but any number can be
                                    // placed here.

        // prepare P
        if (P.size() > options_.sample_size){
            std::vector<Point3D> uniform_P;
            Super4PCS::Sampling::DistUniformSampling(P, options_.delta, &uniform_P);

            // Sample the sets P and Q uniformly.
            for (int i = 0; i < uniform_P.size(); ++i) {
              if (rand() % sample_fraction_P == 0) {
                sampled_P_3D_.push_back(uniform_P[i]);
              }
            }
        }
        else
        {
            std::cout << "(P) More samples requested than available: use whole cloud" << std::endl;
            sampled_P_3D_ = P;
        }



        // prepare Q
        if (Q.size() > options_.sample_size){
            std::vector<Point3D> uniform_Q;
            Super4PCS::Sampling::DistUniformSampling(Q, options_.delta, &uniform_Q);
            int sample_fraction_Q =
                std::max(1, static_cast<int>(uniform_Q.size() / options_.sample_size));

            for (int i = 0; i < uniform_Q.size(); ++i) {
              if (rand() % sample_fraction_Q == 0) {
                sampled_Q_3D_.push_back(uniform_Q[i]);
              }
            }
        }
        else
        {
            std::cout << "(Q) More samples requested than available: use whole cloud" << std::endl;
            sampled_Q_3D_ = Q;
        }


        // Compute the centroids.
        for (int i = 0; i < sampled_P_3D_.size(); ++i) {
          centroid_P_.x += sampled_P_3D_[i].x;
          centroid_P_.y += sampled_P_3D_[i].y;
          centroid_P_.z += sampled_P_3D_[i].z;
        }

        centroid_P_.x /= sampled_P_3D_.size();
        centroid_P_.y /= sampled_P_3D_.size();
        centroid_P_.z /= sampled_P_3D_.size();

        for (int i = 0; i < sampled_Q_3D_.size(); ++i) {
          centroid_Q_.x += sampled_Q_3D_[i].x;
          centroid_Q_.y += sampled_Q_3D_[i].y;
          centroid_Q_.z += sampled_Q_3D_[i].z;
        }

        centroid_Q_.x /= sampled_Q_3D_.size();
        centroid_Q_.y /= sampled_Q_3D_.size();
        centroid_Q_.z /= sampled_Q_3D_.size();

        // Move the samples to the centroids to allow robustness in rotation.
        for (int i = 0; i < sampled_P_3D_.size(); ++i) {
          sampled_P_3D_[i].x -= centroid_P_.x;
          sampled_P_3D_[i].y -= centroid_P_.y;
          sampled_P_3D_[i].z -= centroid_P_.z;
        }
        for (int i = 0; i < sampled_Q_3D_.size(); ++i) {
          sampled_Q_3D_[i].x -= centroid_Q_.x;
          sampled_Q_3D_[i].y -= centroid_Q_.y;
          sampled_Q_3D_[i].z -= centroid_Q_.z;
        }


        initKdTree();
        // Compute the diameter of P approximately (randomly). This is far from being
        // Guaranteed close to the diameter but gives good results for most common
        // objects if they are densely sampled.
        P_diameter_ = 0.0;
        for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
          int at = rand() % sampled_Q_3D_.size();
          int bt = rand() % sampled_Q_3D_.size();
          cv::Point3f u(sampled_Q_3D_[bt].x - sampled_Q_3D_[at].x,
                        sampled_Q_3D_[bt].y - sampled_Q_3D_[at].y,
                        sampled_Q_3D_[bt].z - sampled_Q_3D_[at].z);
          double l = cv::norm(u);
          if (l > P_diameter_) {
            P_diameter_ = l;
          }
        }

        // Mean distance and a bit more... We increase the estimation to allow for
        // noise, wrong estimation and non-uniform sampling.
        P_mean_distance_ = MeanDistance();

        // Normalize the delta (See the paper) and the maximum base distance.
        // delta = P_mean_distance_ * delta;
        max_base_diameter_ = P_diameter_;  // * estimated_overlap_;

        // RANSAC probability and number of needed trials.
        double first_estimation =
            log(kSmallError) / log(1.0 - pow(options_.overlap_estimation,
                                             static_cast<float>(kMinNumberOfTrials)));
        // We use a simple heuristic to elevate the probability to a reasonable value
        // given that we don't simply sample from P, but instead, we bound the
        // distance between the points in the base as a fraction of the diameter.
        number_of_trials_ =
            static_cast<int>(first_estimation * (P_diameter_ / kDiameterFraction) /
                             max_base_diameter_);
        if (options_.terminate_threshold < 0)
          options_.terminate_threshold = options_.overlap_estimation;
        if (number_of_trials_ < kMinNumberOfTrials)
          number_of_trials_ = kMinNumberOfTrials;

        printf("norm_max_dist: %f\n", options_.delta);
        current_trial_ = 0;
        best_LCP_ = 0.0;

        Q_copy_ = Q;
        for (int i = 0; i < 4; ++i) {
          base_[i] = 0;
          current_congruent_[i] = 0;
        }
    }


bool Match4PCSBase::SelectRandomTriangle(int* base1, int* base2, int* base3) {
      if (base1 == NULL || base2 == NULL || base3 == NULL) return false;

      int number_of_points = sampled_P_3D_.size();
      *base1 = *base2 = *base3 = -1;

      // Pick the first point at random.
      int first_point = rand() % number_of_points;

      // Try fixed number of times retaining the best other two.
      float best_wide = 0.0;
      for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
        // Pick and compute
        int second_point = rand() % number_of_points;
        int third_point = rand() % number_of_points;
        cv::Point3f u = sampled_P_3D_[second_point] - sampled_P_3D_[first_point];
        cv::Point3f w = sampled_P_3D_[third_point] - sampled_P_3D_[first_point];
        // We try to have wide triangles but still not too large.
        float how_wide = cv::norm(u.cross(w));
        if (how_wide > best_wide && cv::norm(u) < max_base_diameter_ &&
            cv::norm(w) < max_base_diameter_) {
          best_wide = how_wide;
          *base1 = first_point;
          *base2 = second_point;
          *base3 = third_point;
        }
      }
      if (*base1 == -1 || *base2 == -1 || *base3 == -1)
        return false;
      else
        return true;
}



// Try the current base in P and obtain the best pairing, i.e. the one that
// gives the smaller distance between the two closest points. The invariants
// corresponding the the base pairing are computed.
bool Match4PCSBase::TryQuadrilateral(double* invariant1, double* invariant2,
                                     int& id1, int& id2, int& id3, int& id4) {
  if (invariant1 == NULL || invariant2 == NULL) return false;

  float min_distance = FLT_MAX;
  int best1, best2, best3, best4;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (i == j) continue;
      int k = 0;
      while (k == i || k == j) k++;
      int l = 0;
      while (l == i || l == j || l == k) l++;
      double local_invariant1;
      double local_invariant2;
      // Compute the closest points on both segments, the corresponding
      // invariants and the distance between the closest points.
      float segment_distance = distSegmentToSegment(
          base_3D_[i], base_3D_[j], base_3D_[k], base_3D_[l], &local_invariant1,
          &local_invariant2);
      // Retail the smallest distance and the best order so far.
      if (segment_distance < min_distance) {
        min_distance = segment_distance;
        best1 = i;
        best2 = j;
        best3 = k;
        best4 = l;
        *invariant1 = local_invariant1;
        *invariant2 = local_invariant2;
      }
    }
  }
  std::vector<Point3D> tmp = base_3D_;
  base_3D_[0] = tmp[best1];
  base_3D_[1] = tmp[best2];
  base_3D_[2] = tmp[best3];
  base_3D_[3] = tmp[best4];

  std::array<int, 4> tmpId = {id1, id2, id3, id4};
  id1 = tmpId[best1];
  id2 = tmpId[best2];
  id3 = tmpId[best3];
  id4 = tmpId[best4];

  return true;
}


// Selects a good base from P and computes its invariants. Returns false if
// a good planar base cannot can be found.
bool Match4PCSBase::SelectQuadrilateral(double* invariant1, double* invariant2,
                                        int* base1, int* base2, int* base3,
                                        int* base4) {
  if (invariant1 == NULL || invariant2 == NULL || base1 == NULL ||
      base2 == NULL || base3 == NULL || base4 == NULL)
    return false;

  const float kBaseTooSmall = 0.2;
  int current_trial = 0;

  // Try fix number of times.
  while (current_trial < kNumberOfDiameterTrials) {
    // Select a triangle if possible. otherwise fail.
    if (!SelectRandomTriangle(base1, base2, base3)){
      return false;
    }

    base_3D_[0] = sampled_P_3D_[*base1];
    base_3D_[1] = sampled_P_3D_[*base2];
    base_3D_[2] = sampled_P_3D_[*base3];

    // The 4th point will be a one that is close to be planar to the other 3
    // while still not too close to them.
    const double& x1 = base_3D_[0].x;
    const double& y1 = base_3D_[0].y;
    const double& z1 = base_3D_[0].z;
    const double& x2 = base_3D_[1].x;
    const double& y2 = base_3D_[1].y;
    const double& z2 = base_3D_[1].z;
    const double& x3 = base_3D_[2].x;
    const double& y3 = base_3D_[2].y;
    const double& z3 = base_3D_[2].z;

    // Fit a plan.
    double denom = (-x3 * y2 * z1 + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 -
                    x2 * y1 * z3 + x1 * y2 * z3);

    if (denom != 0) {
      double A =
          (-y2 * z1 + y3 * z1 + y1 * z2 - y3 * z2 - y1 * z3 + y2 * z3) / denom;
      double B =
          (x2 * z1 - x3 * z1 - x1 * z2 + x3 * z2 + x1 * z3 - x2 * z3) / denom;
      double C =
          (-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3) / denom;
      *base4 = -1;
      double best_distance = FLT_MAX;
      // Go over all points in P.
      for (unsigned int i = 0; i < sampled_P_3D_.size(); ++i) {
        double d1 = PointsDistance(sampled_P_3D_[i], sampled_P_3D_[*base1]);
        double d2 = PointsDistance(sampled_P_3D_[i], sampled_P_3D_[*base2]);
        double d3 = PointsDistance(sampled_P_3D_[i], sampled_P_3D_[*base3]);
        float too_small = max_base_diameter_ * kBaseTooSmall;
        if (d1 >= too_small && d2 >= too_small && d3 >= too_small) {
          // Not too close to any of the first 3.
          double distance =
              std::abs(A * sampled_P_3D_[i].x + B * sampled_P_3D_[i].y +
                   C * sampled_P_3D_[i].z - 1.0);
          // Search for the most planar.
          if (distance < best_distance) {
            best_distance = distance;
            *base4 = int(i);
          }
        }
      }
      // If we have a good one we can quit.
      if (*base4 != -1) {
        base_3D_[3] = sampled_P_3D_[*base4];
        TryQuadrilateral(invariant1, invariant2, *base1, *base2, *base3, *base4);
        return true;
      }
    }
    current_trial++;
  }

  // We failed to find good enough base..
  return false;
}

void Match4PCSBase::initKdTree(){
  int number_of_points = sampled_P_3D_.size();

  // Build the kdtree.
  kd_tree_ = Super4PCS::KdTree<Scalar>(number_of_points);

  Super4PCS::KdTree<Scalar>::VectorType p;
  for (int i = 0; i < number_of_points; ++i) {
    p << sampled_P_3D_[i].x,
        sampled_P_3D_[i].y,
        sampled_P_3D_[i].z;

    kd_tree_.add(p);
  }
  kd_tree_.finalize();
}

} // namespace Super4PCS

