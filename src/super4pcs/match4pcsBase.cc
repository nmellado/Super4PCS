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
#include <opencv2/core/eigen.hpp>

#include "Eigen/Core"
#include "Eigen/Geometry"                 // MatrixBase.homogeneous()
#include "Eigen/SVD"                      // Transform.computeRotationScaling()

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
  const Scalar kDiameterFraction = 0.2;
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
          cv::Point3d u(sampled_Q_3D_[bt].x - sampled_Q_3D_[at].x,
                        sampled_Q_3D_[bt].y - sampled_Q_3D_[at].y,
                        sampled_Q_3D_[bt].z - sampled_Q_3D_[at].z);
          Scalar l = cv::norm(u);
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
        Scalar first_estimation =
            log(kSmallError) / log(1.0 - pow(options_.overlap_estimation,
                                             static_cast<Scalar>(kMinNumberOfTrials)));
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
        transform_ = Eigen::Matrix<Scalar, 4, 4>::Identity();
    }


bool Match4PCSBase::SelectRandomTriangle(int* base1, int* base2, int* base3) {
      if (base1 == NULL || base2 == NULL || base3 == NULL) return false;

      int number_of_points = sampled_P_3D_.size();
      *base1 = *base2 = *base3 = -1;

      // Pick the first point at random.
      int first_point = rand() % number_of_points;

      // Try fixed number of times retaining the best other two.
      Scalar best_wide = 0.0;
      for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
        // Pick and compute
        int second_point = rand() % number_of_points;
        int third_point = rand() % number_of_points;
        cv::Point3d u = sampled_P_3D_[second_point] - sampled_P_3D_[first_point];
        cv::Point3d w = sampled_P_3D_[third_point] - sampled_P_3D_[first_point];
        // We try to have wide triangles but still not too large.
        Scalar how_wide = cv::norm(u.cross(w));
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
bool Match4PCSBase::TryQuadrilateral(Scalar *invariant1, Scalar *invariant2,
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
bool Match4PCSBase::SelectQuadrilateral(Scalar* invariant1, Scalar* invariant2,
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

bool Match4PCSBase::TryCongruentSet(
        int base_id1, int base_id2,
        int base_id3, int base_id4,
        const std::vector<Super4PCS::Quadrilateral>& congruent_quads){
    std::array<std::pair<Point3D, Point3D>,4> congruent_points;

    // get references to the basis coordinates
    const Point3D& b1 = sampled_P_3D_[base_id1];
    const Point3D& b2 = sampled_P_3D_[base_id2];
    const Point3D& b3 = sampled_P_3D_[base_id3];
    const Point3D& b4 = sampled_P_3D_[base_id4];


    // Centroid of the basis, computed once and using only the three first points
    Eigen::Matrix<Scalar, 3, 1> centroid1;
    // Must be improved when running without opencv
    centroid1 << (b1.x + b2.x + b3.x) / Scalar(3.),
                 (b1.y + b2.y + b3.y) / Scalar(3.),
                 (b1.z + b2.z + b3.z) / Scalar(3.);

    // Centroid of the sets, computed in the loop using only the three first points
    Eigen::Matrix<Scalar, 3, 1> centroid2;

    // set the basis coordinates in the congruent quad array
    congruent_points[0].first = b1;
    congruent_points[1].first = b2;
    congruent_points[2].first = b3;
    congruent_points[3].first = b4;

    Eigen::Matrix<Scalar, 4, 4> transform;
    for (int i = 0; i < congruent_quads.size(); ++i) {
      int a = congruent_quads[i].vertices[0];
      int b = congruent_quads[i].vertices[1];
      int c = congruent_quads[i].vertices[2];
      int d = congruent_quads[i].vertices[3];
      congruent_points[0].second = sampled_Q_3D_[a];
      congruent_points[1].second = sampled_Q_3D_[b];
      congruent_points[2].second = sampled_Q_3D_[c];
      congruent_points[3].second = sampled_Q_3D_[d];

  #ifdef STATIC_BASE
      std::cout << "Ids:" << std::endl;
      std::cout << base_id1 << "\t"
                << base_id2 << "\t"
                << base_id3 << "\t"
                << base_id4 << std::endl;
      std::cout << a << "\t"
                << b << "\t"
                << c << "\t"
                << d << std::endl;
  #endif

      centroid2 << (congruent_points[0].second.x + congruent_points[1].second.x + congruent_points[2].second.x) / Scalar(3.),
                   (congruent_points[0].second.y + congruent_points[1].second.y + congruent_points[2].second.y) / Scalar(3.),
                   (congruent_points[0].second.z + congruent_points[1].second.z + congruent_points[2].second.z) / Scalar(3.);

      Scalar rms = -1;

      bool ok =
      ComputeRigidTransformation(congruent_points,   // input congruent quads
                                 centroid1,          // input: basis centroid
                                 centroid2,          // input: candidate quad centroid
                                 options_.max_angle * M_PI / 180.0, // maximum per-dimension angle, check return value to detect invalid cases
                                 transform,          // output: transformation
                                 rms,                // output: rms error of the transformation between the basis and the congruent quad
                             #ifdef MULTISCALE
                                 true
                             #else
                                 false
                             #endif
                                 );             // state: compute scale ratio ?

      if (ok && rms >= Scalar(0.)) {

        // We give more tolerantz in computing the best rigid transformation.
        if (rms < distance_factor * options_.delta) {
          // The transformation is computed from the point-clouds centered inn [0,0,0]

          // Verify the rest of the points in Q against P.
          Scalar lcp = Verify(transform);
          if (lcp > best_LCP_) {
            // Retain the best LCP and transformation.
            base_[0] = base_id1;
            base_[1] = base_id2;
            base_[2] = base_id3;
            base_[3] = base_id4;

            current_congruent_[0] = a;
            current_congruent_[1] = b;
            current_congruent_[2] = c;
            current_congruent_[3] = d;

            best_LCP_    = lcp;
            transform_   = transform;
            qcentroid1_  = centroid1;
            qcentroid2_  = centroid2;
          }
          // Terminate if we have the desired LCP already.
          if (best_LCP_ > options_.terminate_threshold){
            return true;
          }
        }
      }
    }

    // If we reached here we do not have yet the desired LCP.
    return false;
}


bool Match4PCSBase::ComputeRigidTransformation(
        const std::array< std::pair<Point3D, Point3D>,4>& pairs,
        const Eigen::Matrix<Scalar, 3, 1>& centroid1,
        Eigen::Matrix<Scalar, 3, 1> centroid2,
        Scalar max_angle,
        Eigen::Matrix<Scalar, 4, 4> &transform,
        Scalar& rms_,
        bool computeScale ) {

  rms_ = match_4pcs::kLargeNumber;

  if (pairs.size() == 0 || pairs.size() % 2 != 0)
      return false;


  Scalar kSmallNumber = 1e-6;
  cv::Mat rotation = cv::Mat::eye(3, 3, CV_64F);

  // We only use the first 3 pairs. This simplifies the process considerably
  // because it is the planar case.

  const cv::Point3d& p0 = pairs[0].first;
  const cv::Point3d& p1 = pairs[1].first;
  const cv::Point3d& p2 = pairs[2].first;
        cv::Point3d  q0 = pairs[0].second;
        cv::Point3d  q1 = pairs[1].second;
        cv::Point3d  q2 = pairs[2].second;

  Scalar scaleEst (1.);

  // Compute scale factor if needed
  if (computeScale){
      const cv::Point3d& p3 = pairs[3].first;
      const cv::Point3d& q3 = pairs[3].second;

      Scalar ratio1 = cv::norm(p1 - p0) / cv::norm(q1 - q0);
      Scalar ratio2 = cv::norm(p3 - p2) / cv::norm(q3 - q2);

      Scalar ratioDev  = std::abs(ratio1/ratio2 - Scalar(1.));  // deviation between the two
      Scalar ratioMean = (ratio1+ratio2)/Scalar(2.);            // mean of the two

      if ( ratioDev > Scalar(0.1) )
          return match_4pcs::kLargeNumber;

//      std::cout << ratio1 << " "
//                << ratio2 << " "
//                << ratioDev << " "
//                << ratioMean << std::endl;

      scaleEst = ratioMean;

      // apply scale factor to q
      q0 = q0*scaleEst;
      q1 = q1*scaleEst;
      q2 = q2*scaleEst;
      centroid2 *= scaleEst;
  }

  cv::Point3d vector_p1 = p1 - p0;
  if (cv::norm(vector_p1) == 0) return match_4pcs::kLargeNumber;
  vector_p1 = vector_p1 * (1.0 / cv::norm(vector_p1));
  cv::Point3d vector_p2 = (p2 - p0) - ((p2 - p0).dot(vector_p1)) * vector_p1;
  if (cv::norm(vector_p2) == 0) return match_4pcs::kLargeNumber;
  vector_p2 = vector_p2 * (1.0 / cv::norm(vector_p2));
  cv::Point3d vector_p3 = vector_p1.cross(vector_p2);

  cv::Point3d vector_q1 = q1 - q0;
  if (cv::norm(vector_q1) == 0) return match_4pcs::kLargeNumber;
  vector_q1 = vector_q1 * (1.0 / cv::norm(vector_q1));
  cv::Point3d vector_q2 = (q2 - q0) - ((q2 - q0).dot(vector_q1)) * vector_q1;
  if (cv::norm(vector_q2) == 0) return match_4pcs::kLargeNumber;
  vector_q2 = vector_q2 * (1.0 / cv::norm(vector_q2));
  cv::Point3d vector_q3 = vector_q1.cross(vector_q2);

  cv::Mat rotate_p(3, 3, CV_64F);
  rotate_p.at<double>(0, 0) = vector_p1.x;
  rotate_p.at<double>(0, 1) = vector_p1.y;
  rotate_p.at<double>(0, 2) = vector_p1.z;
  rotate_p.at<double>(1, 0) = vector_p2.x;
  rotate_p.at<double>(1, 1) = vector_p2.y;
  rotate_p.at<double>(1, 2) = vector_p2.z;
  rotate_p.at<double>(2, 0) = vector_p3.x;
  rotate_p.at<double>(2, 1) = vector_p3.y;
  rotate_p.at<double>(2, 2) = vector_p3.z;

  cv::Mat rotate_q(3, 3, CV_64F);
  rotate_q.at<double>(0, 0) = vector_q1.x;
  rotate_q.at<double>(0, 1) = vector_q1.y;
  rotate_q.at<double>(0, 2) = vector_q1.z;
  rotate_q.at<double>(1, 0) = vector_q2.x;
  rotate_q.at<double>(1, 1) = vector_q2.y;
  rotate_q.at<double>(1, 2) = vector_q2.z;
  rotate_q.at<double>(2, 0) = vector_q3.x;
  rotate_q.at<double>(2, 1) = vector_q3.y;
  rotate_q.at<double>(2, 2) = vector_q3.z;

  rotation = rotate_p.t() * rotate_q;

  // Discard singular solutions. The rotation should be orthogonal.
  cv::Mat unit = rotation * rotation.t();
  if (std::abs(unit.at<double>(0, 0) - 1.0) > kSmallNumber ||
      std::abs(unit.at<double>(1, 1) - 1.0) > kSmallNumber ||
      std::abs(unit.at<double>(2, 2) - 1.0) > kSmallNumber){
      return false;
  }

  // Discard too large solutions (todo: lazy evaluation during boolean computation
  Scalar theta_x = std::abs(std::atan2(rotation.at<double>(2, 1), rotation.at<double>(2, 2)));
  Scalar theta_y = std::abs(std::atan2(-rotation.at<double>(2, 0),
                             std::sqrt(std::pow(rotation.at<double>(2, 1),2) +
                                       std::pow(rotation.at<double>(2, 2),2))));
  Scalar theta_z = std::abs(atan2(rotation.at<double>(1, 0), rotation.at<double>(0, 0)));
  if (theta_x > max_angle ||
      theta_y > max_angle ||
      theta_z > max_angle)
      return false;


  // Compute rms and return it.
  rms_ = Scalar(0.0);
  {
      cv::Mat first(3, 1, CV_64F), transformed;
      for (int i = 0; i < 3; ++i) {
          first.at<double>(0, 0) = scaleEst*pairs[i].second.x - centroid2(0);
          first.at<double>(1, 0) = scaleEst*pairs[i].second.y - centroid2(1);
          first.at<double>(2, 0) = scaleEst*pairs[i].second.z - centroid2(2);
          transformed = rotation * first;
          rms_ += sqrt(std::pow(transformed.at<double>(0, 0) -
                              (pairs[i].first.x - centroid1(0)), 2) +
                       std::pow(transformed.at<double>(1, 0) -
                              (pairs[i].first.y - centroid1(1)), 2) +
                       std::pow(transformed.at<double>(2, 0) -
                              (pairs[i].first.z - centroid1(2)), 2));
      }
  }

  rms_ /= Scalar(pairs.size());

  Eigen::Transform<Scalar, 3, Eigen::Affine> etrans (Eigen::Transform<Scalar, 3, Eigen::Affine>::Identity());

  // compute rotation and translation
  {
      Eigen::Matrix<Scalar, 3, 3> rot;
      cv::cv2eigen(rotation, rot);

      //std::cout << scaleEst << endl;

      etrans.scale(scaleEst);       // apply scale factor
      etrans.translate(centroid1);  // translation between quads
      etrans.rotate(rot);           // rotate to align frames
      etrans.translate(-centroid2); // move to congruent quad frame

      transform = etrans.matrix();
  }

  return true;
}



// Verify a given transformation by computing the number of points in P at
// distance at most (normalized) delta from some point in Q. In the paper
// we describe randomized verification. We apply deterministic one here with
// early termination. It was found to be fast in practice.
Match4PCSBase::Scalar
Match4PCSBase::Verify(const Eigen::Matrix<Scalar, 4, 4>& mat) {

#ifdef TEST_GLOBAL_TIMINGS
    Timer t_verify (true);
#endif

  // We allow factor 2 scaling in the normalization.
  Scalar epsilon = options_.delta;
  int good_points = 0;
  int number_of_points = sampled_Q_3D_.size();
  int terminate_value = best_LCP_ * number_of_points;

  Scalar sq_eps = epsilon*epsilon;
  typename Point3D::VectorType p;

  for (int i = 0; i < number_of_points; ++i) {

    // Use the kdtree to get the nearest neighbor
#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif

    auto pos = sampled_Q_3D_[i];
    p << pos.x, pos.y, pos.z;

    Super4PCS::KdTree<Scalar>::Index resId =
    kd_tree_.doQueryRestrictedClosestIndex(
                (mat * p.homogeneous()).head<3>(),
                sq_eps);


#ifdef TEST_GLOBAL_TIMINGS
    kdTreeTime += Scalar(t.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif

    if ( resId != Super4PCS::KdTree<Scalar>::invalidIndex() ) {
//      Point3D& q = sampled_P_3D_[near_neighbor_index[0]];
//      bool rgb_good =
//          (p.rgb()[0] >= 0 && q.rgb()[0] >= 0)
//              ? cv::norm(p.rgb() - q.rgb()) < options_.max_color_distance
//              : true;
//      bool norm_good = norm(p.normal()) > 0 && norm(q.normal()) > 0
//                           ? fabs(p.normal().ddot(q.normal())) >= cos_dist
//                           : true;
//      if (rgb_good && norm_good) {
        good_points++;
//      }
    }

    // We can terminate if there is no longer chance to get better than the
    // current best LCP.
    if (number_of_points - i + good_points < terminate_value) {
      break;
    }
  }

#ifdef TEST_GLOBAL_TIMINGS
  verifyTime += Scalar(t_verify.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif
  return Scalar(good_points) / number_of_points;
}




// Performs N RANSAC iterations and compute the best transformation. Also,
// transforms the set Q by this optimal transformation.
bool Match4PCSBase::Perform_N_steps(int n, cv::Mat* transformation,
                                    std::vector<Point3D>* Q) {
  if (transformation == NULL || Q == NULL) return false;

#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif

  Scalar last_best_LCP = best_LCP_;
  bool ok;
  int64 t0 = clock();
  for (int i = current_trial_; i < current_trial_ + n; ++i) {
    ok = TryOneBase();

    Scalar fraction_try  = Scalar(i) / Scalar(number_of_trials_);
    Scalar fraction_time = Scalar(clock() - t0) / CLOCKS_PER_SEC /
                          options_.max_time_seconds;
    Scalar fraction = std::max(fraction_time, fraction_try);
    printf("done: %d%c best: %f                  \r",
           static_cast<int>(fraction * 100), '%', best_LCP_);
    fflush(stdout);
    // ok means that we already have the desired LCP.
    if (ok || i > number_of_trials_ || fraction >= 0.99 || best_LCP_ == 1.0) break;
  }

  current_trial_ += n;
  if (best_LCP_ > last_best_LCP) {
    *Q = Q_copy_;

      // The transformation has been computed between the two point clouds centered
    // at the origin, we need to recompute the translation to apply it to the original clouds
    {
        Eigen::Matrix<Scalar, 3,1> centroid_P,centroid_Q;
        cv::Mat first(3, 1, CV_64F);
        first.at<double>(0, 0) = centroid_P_.x;
        first.at<double>(1, 0) = centroid_P_.y;
        first.at<double>(2, 0) = centroid_P_.z;
        cv::cv2eigen(first, centroid_P);
        first.at<double>(0, 0) = centroid_Q_.x;
        first.at<double>(1, 0) = centroid_Q_.y;
        first.at<double>(2, 0) = centroid_Q_.z;
        cv::cv2eigen(first, centroid_Q);

        Eigen::Matrix<Scalar, 3, 3> rot, scale;
        Eigen::Transform<Scalar, 3, Eigen::Affine> (transform_).computeRotationScaling(&rot, &scale);
        transform_.col(3) = (qcentroid1_ + centroid_P - ( rot * scale * (qcentroid2_ + centroid_Q))).homogeneous();
    }

    cv::eigen2cv( transform_, *transformation );

    // Transforms Q by the new transformation.
    for (unsigned int i = 0; i < Q->size(); ++i) {
      cv::Mat first(4, 1, CV_64F), transformed;
      first.at<double>(0, 0) = (*Q)[i].x;
      first.at<double>(1, 0) = (*Q)[i].y;
      first.at<double>(2, 0) = (*Q)[i].z;
      first.at<double>(3, 0) = 1;
      transformed = *transformation * first;
      (*Q)[i].x = transformed.at<double>(0, 0);
      (*Q)[i].y = transformed.at<double>(1, 0);
      (*Q)[i].z = transformed.at<double>(2, 0);
    }
  }
#ifdef TEST_GLOBAL_TIMINGS
    totalTime += Scalar(t.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif

  return ok || current_trial_ >= number_of_trials_;
}

} // namespace Super4PCS

