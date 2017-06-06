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

using namespace std;

namespace {


// Computes the best rigid transformation between three corresponding pairs.
// The transformation is characterized by rotation matrix, translation vector
// and a center about which we rotate. The set of pairs is potentially being
// updated by the best permutation of the second set. Returns the RMS of the
// fit. The method is being called with 4 points but it applies the fit for
// only 3 after the best permutation is selected in the second set (see
// bellow). This is done because the solution for planar points is much
// simpler.
// The method is the closed-form solution by Horn:
// people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
double ComputeRigidTransformation(vector<pair<Point3D, Point3D>>* pairs,
                                  cv::Mat* rotation, cv::Point3f* translation,
                                  cv::Point3f* center) {
  if (pairs->size() == 0 || rotation == NULL || translation == NULL ||
      center == NULL)
    return kLargeNumber;
  float kSmallNumber = 1e-6;
  *rotation = cv::Mat::eye(3, 3, CV_64F);

  // Compute the centroid of the sets.
  cv::Point3f centroid1(0, 0, 0);
  cv::Point3f centroid2(0, 0, 0);
  for (int i = 0; i < pairs->size(); ++i) {
    centroid1 += (*pairs)[i].first;
    centroid2 += (*pairs)[i].second;
  }
  centroid1 *= 1.0 / static_cast<float>(pairs->size());
  centroid2 *= 1.0 / static_cast<float>(pairs->size());

  // Search for the best mapping by comparing the distances of a points to the
  // center of gravity in both sets.
  vector<pair<Point3D, Point3D>> temp_pairs = *pairs;
  for (int i = 0; i < pairs->size(); ++i) {
    float distance1 = cv::norm((*pairs)[i].first - centroid1);
    double best_distance = FLT_MAX;
    int best_id;
    for (int j = 0; j < pairs->size(); ++j) {
      float distance2 = cv::norm((*pairs)[j].second - centroid2);
      float t = fabs(distance1 - distance2);
      if (t < best_distance) {
        best_distance = t;
        best_id = j;
      }
    }
    (*pairs)[i].second = temp_pairs[best_id].second;
  }

  // We only use the first 3 pairs. This simplifies the process considerably
  // because it is the planar case.

  const cv::Point3f& p0 = (*pairs)[0].first;
  const cv::Point3f& p1 = (*pairs)[1].first;
  const cv::Point3f& p2 = (*pairs)[2].first;
  const cv::Point3f& q0 = (*pairs)[0].second;
  const cv::Point3f& q1 = (*pairs)[1].second;
  const cv::Point3f& q2 = (*pairs)[2].second;

  centroid1 = p0 + p1 + p2;
  centroid1 *= 1.0 / 3.0;
  centroid2 = q0 + q1 + q2;
  centroid2 *= 1.0 / 3.0;

  cv::Point3f vector_p1 = p1 - p0;
  if (cv::norm(vector_p1) == 0) return kLargeNumber;
  vector_p1 = vector_p1 * (1.0 / cv::norm(vector_p1));
  cv::Point3f vector_p2 = (p2 - p0) - ((p2 - p0).dot(vector_p1)) * vector_p1;
  if (cv::norm(vector_p2) == 0) return kLargeNumber;
  vector_p2 = vector_p2 * (1.0 / cv::norm(vector_p2));
  cv::Point3f vector_p3 = vector_p1.cross(vector_p2);

  cv::Point3f vector_q1 = q1 - q0;
  if (cv::norm(vector_q1) == 0) return kLargeNumber;
  vector_q1 = vector_q1 * (1.0 / cv::norm(vector_q1));
  cv::Point3f vector_q2 = (q2 - q0) - ((q2 - q0).dot(vector_q1)) * vector_q1;
  if (cv::norm(vector_q2) == 0) return kLargeNumber;
  vector_q2 = vector_q2 * (1.0 / cv::norm(vector_q2));
  cv::Point3f vector_q3 = vector_q1.cross(vector_q2);

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

  *rotation = rotate_p.t() * rotate_q;

  // Discard singular solutions. The rotation should be orthogonal.
  cv::Mat unit = *rotation * rotation->t();
  if (fabs(unit.at<double>(0, 0) - 1.0) > kSmallNumber ||
      fabs(unit.at<double>(1, 1) - 1.0) > kSmallNumber ||
      fabs(unit.at<double>(2, 2) - 1.0) > kSmallNumber)
    return kLargeNumber;

  *center = centroid2;
  *translation = centroid1 - centroid2;

  cv::Mat first(3, 1, CV_64F), transformed;
  // Compute rms and return it.
  double rms = 0.0;
  for (int i = 0; i < 3; ++i) {
    first.at<double>(0, 0) = (*pairs)[i].second.x - centroid2.x;
    first.at<double>(1, 0) = (*pairs)[i].second.y - centroid2.y;
    first.at<double>(2, 0) = (*pairs)[i].second.z - centroid2.z;
    transformed = *rotation * first;
    rms += sqrt(Square(transformed.at<double>(0, 0) -
                       ((*pairs)[i].first.x - centroid1.x)) +
                Square(transformed.at<double>(1, 0) -
                       ((*pairs)[i].first.y - centroid1.y)) +
                Square(transformed.at<double>(2, 0) -
                       ((*pairs)[i].first.z - centroid1.z)));
  }

  return rms / pairs->size();
}

// Transforms a point by a given transformation.
void Transform(const cv::Mat& rotation, const cv::Point3f& center,
               const cv::Point3f& translate, Point3D* point) {
  if (point == NULL) return;
  cv::Mat first(3, 1, CV_64F), transformed;
  first.at<double>(0, 0) = (point->x - center.x);
  first.at<double>(1, 0) = (point->y - center.y);
  first.at<double>(2, 0) = (point->z - center.z);
  transformed = rotation * first;
  point->x = transformed.at<double>(0, 0) + center.x + translate.x;
  point->y = transformed.at<double>(1, 0) + center.y + translate.y;
  point->z = transformed.at<double>(2, 0) + center.z + translate.z;

  first.at<double>(0, 0) = (point->normal()(0));
  first.at<double>(1, 0) = (point->normal()(1));
  first.at<double>(2, 0) = (point->normal()(2));
  transformed = rotation * first;
  typename Point3D::VectorType normal;
  normal << transformed.at<double>(0, 0), transformed.at<double>(1, 0), transformed.at<double>(2, 0);
  point->set_normal(normal);
}


}  // namespace

class Match4PCSImpl : public Super4PCS::Match4PCSBase {
 public:
  using Base = Super4PCS::Match4PCSBase;
  explicit Match4PCSImpl(const Match4PCSOptions& options)
        : Base (options),
        ann_tree_(0),
        data_points_() {}

  ~Match4PCSImpl() {
    // Release the ANN data structure and points.
    Clear();
  }

  void Clear() {
    delete ann_tree_;
    ann_tree_ = NULL;
    if (data_points_) annDeallocPts(data_points_);
    data_points_ = NULL;
    annClose();
  }
  // Computes an approximation of the best LCP (directional) from Q to P
  // and the rigid transformation that realizes it. The input sets may or may
  // not contain normal information for any point.
  // @param [in] P The first input set.
  // @param [in] Q The second input set.
  // as a fraction of the size of P ([0..1]).
  // @param [out] transformation Rigid transformation matrix (4x4) that brings
  // Q to the (approximate) optimal LCP.
  // @return the computed LCP measure.
  // The method updates the coordinates of the second set, Q, applying
  // the found transformation.
  float ComputeTransformation(const std::vector<Point3D>& P,
                              std::vector<Point3D>* Q, cv::Mat* transformation);

 private:
  // Private data contains parameters and internal variables that are computed
  // and change during the match computation. All parameters have default
  // values.

  // Internal data members.

  // ANN structure allows to query arbitrary point for range searching.
  ANNkd_tree* ann_tree_;
  // Holds the ANN data points.
  ANNpointArray data_points_;

  // Private member functions.

  // Tries one base and finds the best transformation for this base.
  // Returns true if the achieved LCP is greater than terminate_threshold_,
  // else otherwise.
  bool TryOneBase();

  // Constructs pairs of points in Q, corresponding to a single pair in the
  // in basein P.
  // @param [in] pair_distance The distance between the pairs in P that we have
  // to match in the pairs we select from Q.
  // @param [in] pair_normal_distance The angle between the normals of the pair
  // in P.
  // @param [in] pair_distance_epsilon Tolerance on the pair distance. We allow
  // candidate pair in Q to have distance of
  // pair_distance+-pair_distance_epsilon.
  // @param [in] base_point1 The index of the first point in P.
  // @param [in] base_point2 The index of the second point in P.
  // @param [out] pairs A set of pairs in Q that match the pair in P with
  // respect to distance and normals, up to the given tolerance.
  void BruteForcePairs(double pair_distance, double pair_normals_angle,
                       double pair_distance_epsilon, int base_point1,
                       int base_point2, Base::PairsVector* pairs);

  // For each randomly picked base, verifies the computed transformation by
  // computing the number of points that this transformation brings near points
  // in Q. Returns the current LCP. R is the rotation matrix, (tx,ty,tz) is
  // the translation vector and (cx,cy,cz) is the center of transformation.
  double Verify(const cv::Mat& rotation, const cv::Point3f& center,
                const cv::Point3f& translate);

  // Finds congruent candidates in the set Q, given the invariants and threshold
  // distances. Returns true if a non empty set can be found, false otherwise.
  // @param invariant1 [in] The first invariant corresponding to the set P_pairs
  // of pairs, previously extracted from Q.
  // @param invariant2 [in] The second invariant corresponding to the set
  // Q_pairs of pairs, previously extracted from Q.
  // @param [in] distance_threshold1 The distance for verification.
  // @param [in] distance_threshold2 The distance for matching middle points due
  // to the invariants (See the paper for e1, e2).
  // @param [in] P_pairs The first set of pairs.
  // @param [in] Q_pairs The second set of pairs.
  // @param [out] quadrilaterals The set of congruent quadrilateral. In fact,
  // it's a super set from which we extract the real congruent set.
  bool FindCongruentQuadrilaterals(double invariant1, double invariant2,
                                   double distance_threshold1,
                                   double distance_threshold2,
                                   const PairsVector& P_pairs,
                                   const PairsVector& Q_pairs,
                                   std::vector<Super4PCS::Quadrilateral>* quadrilaterals);

  // Initializes the data structures and needed values before the match
  // computation.
  // @param [in] point_P First input set.
  // @param [in] point_Q Second input set.
  // expected to be in the inliers.
  void Initialize(const std::vector<Point3D>& P, const std::vector<Point3D>& Q);

  // Performs n RANSAC iterations, each one of them containing base selection,
  // finding congruent sets and verification. Returns true if the process can be
  // terminated (the target LCP was obtained or the maximum number of trials has
  // been reached), false otherwise.
  bool Perform_N_steps(int n, cv::Mat* transformation, std::vector<Point3D>* Q);
};

// Finds congruent candidates in the set Q, given the invariants and threshold
// distances.
bool Match4PCSImpl::FindCongruentQuadrilaterals(
    double invariant1, double invariant2, double distance_threshold1,
    double distance_threshold2, const std::vector<std::pair<int, int>>& P_pairs,
    const std::vector<std::pair<int, int>>& Q_pairs,
    std::vector<Super4PCS::Quadrilateral>* quadrilaterals) {
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

// Verify a given transformation by computing the number of points in P at
// distance at most (normalized) delta from some point in Q. In the paper
// we describe randomized verification. We apply deterministic one here with
// early termination. It was found to be fast in practice.
double Match4PCSImpl::Verify(const cv::Mat& rotation, const cv::Point3f& center,
                             const cv::Point3f& translate) {

  // We allow factor 2 scaling in the normalization.
  float epsilon = options_.delta;
  int good_points = 0;
  int number_of_points = sampled_Q_3D_.size();
  int terminate_value = best_LCP_ * number_of_points;

  ANNpoint query_point;
  ANNidxArray near_neighbor_index;
  ANNdistArray distances;
  query_point = annAllocPt(3);
  near_neighbor_index = new ANNidx[1];
  distances = new ANNdist[1];

  const float cos_dist =
      options_.max_normal_difference > 90
          ? 0
          : cos(min(90.0, options_.max_normal_difference) * M_PI / 180.0);

  for (int i = 0; i < number_of_points; ++i) {
    Point3D p = sampled_Q_3D_[i];
    Transform(rotation, center, translate, &p);

    query_point[0] = p.x;
    query_point[1] = p.y;
    query_point[2] = p.z;
    // Use the ANN tree to get the nearest neighbor (we use the exact version).
    ann_tree_->annkSearch(query_point, 1, near_neighbor_index, distances, 0);
    if (sqrt(distances[0]) < epsilon) {
      Point3D& q = sampled_P_3D_[near_neighbor_index[0]];
      bool rgb_good =
          (p.rgb()[0] >= 0 && q.rgb()[0] >= 0)
              ? (p.rgb() - q.rgb()).norm() < options_.max_color_distance
              : true;
      bool norm_good = p.normal().squaredNorm() > 0 && q.normal().squaredNorm() > 0
                           ? std::abs(p.normal().dot(q.normal())) >= cos_dist
                           : true;
      if (rgb_good && norm_good) {
        good_points++;
      }
    }
    // We can terminate if there is no longer chance to get better than the
    // current best LCP.
    if (number_of_points - i + good_points < terminate_value) {
      break;
    }
  }

  annDeallocPt(query_point);
  delete[] near_neighbor_index;
  delete[] distances;

  return static_cast<float>(good_points) / number_of_points;
}

// Constructs two sets of pairs in Q, each corresponds to one pair in the base
// in P, by having the same distance (up to some tolerantz) and optionally the
// same angle between normals and same color.
void Match4PCSImpl::BruteForcePairs(double pair_distance,
                                    double pair_normals_angle,
                                    double pair_distance_epsilon,
                                    int base_point1, int base_point2,
                                    vector<pair<int, int>>* pairs) {
  if (pairs == NULL) return;

  pairs->clear();
  pairs->reserve(2 * sampled_Q_3D_.size());

  cv::Point3f segment1 = base_3D_[base_point2] - base_3D_[base_point1];
  segment1 *= 1.0 / cv::norm(segment1);

  const float norm_threshold =
      0.5 * options_.max_normal_difference * M_PI / 180.0;

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
      const float distance = cv::norm(q - p);
      if (fabs(distance - pair_distance) > pair_distance_epsilon) continue;
      const bool use_normals = q.normal().squaredNorm() > 0 && p.normal().squaredNorm() > 0;
      bool normals_good = true;
      if (use_normals) {
        const double first_normal_angle = (q.normal() - p.normal()).norm();
        const double second_normal_angle = (q.normal() + p.normal()).norm();
        // Take the smaller normal distance.
        const float first_norm_distance =
            min(fabs(first_normal_angle - pair_normals_angle),
                fabs(second_normal_angle - pair_normals_angle));
        // Verify appropriate angle between normals and distance.
        normals_good = first_norm_distance < norm_threshold;
      }
      if (!normals_good) continue;
      cv::Point3f segment2 = q - p;
      segment2 *= 1.0 / cv::norm(segment2);
      // Verify restriction on the rotation angle, translation and colors.
      const bool use_rgb = (p.rgb()[0] >= 0 && q.rgb()[0] >= 0 &&
                            base_3D_[base_point1].rgb()[0] >= 0 &&
                            base_3D_[base_point2].rgb()[0] >= 0);
      const bool rgb_good =
          use_rgb ? (p.rgb() - base_3D_[base_point1].rgb()).norm() <
                            options_.max_color_distance &&
                    (q.rgb() - base_3D_[base_point2].rgb()).norm() <
                            options_.max_color_distance
                  : true;
      const bool dist_good = cv::norm(p - base_3D_[base_point1]) <
                                 options_.max_translation_distance &&
                             cv::norm(q - base_3D_[base_point2]) <
                                 options_.max_translation_distance;
      if (acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0 &&
          dist_good && rgb_good) {
        // Add ordered pair.
        pairs->push_back(pair<int, int>(j, i));
      }
      // The same for the second order.
      segment2 = p - q;
      segment2 *= 1.0 / cv::norm(segment2);
      if (acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0 &&
          dist_good && rgb_good) {
        // Add ordered pair.
        pairs->push_back(pair<int, int>(i, j));
      }
    }
  }
}

// Pick one base, finds congruent 4-points in Q, verifies for all
// transformations, and retains the best transformation and LCP. This is
// a complete RANSAC iteration.
bool Match4PCSImpl::TryOneBase() {
  vector<pair<Point3D, Point3D>> congruent_points(4);
  double invariant1, invariant2;
  int base_id1, base_id2, base_id3, base_id4;
  float distance_factor = 2.0;

  if (!SelectQuadrilateral(&invariant1, &invariant2, &base_id1, &base_id2,
                           &base_id3, &base_id4)) {
    return false;
  }

  // Computes distance between pairs.
  double distance1 = PointsDistance(base_3D_[0], base_3D_[1]);
  double distance2 = PointsDistance(base_3D_[2], base_3D_[3]);

  vector<pair<int, int>> pairs1, pairs2;
  vector<Super4PCS::Quadrilateral> congruent_quads;

  // Compute normal angles.
  double normal_angle1 = (base_3D_[0].normal() - base_3D_[1].normal()).norm();
  double normal_angle2 = (base_3D_[2].normal() - base_3D_[3].normal()).norm();

  BruteForcePairs(distance1, normal_angle1, distance_factor * options_.delta, 0,
                  1, &pairs1);
  BruteForcePairs(distance2, normal_angle2, distance_factor * options_.delta, 2,
                  3, &pairs2);
  if (pairs1.size() == 0 || pairs2.size() == 0) {
    return false;
  }
  //cout << pairs1.size() << " " << pairs1.size() << endl;
  if (!FindCongruentQuadrilaterals(invariant1, invariant2,
                                   distance_factor * options_.delta,
                                   distance_factor * options_.delta, pairs1,
                                   pairs2, &congruent_quads)) {
    return false;
  }

  //cout << "congruent_quads.size() = " << congruent_quads.size()  << endl;

  cv::Mat rotation(3, 3, CV_64F);
  for (int i = 0; i < congruent_quads.size(); ++i) {
    congruent_points.resize(4);
    int a = congruent_quads[i].vertices[0];
    int b = congruent_quads[i].vertices[1];
    int c = congruent_quads[i].vertices[2];
    int d = congruent_quads[i].vertices[3];
    congruent_points[0].first = sampled_P_3D_[base_id1];
    congruent_points[0].second = sampled_Q_3D_[a];
    congruent_points[1].first = sampled_P_3D_[base_id2];
    congruent_points[1].second = sampled_Q_3D_[b];
    congruent_points[2].first = sampled_P_3D_[base_id3];
    congruent_points[2].second = sampled_Q_3D_[c];
    congruent_points[3].first = sampled_P_3D_[base_id4];
    congruent_points[3].second = sampled_Q_3D_[d];


    cv::Point3f center;
    cv::Point3f translate;
    double f = ComputeRigidTransformation(&congruent_points, &rotation,
                                          &translate, &center);

    float theta_x =
        std::abs(atan2(rotation.at<double>(2, 1), rotation.at<double>(2, 2)));
    float theta_y = std::abs(atan2(-rotation.at<double>(2, 0),
                               sqrt(Square(rotation.at<double>(2, 1)) +
                                    Square(rotation.at<double>(2, 2)))));
    float theta_z =
        std::abs(atan2(rotation.at<double>(1, 0), rotation.at<double>(0, 0)));

    // Check angle limitation.
    if (theta_x <= options_.max_angle * M_PI / 180.0 &&
        theta_y <= options_.max_angle * M_PI / 180.0 &&
        theta_z <= options_.max_angle * M_PI / 180.0) {


      // We give more tolerantz in computing the best rigid transformation.
      if (f < distance_factor * options_.delta) {
        // Verify the rest of the points in Q against P.
        f = Verify(rotation, center, translate);
        if (f > best_LCP_) {
          // Retain the best LCP and transformation.
          base_[0] = base_id1;
          base_[1] = base_id2;
          base_[2] = base_id3;
          base_[3] = base_id4;

          current_congruent_[0] = a;
          current_congruent_[1] = b;
          current_congruent_[2] = c;
          current_congruent_[3] = d;

          best_LCP_ = f;
          rotate_ = rotation.clone();
          centroid_ = center;
          translate_ = translate;
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

struct eqstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1, s2) == 0;
  }
};

// Initialize all internal data structures and data members.
void Match4PCSImpl::Initialize(const std::vector<Point3D>& P,
                               const std::vector<Point3D>& Q) {

  Base::init(P, Q);


  // Build the ANN tree.
  int number_of_points = sampled_P_3D_.size();
  if (data_points_) {
    annDeallocPts(data_points_);
  }
  if (ann_tree_) delete ann_tree_;

  data_points_ = annAllocPts(number_of_points, 3);
  for (int i = 0; i < sampled_P_3D_.size(); ++i) {
    data_points_[i][0] = sampled_P_3D_[i].x;
    data_points_[i][1] = sampled_P_3D_[i].y;
    data_points_[i][2] = sampled_P_3D_[i].z;
  }
  ann_tree_ = new ANNkd_tree(data_points_, number_of_points, 3);

  cv::Mat rotation = cv::Mat::eye(3, 3, CV_64F);
  cv::Vec3f translate(0, 0, 0);
  cv::Vec3f center(0, 0, 0);
  best_LCP_ = Verify(rotation, center, translate);
  printf("Initial LCP: %f\n", best_LCP_);
}

// Performs N RANSAC iterations and compute the best transformation. Also,
// transforms the set Q by this optimal transformation.
bool Match4PCSImpl::Perform_N_steps(int n, cv::Mat* transformation,
                                    std::vector<Point3D>* Q) {

#ifdef TEST_RECORD_CONVERGENCE
  ofstream myfile;
  myfile.open ("convergence.txt", ios::out | ios::app);
  float last_best_LCPTrace = best_LCP_;

  myfile << 0 << " " << best_LCP_ << " " << 0 << endl;
#endif

  float last_best_LCP = best_LCP_;
  //vector<Point3D> cpy = sampled_Q_3D_;
  bool ok;
  int64 t0 = clock();
  for (int i = current_trial_; i < current_trial_ + n; ++i) {
    ok = TryOneBase();

    float fraction_try =
        static_cast<float>(i) / static_cast<float>(number_of_trials_);
    float fraction_time = static_cast<float>(clock() - t0) / 1000000.0 /
                          options_.max_time_seconds;
    float fraction = max(fraction_time, fraction_try);
    printf("done: %d%c best: %f                  \r",
           static_cast<int>(fraction * 100), '%', best_LCP_);
    fflush(stdout);
#ifdef TEST_RECORD_CONVERGENCE
    if (best_LCP_ > last_best_LCPTrace){
      last_best_LCPTrace = best_LCP_;

      // We are better than the last LCP. Update the matrix and transform Q.
      Point3D p(centroid_.x + centroid_Q_.x, centroid_.y + centroid_Q_.y,
                centroid_.z + centroid_Q_.z);

      cv::Mat t (4, 4, CV_64F, cv::Scalar(0.0));
      Transform(rotate_, cv::Point3f(0, 0, 0), cv::Point3f(0, 0, 0), &p);
      t.at<double>(0, 0) = rotate_.at<double>(0, 0);
      t.at<double>(0, 1) = rotate_.at<double>(0, 1);
      t.at<double>(0, 2) = rotate_.at<double>(0, 2);
      t.at<double>(0, 3) =
          centroid_.x - p.x + translate_.x + centroid_P_.x;
      t.at<double>(1, 0) = rotate_.at<double>(1, 0);
      t.at<double>(1, 1) = rotate_.at<double>(1, 1);
      t.at<double>(1, 2) = rotate_.at<double>(1, 2);
      t.at<double>(1, 3) =
          centroid_.y - p.y + translate_.y + centroid_P_.y;
      t.at<double>(2, 0) = rotate_.at<double>(2, 0);
      t.at<double>(2, 1) = rotate_.at<double>(2, 1);
      t.at<double>(2, 2) = rotate_.at<double>(2, 2);
      t.at<double>(2, 3) =
          centroid_.z - p.z + translate_.z + centroid_P_.z;
      t.at<double>(3, 0) = 0;
      t.at<double>(3, 1) = 0;
      t.at<double>(3, 2) = 0;
      t.at<double>(3, 3) = 1;

      myfile << i << " "
             << best_LCP_ << " "
             << (clock() - t0)/ 1000000.0 << " "
             // We suppose here that we look for the identity tranformation
             << std::abs(cv::sum(t)[0] - 4.) << " "
             << endl;
    }
#endif
    // ok means that we already have the desired LCP.
    if (ok || i > number_of_trials_ || fraction > 0.99) break;
  }

  current_trial_ += n;
  if (best_LCP_ > last_best_LCP) {
    // We are better than the last LCP. Update the matrix and transform Q.
    Point3D p(centroid_.x + centroid_Q_.x, centroid_.y + centroid_Q_.y,
              centroid_.z + centroid_Q_.z);
    Transform(rotate_, cv::Point3f(0, 0, 0), cv::Point3f(0, 0, 0), &p);
    *Q = Q_copy_;
    transformation->at<double>(0, 0) = rotate_.at<double>(0, 0);
    transformation->at<double>(0, 1) = rotate_.at<double>(0, 1);
    transformation->at<double>(0, 2) = rotate_.at<double>(0, 2);
    transformation->at<double>(0, 3) =
        centroid_.x - p.x + translate_.x + centroid_P_.x;
    transformation->at<double>(1, 0) = rotate_.at<double>(1, 0);
    transformation->at<double>(1, 1) = rotate_.at<double>(1, 1);
    transformation->at<double>(1, 2) = rotate_.at<double>(1, 2);
    transformation->at<double>(1, 3) =
        centroid_.y - p.y + translate_.y + centroid_P_.y;
    transformation->at<double>(2, 0) = rotate_.at<double>(2, 0);
    transformation->at<double>(2, 1) = rotate_.at<double>(2, 1);
    transformation->at<double>(2, 2) = rotate_.at<double>(2, 2);
    transformation->at<double>(2, 3) =
        centroid_.z - p.z + translate_.z + centroid_P_.z;
    transformation->at<double>(3, 0) = 0;
    transformation->at<double>(3, 1) = 0;
    transformation->at<double>(3, 2) = 0;
    transformation->at<double>(3, 3) = 1;

    // Transforms Q by the new transformation.
    for (int i = 0; i < Q->size(); ++i) {
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
#ifdef TEST_RECORD_CONVERGENCE
  myfile.close();
#endif

  return ok || current_trial_ >= number_of_trials_;
}

// The main 4PCS function. Computes the best rigid transformation and transfoms
// Q toward P by this transformation.
float Match4PCSImpl::ComputeTransformation(const std::vector<Point3D>& P,
                                           std::vector<Point3D>* Q,
                                           cv::Mat* transformation) {
  if (Q == nullptr || transformation == nullptr) return kLargeNumber;
  Initialize(P, *Q);
  *transformation = cv::Mat(4, 4, CV_64F, cv::Scalar(0.0));
  for (int i = 0; i < 4; ++i) transformation->at<double>(i, i) = 1.0;
  Perform_N_steps(number_of_trials_, transformation, Q);

  return best_LCP_;
}

Match4PCS::Match4PCS(const Match4PCSOptions& options)
    : pimpl_{new Match4PCSImpl{options}} {}

Match4PCS::~Match4PCS() {}

float Match4PCS::ComputeTransformation(const std::vector<Point3D>& P,
                                       std::vector<Point3D>* Q,
                                       cv::Mat* transformation) {
  return pimpl_->ComputeTransformation(P, Q, transformation);
}
}
