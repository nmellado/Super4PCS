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
#include "shared4pcs.h"

#include <vector>
#include <chrono>

#include "Eigen/Core"
#include "Eigen/Geometry"                 // MatrixBase.homogeneous()
#include "Eigen/SVD"                      // Transform.computeRotationScaling()

#include "shared4pcs.h"
#include "sampling.h"
#include "accelerators/kdtree.h"

const double pi = std::acos(-1);



// Compute the closest points between two 3D line segments and obtain the two
// invariants corresponding to the closet points. This is the "intersection"
// point that determines the invariants. Since the 4 points are not exactly
// planar, we use the center of the line segment connecting the two closest
// points as the "intersection".
template < typename VectorType, typename Scalar>
static Scalar
distSegmentToSegment(const VectorType& p1, const VectorType& p2,
                     const VectorType& q1, const VectorType& q2,
                     Scalar& invariant1, Scalar& invariant2) {

  static const Scalar kSmallNumber = 0.0001;
  VectorType u = p2 - p1;
  VectorType v = q2 - q1;
  VectorType w = p1 - q1;
  Scalar a = u.dot(u);
  Scalar b = u.dot(v);
  Scalar c = v.dot(v);
  Scalar d = u.dot(w);
  Scalar e = v.dot(w);
  Scalar f = a * c - b * b;
  // s1,s2 and t1,t2 are the parametric representation of the intersection.
  // they will be the invariants at the end of this simple computation.
  Scalar s1 = 0.0;
  Scalar s2 = f;
  Scalar t1 = 0.0;
  Scalar t2 = f;

  if (f < kSmallNumber) {
    s1 = 0.0;
    s2 = 1.0;
    t1 = e;
    t2 = c;
  } else {
    s1 = (b * e - c * d);
    t1 = (a * e - b * d);
    if (s1 < 0.0) {
      s1 = 0.0;
      t1 = e;
      t2 = c;
    } else if (s1 > s2) {
      s1 = s2;
      t1 = e + b;
      t2 = c;
    }
  }

  if (t1 < 0.0) {
    t1 = 0.0;
    if (-d < 0.0)
      s1 = 0.0;
    else if (-d > a)
      s1 = s2;
    else {
      s1 = -d;
      s2 = a;
    }
  } else if (t1 > t2) {
    t1 = t2;
    if ((-d + b) < 0.0)
      s1 = 0;
    else if ((-d + b) > a)
      s1 = s2;
    else {
      s1 = (-d + b);
      s2 = a;
    }
  }
  invariant1 = (std::abs(s1) < kSmallNumber ? 0.0 : s1 / s2);
  invariant2 = (std::abs(t1) < kSmallNumber ? 0.0 : t1 / t2);

  return ( w + (invariant1 * u) - (invariant2 * v)).norm();
}


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
    query_point = sampled_P_3D_[i].pos().cast<Scalar>() ;

    Super4PCS::KdTree<Scalar>::Index resId =
        kd_tree_.doQueryRestrictedClosestIndex(query_point, P_diameter_ * kDiameterFraction, i);

    if (resId != Super4PCS::KdTree<Scalar>::invalidIndex()) {
      distance += (sampled_P_3D_[i].pos() - sampled_P_3D_[resId].pos()).norm();
      number_of_samples++;
    }
  }

  return distance / number_of_samples;
}

void Match4PCSBase::init(const std::vector<Point3D>& P,
                         const std::vector<Point3D>& Q){

#ifdef TEST_GLOBAL_TIMINGS
    kdTreeTime = 0;
    totalTime  = 0;
    verifyTime = 0;
#endif

    const Scalar kSmallError = 0.00001;
    const int kMinNumberOfTrials = 4;
    const Scalar kDiameterFraction = 0.3;

    centroid_P_ = VectorType::Zero();
    centroid_Q_ = VectorType::Zero();

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
        centroid_P_ += sampled_P_3D_[i].pos();
    }

    centroid_P_ /= Scalar(sampled_P_3D_.size());

    for (int i = 0; i < sampled_Q_3D_.size(); ++i) {
        centroid_Q_ += sampled_Q_3D_[i].pos();
    }

    centroid_Q_ /= Scalar(sampled_Q_3D_.size());

    // Move the samples to the centroids to allow robustness in rotation.
    for (int i = 0; i < sampled_P_3D_.size(); ++i) {
        sampled_P_3D_[i].pos() -= centroid_P_;
    }
    for (int i = 0; i < sampled_Q_3D_.size(); ++i) {
        sampled_Q_3D_[i].pos() -= centroid_Q_;
    }


    initKdTree();
    // Compute the diameter of P approximately (randomly). This is far from being
    // Guaranteed close to the diameter but gives good results for most common
    // objects if they are densely sampled.
    P_diameter_ = 0.0;
    for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
        int at = rand() % sampled_Q_3D_.size();
        int bt = rand() % sampled_Q_3D_.size();

        Scalar l = (sampled_Q_3D_[bt].pos() - sampled_Q_3D_[at].pos()).norm();
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

    // call Virtual handler
    Initialize(P, Q);

    best_LCP_ = Verify(transform_);
    printf("Initial LCP: %f\n", best_LCP_);
}


bool Match4PCSBase::SelectRandomTriangle(int &base1, int &base2, int &base3) {
      int number_of_points = sampled_P_3D_.size();
      base1 = base2 = base3 = -1;

      // Pick the first point at random.
      int first_point = rand() % number_of_points;

      const Scalar sq_max_base_diameter_ = max_base_diameter_*max_base_diameter_;

      // Try fixed number of times retaining the best other two.
      Scalar best_wide = 0.0;
      for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
        // Pick and compute
        const int second_point = rand() % number_of_points;
        const int third_point = rand() % number_of_points;
        const VectorType u =
                sampled_P_3D_[second_point].pos() -
                sampled_P_3D_[first_point].pos();
        const VectorType w =
                sampled_P_3D_[third_point].pos() -
                sampled_P_3D_[first_point].pos();
        // We try to have wide triangles but still not too large.
        Scalar how_wide = (u.cross(w)).norm();
        if (how_wide > best_wide &&
                u.squaredNorm() < sq_max_base_diameter_ &&
                w.squaredNorm() < sq_max_base_diameter_) {
          best_wide = how_wide;
          base1 = first_point;
          base2 = second_point;
          base3 = third_point;
        }
      }
      return base1 != -1 && base2 != -1 && base3 != -1;
}



// Try the current base in P and obtain the best pairing, i.e. the one that
// gives the smaller distance between the two closest points. The invariants
// corresponding the the base pairing are computed.
bool Match4PCSBase::TryQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                     int& id1, int& id2, int& id3, int& id4) {

  Scalar min_distance = std::numeric_limits<Scalar>::max();
  int best1, best2, best3, best4;
  best1 = best2 = best3 = best4 = -1;
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
      Scalar segment_distance = distSegmentToSegment(
                  base_3D_[i].pos(), base_3D_[j].pos(),
                  base_3D_[k].pos(), base_3D_[l].pos(),
                  local_invariant1, local_invariant2);
      // Retail the smallest distance and the best order so far.
      if (segment_distance < min_distance) {
        min_distance = segment_distance;
        best1 = i;
        best2 = j;
        best3 = k;
        best4 = l;
        invariant1 = local_invariant1;
        invariant2 = local_invariant2;
      }
    }
  }

  if(best1 < 0 || best2 < 0 || best3 < 0 || best4 < 0 ) return false;

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
bool Match4PCSBase::SelectQuadrilateral(Scalar& invariant1, Scalar& invariant2,
                                        int& base1, int& base2, int& base3,
                                        int& base4) {

  const Scalar kBaseTooSmall (0.2);
  int current_trial = 0;

  // Try fix number of times.
  while (current_trial < kNumberOfDiameterTrials) {
    // Select a triangle if possible. otherwise fail.
    if (!SelectRandomTriangle(base1, base2, base3)){
      return false;
    }

    base_3D_[0] = sampled_P_3D_[base1];
    base_3D_[1] = sampled_P_3D_[base2];
    base_3D_[2] = sampled_P_3D_[base3];

    // The 4th point will be a one that is close to be planar to the other 3
    // while still not too close to them.
    const double x1 = base_3D_[0].x();
    const double y1 = base_3D_[0].y();
    const double z1 = base_3D_[0].z();
    const double x2 = base_3D_[1].x();
    const double y2 = base_3D_[1].y();
    const double z2 = base_3D_[1].z();
    const double x3 = base_3D_[2].x();
    const double y3 = base_3D_[2].y();
    const double z3 = base_3D_[2].z();

    // Fit a plan.
    Scalar denom = (-x3 * y2 * z1 + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 -
                    x2 * y1 * z3 + x1 * y2 * z3);

    if (denom != 0) {
      Scalar A =
          (-y2 * z1 + y3 * z1 + y1 * z2 - y3 * z2 - y1 * z3 + y2 * z3) / denom;
      Scalar B =
          (x2 * z1 - x3 * z1 - x1 * z2 + x3 * z2 + x1 * z3 - x2 * z3) / denom;
      Scalar C =
          (-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3) / denom;
      base4 = -1;
      Scalar best_distance = std::numeric_limits<Scalar>::max();
      // Go over all points in P.
      const Scalar too_small = std::pow(max_base_diameter_ * kBaseTooSmall, 2);
      for (unsigned int i = 0; i < sampled_P_3D_.size(); ++i) {
        if ((sampled_P_3D_[i].pos()- sampled_P_3D_[base1].pos()).squaredNorm() >= too_small &&
            (sampled_P_3D_[i].pos()- sampled_P_3D_[base2].pos()).squaredNorm() >= too_small &&
            (sampled_P_3D_[i].pos()- sampled_P_3D_[base3].pos()).squaredNorm() >= too_small) {
          // Not too close to any of the first 3.
          const Scalar distance =
              std::abs(A * sampled_P_3D_[i].x() + B * sampled_P_3D_[i].y() +
                   C * sampled_P_3D_[i].z() - 1.0);
          // Search for the most planar.
          if (distance < best_distance) {
            best_distance = distance;
            base4 = int(i);
          }
        }
      }
      // If we have a good one we can quit.
      if (base4 != -1) {
        base_3D_[3] = sampled_P_3D_[base4];
        if(TryQuadrilateral(invariant1, invariant2, base1, base2, base3, base4))
            return true;
      }
    }
    current_trial++;
  }

  // We failed to find good enough base..
  return false;
}

void Match4PCSBase::initKdTree(){
  size_t number_of_points = sampled_P_3D_.size();

  // Build the kdtree.
  kd_tree_ = Super4PCS::KdTree<Scalar>(number_of_points);

  for (size_t i = 0; i < number_of_points; ++i) {
    kd_tree_.add(sampled_P_3D_[i].pos());
  }
  kd_tree_.finalize();
}

bool Match4PCSBase::TryCongruentSet(
        int base_id1,
        int base_id2,
        int base_id3,
        int base_id4,
        const std::vector<match_4pcs::Quadrilateral>& congruent_quads,
        size_t &nbCongruent){
    std::array<std::pair<Point3D, Point3D>,4> congruent_points;

    // get references to the basis coordinates
    const Point3D& b1 = sampled_P_3D_[base_id1];
    const Point3D& b2 = sampled_P_3D_[base_id2];
    const Point3D& b3 = sampled_P_3D_[base_id3];
    const Point3D& b4 = sampled_P_3D_[base_id4];


    // Centroid of the basis, computed once and using only the three first points
    Eigen::Matrix<Scalar, 3, 1> centroid1 = (b1.pos() + b2.pos() + b3.pos()) / Scalar(3);

    // Centroid of the sets, computed in the loop using only the three first points
    Eigen::Matrix<Scalar, 3, 1> centroid2;

    // set the basis coordinates in the congruent quad array
    congruent_points[0].first = b1;
    congruent_points[1].first = b2;
    congruent_points[2].first = b3;
    congruent_points[3].first = b4;

    nbCongruent = 0;

    Eigen::Matrix<Scalar, 4, 4> transform;
    for (size_t i = 0; i < congruent_quads.size(); ++i) {
      const int a = congruent_quads[i].vertices[0];
      const int b = congruent_quads[i].vertices[1];
      const int c = congruent_quads[i].vertices[2];
      const int d = congruent_quads[i].vertices[3];
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

      centroid2 = (congruent_points[0].second.pos() +
                   congruent_points[1].second.pos() +
                   congruent_points[2].second.pos()) / Scalar(3.);

      Scalar rms = -1;

      const bool ok =
      ComputeRigidTransformation(congruent_points,   // input congruent quads
                                 centroid1,          // input: basis centroid
                                 centroid2,          // input: candidate quad centroid
                                 options_.max_angle * pi / 180.0, // maximum per-dimension angle, check return value to detect invalid cases
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

          nbCongruent++;
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


bool Match4PCSBase::ComputeRigidTransformation(const std::array< std::pair<Point3D, Point3D>,4>& pairs,
        const Eigen::Matrix<Scalar, 3, 1>& centroid1,
        Eigen::Matrix<Scalar, 3, 1> centroid2,
        Scalar max_angle,
        Eigen::Ref<MatrixType> transform,
        Scalar& rms_,
        bool computeScale ) {

  rms_ = kLargeNumber;

  if (pairs.size() == 0 || pairs.size() % 2 != 0)
      return false;


  Scalar kSmallNumber = 1e-6;

  // We only use the first 3 pairs. This simplifies the process considerably
  // because it is the planar case.

  const VectorType& p0 = pairs[0].first.pos();
  const VectorType& p1 = pairs[1].first.pos();
  const VectorType& p2 = pairs[2].first.pos();
        VectorType  q0 = pairs[0].second.pos();
        VectorType  q1 = pairs[1].second.pos();
        VectorType  q2 = pairs[2].second.pos();

  Scalar scaleEst (1.);

  // Compute scale factor if needed
  if (computeScale){
      const VectorType& p3 = pairs[3].first.pos();
      const VectorType& q3 = pairs[3].second.pos();

      const Scalar ratio1 = (p1 - p0).norm() / (q1 - q0).norm();
      const Scalar ratio2 = (p3 - p2).norm() / (q3 - q2).norm();

      const Scalar ratioDev  = std::abs(ratio1/ratio2 - Scalar(1.));  // deviation between the two
      const Scalar ratioMean = (ratio1+ratio2)/Scalar(2.);            // mean of the two

      if ( ratioDev > Scalar(0.1) )
          return kLargeNumber;

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

  VectorType vector_p1 = p1 - p0;
  if (vector_p1.squaredNorm() == 0) return kLargeNumber;
  vector_p1.normalize();
  VectorType vector_p2 = (p2 - p0) - ((p2 - p0).dot(vector_p1)) * vector_p1;
  if (vector_p2.squaredNorm() == 0) return kLargeNumber;
  vector_p2.normalize();
  VectorType vector_p3 = vector_p1.cross(vector_p2);

  VectorType vector_q1 = q1 - q0;
  if (vector_q1.squaredNorm() == 0) return kLargeNumber;
  vector_q1.normalize();
  VectorType vector_q2 = (q2 - q0) - ((q2 - q0).dot(vector_q1)) * vector_q1;
  if (vector_q2.squaredNorm() == 0) return kLargeNumber;
  vector_q2.normalize();
  VectorType vector_q3 = vector_q1.cross(vector_q2);

  //cv::Mat rotation = cv::Mat::eye(3, 3, CV_64F);
  Eigen::Matrix<Scalar, 3, 3> rotation = Eigen::Matrix<Scalar, 3, 3>::Identity();

  Eigen::Matrix<Scalar, 3, 3> rotate_p;
  rotate_p.row(0) = vector_p1;
  rotate_p.row(1) = vector_p2;
  rotate_p.row(2) = vector_p3;

  Eigen::Matrix<Scalar, 3, 3> rotate_q;
  rotate_q.row(0) = vector_q1;
  rotate_q.row(1) = vector_q2;
  rotate_q.row(2) = vector_q3;

  rotation = rotate_p.transpose() * rotate_q;

  /*cv::Mat rotate_p(3, 3, CV_64F);
  rotate_p.at<double>(0, 0) = vector_p1.x();
  rotate_p.at<double>(0, 1) = vector_p1.y();
  rotate_p.at<double>(0, 2) = vector_p1.z();
  rotate_p.at<double>(1, 0) = vector_p2.x();
  rotate_p.at<double>(1, 1) = vector_p2.y();
  rotate_p.at<double>(1, 2) = vector_p2.z();
  rotate_p.at<double>(2, 0) = vector_p3.x();
  rotate_p.at<double>(2, 1) = vector_p3.y();
  rotate_p.at<double>(2, 2) = vector_p3.z();*/

  /*cv::Mat rotate_q(3, 3, CV_64F);
  rotate_q.at<double>(0, 0) = vector_q1.x();
  rotate_q.at<double>(0, 1) = vector_q1.y();
  rotate_q.at<double>(0, 2) = vector_q1.z();
  rotate_q.at<double>(1, 0) = vector_q2.x();
  rotate_q.at<double>(1, 1) = vector_q2.y();
  rotate_q.at<double>(1, 2) = vector_q2.z();
  rotate_q.at<double>(2, 0) = vector_q3.x();
  rotate_q.at<double>(2, 1) = vector_q3.y();
  rotate_q.at<double>(2, 2) = vector_q3.z();*/

  //rotation = rotate_p.t() * rotate_q;

  // Discard singular solutions. The rotation should be orthogonal.
  if (((rotation * rotation).diagonal().array() - Scalar(1) > kSmallNumber).any())
      return false;


  /*cv::Mat unit = rotation * rotation.t();
  if (std::abs(unit.at<double>(0, 0) - 1.0) > kSmallNumber ||
      std::abs(unit.at<double>(1, 1) - 1.0) > kSmallNumber ||
          std::abs(unit.at<double>(2, 2) - 1.0) > kSmallNumber){
      return false;
  }*/

  //FIXME
  if (max_angle >= 0) {
      // Discard too large solutions (todo: lazy evaluation during boolean computation
      if (! (
                  std::abs(std::atan2(rotation(2, 1), rotation(2, 2)))
                  <= max_angle &&

                  std::abs(std::atan2(-rotation(2, 0),
                                      std::sqrt(std::pow(rotation(2, 1),2) +
                                                std::pow(rotation(2, 2),2))))
                  <= max_angle &&

                  std::abs(atan2(rotation(1, 0), rotation(0, 0)))
                  <= max_angle
             ))
          return false;
  }


  //FIXME
  // Compute rms and return it.
  rms_ = Scalar(0.0);
  {
      VectorType first, transformed;

      //cv::Mat first(3, 1, CV_64F), transformed;
      for (int i = 0; i < 3; ++i) {
          first = scaleEst*pairs[i].second.pos() - centroid2;
          transformed = rotation * first;

//          first.at<double>(0, 0) = scaleEst*pairs[i].second.x() - centroid2(0);
//          first.at<double>(1, 0) = scaleEst*pairs[i].second.y() - centroid2(1);
//          first.at<double>(2, 0) = scaleEst*pairs[i].second.z() - centroid2(2);
//          transformed = rotation * first;
          rms_ += (transformed - pairs[i].first.pos() + centroid1).norm();
//          rms_ += sqrt(std::pow(transformed.at<double>(0, 0) -
//                              (pairs[i].first.x() - centroid1(0)), 2) +
//                       std::pow(transformed.at<double>(1, 0) -
//                              (pairs[i].first.y() - centroid1(1)), 2) +
//                       std::pow(transformed.at<double>(2, 0) -
//                              (pairs[i].first.z() - centroid1(2)), 2));
      }
  }

  rms_ /= Scalar(pairs.size());

  Eigen::Transform<Scalar, 3, Eigen::Affine> etrans (Eigen::Transform<Scalar, 3, Eigen::Affine>::Identity());

  // compute rotation and translation
  {
      //std::cout << scaleEst << endl;

      etrans.scale(scaleEst);       // apply scale factor
      etrans.translate(centroid1);  // translation between quads
      etrans.rotate(rotation);           // rotate to align frames
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
Match4PCSBase::Verify(const Eigen::Ref<const MatrixType> &mat) {

#ifdef TEST_GLOBAL_TIMINGS
    Timer t_verify (true);
#endif

  // We allow factor 2 scaling in the normalization.
  const Scalar epsilon = options_.delta;
  int good_points = 0;
  const size_t number_of_points = sampled_Q_3D_.size();
  const int terminate_value = best_LCP_ * number_of_points;

  const Scalar sq_eps = epsilon*epsilon;

  for (int i = 0; i < number_of_points; ++i) {

    // Use the kdtree to get the nearest neighbor
#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif

    Super4PCS::KdTree<Scalar>::Index resId =
    kd_tree_.doQueryRestrictedClosestIndex(
                (mat * sampled_Q_3D_[i].pos().homogeneous()).head<3>(),
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
  return Scalar(good_points) / Scalar(number_of_points);
}



// The main 4PCS function. Computes the best rigid transformation and transfoms
// Q toward P by this transformation.
Match4PCSBase::Scalar
Match4PCSBase::ComputeTransformation(const std::vector<Point3D>& P,
                                     std::vector<Point3D>* Q,
                                     Eigen::Ref<MatrixType> transformation) {

  if (Q == nullptr) return kLargeNumber;
  init(P, *Q);

  transformation = MatrixType::Identity();
  Perform_N_steps(number_of_trials_, transformation, Q);

#ifdef TEST_GLOBAL_TIMINGS
  std::cout << "----------- Timings (msec) -------------"          << std::endl;
  std::cout << " Total computation time  : " << totalTime          << std::endl;
  std::cout << " Total verify time       : " << verifyTime         << std::endl;
  std::cout << "    Kdtree query         : " << kdTreeTime         << std::endl;
#endif

  return best_LCP_;
}


// Performs N RANSAC iterations and compute the best transformation. Also,
// transforms the set Q by this optimal transformation.
bool Match4PCSBase::Perform_N_steps(int n,
                                    Eigen::Ref<MatrixType> transformation,
                                    std::vector<Point3D>* Q) {
	using std::chrono::system_clock;
  if (Q == nullptr) return false;

#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif

  Scalar last_best_LCP = best_LCP_;
  bool ok = false;
  std::chrono::time_point<system_clock> t0 = system_clock::now(), end;
  for (int i = current_trial_; i < current_trial_ + n; ++i) {
    ok = TryOneBase();

    Scalar fraction_try  = Scalar(i) / Scalar(number_of_trials_);
    Scalar fraction_time = 
		std::chrono::duration_cast<std::chrono::seconds>
		(system_clock::now() - t0).count() /
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
        centroid_P = centroid_P_;
        centroid_Q = centroid_Q_;

        Eigen::Matrix<Scalar, 3, 3> rot, scale;
        Eigen::Transform<Scalar, 3, Eigen::Affine> (transform_).computeRotationScaling(&rot, &scale);
        transform_.col(3) = (qcentroid1_ + centroid_P - ( rot * scale * (qcentroid2_ + centroid_Q))).homogeneous();
        transformation = transform_;
    }

    // Transforms Q by the new transformation.
    for (size_t i = 0; i < Q->size(); ++i) {
      (*Q)[i].pos() = (transformation * (*Q)[i].pos().homogeneous()).head<3>();

//      cv::Mat first(4, 1, CV_64F), transformed;
//      first.at<double>(0, 0) = (*Q)[i].x();
//      first.at<double>(1, 0) = (*Q)[i].y();
//      first.at<double>(2, 0) = (*Q)[i].z();
//      first.at<double>(3, 0) = 1;
//      transformed = *transformation * first;
//      (*Q)[i].x() = transformed.at<double>(0, 0);
//      (*Q)[i].y() = transformed.at<double>(1, 0);
//      (*Q)[i].z() = transformed.at<double>(2, 0);
    }
  }
#ifdef TEST_GLOBAL_TIMINGS
    totalTime += Scalar(t.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif

  return ok || current_trial_ >= number_of_trials_;
}

// Pick one base, finds congruent 4-points in Q, verifies for all
// transformations, and retains the best transformation and LCP. This is
// a complete RANSAC iteration.
bool Match4PCSBase::TryOneBase() {
  Scalar invariant1, invariant2;
  int base_id1, base_id2, base_id3, base_id4;

//#define STATIC_BASE

#ifdef STATIC_BASE
  static bool first_time = true;

  if (first_time){
      base_id1 = 0;
      base_id2 = 3;
      base_id3 = 1;
      base_id4 = 4;

      base_3D_[0] = sampled_P_3D_ [base_id1];
      base_3D_[1] = sampled_P_3D_ [base_id2];
      base_3D_[2] = sampled_P_3D_ [base_id3];
      base_3D_[3] = sampled_P_3D_ [base_id4];

      TryQuadrilateral(&invariant1, &invariant2, base_id1, base_id2, base_id3, base_id4);

      first_time = false;
  }
  else
      return false;

#else

  if (!SelectQuadrilateral(invariant1, invariant2, base_id1, base_id2,
                           base_id3, base_id4)) {
    return false;
  }
#endif

  // Computes distance between pairs.
  const Scalar distance1 = (base_3D_[0].pos()- base_3D_[1].pos()).norm();
  const Scalar distance2 = (base_3D_[2].pos()- base_3D_[3].pos()).norm();

  std::vector<std::pair<int, int>> pairs1, pairs2;
  std::vector<match_4pcs::Quadrilateral> congruent_quads;

  // Compute normal angles.
  const Scalar normal_angle1 = (base_3D_[0].normal() - base_3D_[1].normal()).norm();
  const Scalar normal_angle2 = (base_3D_[2].normal() - base_3D_[3].normal()).norm();

  ExtractPairs(distance1, normal_angle1, distance_factor * options_.delta, 0,
                  1, &pairs1);
  ExtractPairs(distance2, normal_angle2, distance_factor * options_.delta, 2,
                  3, &pairs2);

//  std::cout << "Pair creation output: \n"
//            << pairs1.size() << " - "
//            << pairs2.size() << std::endl;

  if (pairs1.size() == 0 || pairs2.size() == 0) {
    return false;
  }


  if (!FindCongruentQuadrilaterals(invariant1, invariant2,
                                   distance_factor * options_.delta,
                                   distance_factor * options_.delta,
                                   pairs1,
                                   pairs2,
                                   &congruent_quads)) {
    return false;
  }

  size_t nb = 0;

  bool match = TryCongruentSet(base_id1, base_id2, base_id3, base_id4,
                               congruent_quads,
                               nb);

  //if (nb != 0)
  //    std::cout << "Congruent quads: (" << nb << ")    " << std::endl;

  return match;
}

} // namespace Super4PCS

