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

#include "super4pcs/algorithms/match4pcsBase.h"
#include "super4pcs/shared4pcs.h"
#include "super4pcs/sampling.h"
#include "super4pcs/accelerators/kdtree.h"

#include <vector>
#include <atomic>

#include <Eigen/Core>
#include <Eigen/Geometry>

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


namespace GlobalRegistration{

Match4PCSBase::Match4PCSBase(  const Match4PCSOptions& options
                             , const Utils::Logger& logger
#ifdef SUPER4PCS_USE_OPENMP
                             , const int omp_nthread_congruent
#endif
                               )
  :number_of_trials_(0)
  , max_base_diameter_(-1)
  , P_mean_distance_(1.0)
  , best_LCP_(0.0)
  , options_(options)
  , randomGenerator_(options.randomSeed)
  , logger_(logger)
#ifdef SUPER4PCS_USE_OPENMP
  , omp_nthread_congruent_(omp_nthread_congruent)
#endif
{
  base_3D_.resize(4);
}

Match4PCSBase::~Match4PCSBase(){}

Match4PCSBase::Scalar
Match4PCSBase::MeanDistance() {
  const Scalar kDiameterFraction = 0.2;
  using RangeQuery = GlobalRegistration::KdTree<Scalar>::RangeQuery<>;

  int number_of_samples = 0;
  Scalar distance = 0.0;

  for (size_t i = 0; i < sampled_P_3D_.size(); ++i) {

    RangeQuery query;
    query.sqdist = P_diameter_ * kDiameterFraction;
    query.queryPoint = sampled_P_3D_[i].pos().cast<Scalar>();

    auto result =
        kd_tree_.doQueryRestrictedClosestIndex(query , i);

    if (result.first != GlobalRegistration::KdTree<Scalar>::invalidIndex()) {
      distance += std::sqrt(result.second);
      number_of_samples++;
    }
  }

  return distance / number_of_samples;
}


bool Match4PCSBase::SelectRandomTriangle(int &base1, int &base2, int &base3) {
      int number_of_points = sampled_P_3D_.size();
      base1 = base2 = base3 = -1;

      // Pick the first point at random.
      int first_point = randomGenerator_() % number_of_points;

      const Scalar sq_max_base_diameter_ = max_base_diameter_*max_base_diameter_;

      // Try fixed number of times retaining the best other two.
      Scalar best_wide = 0.0;
      for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
        // Pick and compute
        const int second_point = randomGenerator_() % number_of_points;
        const int third_point = randomGenerator_() % number_of_points;
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
  kd_tree_ = GlobalRegistration::KdTree<Scalar>(number_of_points);

  for (size_t i = 0; i < number_of_points; ++i) {
    kd_tree_.add(sampled_P_3D_[i].pos());
  }
  kd_tree_.finalize();
}

bool Match4PCSBase::ComputeRigidTransformation(
        const std::array<Point3D, 4>& ref,
        const std::array<Point3D, 4>& candidate,
        const Eigen::Matrix<Scalar, 3, 1>& centroid1,
        Eigen::Matrix<Scalar, 3, 1> centroid2,
        Scalar max_angle,
        Eigen::Ref<MatrixType> transform,
        Scalar& rms_,
        bool computeScale ) const {

  rms_ = kLargeNumber;

  Scalar kSmallNumber = 1e-6;

  // We only use the first 3 pairs. This simplifies the process considerably
  // because it is the planar case.

  const VectorType& p0 = ref[0].pos();
  const VectorType& p1 = ref[1].pos();
  const VectorType& p2 = ref[2].pos();
        VectorType  q0 = candidate[0].pos();
        VectorType  q1 = candidate[1].pos();
        VectorType  q2 = candidate[2].pos();

  Scalar scaleEst (1.);

  // Compute scale factor if needed
  if (computeScale){
      const VectorType& p3 = ref[3].pos();
      const VectorType& q3 = candidate[3].pos();

      const Scalar ratio1 = (p1 - p0).norm() / (q1 - q0).norm();
      const Scalar ratio2 = (p3 - p2).norm() / (q3 - q2).norm();

      const Scalar ratioDev  = std::abs(ratio1/ratio2 - Scalar(1.));  // deviation between the two
      const Scalar ratioMean = (ratio1+ratio2)/Scalar(2.);            // mean of the two

      if ( ratioDev > Scalar(0.1) )
          return kLargeNumber;


      //Log<LogLevel::Verbose>( ratio1, " ", ratio2, " ", ratioDev, " ", ratioMean);
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


  // Discard singular solutions. The rotation should be orthogonal.
  if (((rotation * rotation).diagonal().array() - Scalar(1) > kSmallNumber).any())
      return false;

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
          first = scaleEst*candidate[i].pos() - centroid2;
          transformed = rotation * first;
          rms_ += (transformed - ref[i].pos() + centroid1).norm();
      }
  }

  rms_ /= Scalar(ref.size());

  Eigen::Transform<Scalar, 3, Eigen::Affine> etrans (Eigen::Transform<Scalar, 3, Eigen::Affine>::Identity());
  transform = etrans
      .scale(scaleEst)
      .translate(centroid1)
      .rotate(rotation)
      .translate(-centroid2)
      .matrix();

  return true;
}



// Verify a given transformation by computing the number of points in P at
// distance at most (normalized) delta from some point in Q. In the paper
// we describe randomized verification. We apply deterministic one here with
// early termination. It was found to be fast in practice.
Match4PCSBase::Scalar
Match4PCSBase::Verify(const Eigen::Ref<const MatrixType> &mat) const {
  using RangeQuery = GlobalRegistration::KdTree<Scalar>::RangeQuery<>;

#ifdef TEST_GLOBAL_TIMINGS
    Timer t_verify (true);
#endif

  // We allow factor 2 scaling in the normalization.
  const Scalar epsilon = options_.delta;
#ifdef SUPER4PCS_USE_WEIGHTED_LCP
  std::atomic<float> good_points(0);

  auto kernel = [](Scalar x) {
    return std::pow(std::pow(x,4) - Scalar(1), 2);
  };

  auto computeWeight = [kernel](Scalar sqx, Scalar th) {
    return kernel( std::sqrt(sqx) / th );
  };
#else
  std::atomic_uint good_points(0);
#endif
  const size_t number_of_points = sampled_Q_3D_.size();
  const size_t terminate_value = best_LCP_ * number_of_points;

  const Scalar sq_eps = epsilon*epsilon;
#ifdef SUPER4PCS_USE_WEIGHTED_LCP
  const Scalar    eps = std::sqrt(sq_eps);
#endif

  for (size_t i = 0; i < number_of_points; ++i) {

    // Use the kdtree to get the nearest neighbor
#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif

    RangeQuery query;
    query.queryPoint = (mat * sampled_Q_3D_[i].pos().homogeneous()).head<3>();
    query.sqdist     = sq_eps;

    auto result = kd_tree_.doQueryRestrictedClosestIndex( query );

#ifdef TEST_GLOBAL_TIMINGS
    kdTreeTime += Scalar(t.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif

    if ( result.first != GlobalRegistration::KdTree<Scalar>::invalidIndex() ) {
//      Point3D& q = sampled_P_3D_[near_neighbor_index[0]];
//      bool rgb_good =
//          (p.rgb()[0] >= 0 && q.rgb()[0] >= 0)
//              ? cv::norm(p.rgb() - q.rgb()) < options_.max_color_distance
//              : true;
//      bool norm_good = norm(p.normal()) > 0 && norm(q.normal()) > 0
//                           ? fabs(p.normal().ddot(q.normal())) >= cos_dist
//                           : true;
//      if (rgb_good && norm_good) {
#ifdef SUPER4PCS_USE_WEIGHTED_LCP
        assert (result.second <= query.sqdist);
        good_points = good_points + computeWeight(result.second, eps);
#else
        good_points++;
#endif
//      }
    }

    // We can terminate if there is no longer chance to get better than the
    // current best LCP.
    if (number_of_points - i + size_t(good_points) < terminate_value) {
      break;
    }
  }

#ifdef TEST_GLOBAL_TIMINGS
  verifyTime += Scalar(t_verify.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif
  return Scalar(good_points) / Scalar(number_of_points);
}

} // namespace Super4PCS

