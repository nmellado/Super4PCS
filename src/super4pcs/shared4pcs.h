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

#ifndef _SHARED_4PCS_H_
#define _SHARED_4PCS_H_

#include <vector>
#include <memory>

#include <opencv2/core/core.hpp>
#include "Eigen/Core"

#include "utils/disablewarnings.h"

namespace match_4pcs {


static const float kLargeNumber = 1e9;

namespace internal{
inline void RGB2HSV(float r, float g, float b,
                    float &h, float &s, float &v)
{
    float K = 0.f;

    if (g < b)
    {
        std::swap(g, b);
        K = -1.f;
    }

    if (r < g)
    {
        std::swap(r, g);
        K = -2.f / 6.f - K;
    }

    float chroma = r - std::min(g, b);
    h = std::abs(K + (g - b) / (6.f * chroma + 1e-20f));
    s = chroma / (r + 1e-20f);
    v = r;
}
}

// The basic 3D point structure. A point potentially contains also directional
// information and color.
//template <typename _Scalar = float>
class Point3D : public cv::Point3d {
 public:
  using Scalar = double; //_Scalar;
  using VectorType = Eigen::Matrix<Scalar, 3, 1>;

  inline Point3D(double x, double y, double z) : cv::Point3d(x, y, z) {}
  inline Point3D(const cv::Point3d& other): cv::Point3d(other) {}
  inline Point3D(const Point3D& other):
      cv::Point3d(other),
      normal_(other.normal_),
      rgb_(other.rgb_)/*,
      hsv_(other.hsv_)*/ {}
  template<typename Scalar>
  explicit inline Point3D(const Eigen::Matrix<Scalar, 3, 1>& other):
      cv::Point3d(other(0), other(1), other(2)){
  }

  inline Point3D() : cv::Point3d(0.0f, 0.0f, 0.0f) {}
  inline const VectorType& rgb() const { return rgb_; }
//  inline const cv::Vec3f& hsv() const { return hsv_; }
  inline const VectorType& normal() const { return normal_; }
  inline void set_rgb(const VectorType& rgb) {
      rgb_ = rgb;
      hasColor_ = true;
//      internal::RGB2HSV(rgb_[0]/255.f,rgb_[1]/255.f,rgb_[2]/255.f, hsv_[0],hsv_[1],hsv_[2]);
  }
  inline void set_normal(const VectorType& normal) {
      normal_ = normal.normalized();
  }

  inline void normalize() {
    double n = cv::norm(*this);
    x /= n;
    y /= n;
    z /= n;
  }
  inline bool hasColor() const { return hasColor_; }

 private:
  // Normal.
  VectorType normal_{0.0f, 0.0f, 0.0f};
  // Color.
  VectorType rgb_{-1.0f, -1.0f, -1.0f};
//  cv::Vec3f hsv_{-1.0f, -1.0f, -1.0f};
  bool hasColor_ = false;
};

// Throughout this file the first model is called P, the second is Q
// and the measure we adopt for the quality of a given transformation is the
// Largest Common Pointset (LCP), normalized to the size of P thus lies in
// [0,1] (See the paper for details).


// Local static helpers that do not depend on state.

static inline float Square(float x) { return x * x; }

static inline float PointsDistance(const Point3D& p, const Point3D& q) {
  return cv::norm(p - q);
}

// Compute the closest points between two 3D line segments and obtain the two
// invariants corresponding to the closet points. This is the "intersection"
// point that determines the invariants. Since the 4 points are not exactly
// planar, we use the center of the line segment connecting the two closest
// points as the "intersection".
static float distSegmentToSegment(const cv::Point3d& p1, const cv::Point3d& p2,
                           const cv::Point3d& q1, const cv::Point3d& q2,
                           double* invariant1, double* invariant2) {
  if (invariant1 == 0 || invariant2 == 0) return kLargeNumber;
  const float kSmallNumber = 0.0001;
  cv::Point3d u = p2 - p1;
  cv::Point3d v = q2 - q1;
  cv::Point3d w = p1 - q1;
  float a = u.dot(u);
  float b = u.dot(v);
  float c = v.dot(v);
  float d = u.dot(w);
  float e = v.dot(w);
  float f = a * c - b * b;
  // s1,s2 and t1,t2 are the parametric representation of the intersection.
  // they will be the invariants at the end of this simple computation.
  float s1 = 0.0;
  float s2 = f;
  float t1 = 0.0;
  float t2 = f;

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
  *invariant1 = (abs(s1) < kSmallNumber ? 0.0 : s1 / s2);
  *invariant2 = (abs(t1) < kSmallNumber ? 0.0 : t1 / t2);
  cv::Point3d distance = w + (*invariant1 * u) - (*invariant2 * v);
  return cv::norm(distance);
}

// ----- 4PCS Options -----
// delta and overlap_estimation are the application parameters. All other
// parameters are more likely to keep fixed and they can be set via the setters.
struct Match4PCSOptions {
  Match4PCSOptions() {}

  // The delta for the LCP (see the paper).
  double delta = 5.0;
  // Estimated overlap between P and Q. This is the fraction of points in P that
  // may have corresponding point in Q. It's being used to estimate the number
  // of RANSAC iterations needed to guarantee small failure probability.
  double overlap_estimation = 0.2;

  // Maximum normal difference.
  double max_normal_difference = 10.0;
  // Maximum translation distance.
  double max_translation_distance = kLargeNumber;
  // Maximum rotation angle.
  double max_angle = kLargeNumber;
  // Maximum color RGB distance between corresponding vertices.
  double max_color_distance = 100.0;
  // Threshold on the value of the target function (LCP, see the paper).
  // It is used to terminate the process once we reached this value.
  double terminate_threshold = 1.0;
  // The number of points in the sample. We sample this number of points
  // uniformly from P and Q.
  int sample_size = 200;
  // Maximum time we allow the computation to take. This makes the algorithm
  // an ANY TIME algorithm that can be stopped at any time, producing the best
  // solution so far.
  int max_time_seconds = 60;
};

} // namespace match_4pcs



#endif //_SHARED_4PCS_H_
