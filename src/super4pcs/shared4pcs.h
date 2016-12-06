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
class Point3D : public cv::Point3f {
 public:
  inline Point3D(double x, double y, double z) : cv::Point3f(x, y, z) {}
  inline Point3D(const cv::Point3f& other): cv::Point3f(other) {}
  inline Point3D(const Point3D& other):
      cv::Point3f(other),
      normal_(other.normal_),
      rgb_(other.rgb_)/*,
      hsv_(other.hsv_)*/ {}
  template<typename Scalar>
  explicit inline Point3D(const Eigen::Matrix<Scalar, 3, 1>& other):
      cv::Point3f(other(0), other(1), other(2)){
  }

  inline Point3D() : cv::Point3f(0.0f, 0.0f, 0.0f) {}
  inline const cv::Vec3f& rgb() const { return rgb_; }
//  inline const cv::Vec3f& hsv() const { return hsv_; }
  inline const cv::Point3f& normal() const { return normal_; }
  inline void set_rgb(const cv::Vec3f& rgb) {
      rgb_ = rgb;
      hasColor_ = true;
//      internal::RGB2HSV(rgb_[0]/255.f,rgb_[1]/255.f,rgb_[2]/255.f, hsv_[0],hsv_[1],hsv_[2]);
  }
  inline void set_normal(const cv::Point3d& normal) {
      double norm = cv::norm(normal);
      normal_.x = normal.x / norm;
      normal_.y = normal.y / norm;
      normal_.z = normal.z / norm;
  }
  inline void set_normal(const cv::Point3f& normal) {
      float norm = cv::norm(normal);
      normal_.x = normal.x / norm;
      normal_.y = normal.y / norm;
      normal_.z = normal.z / norm;
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
  cv::Point3f normal_{0.0f, 0.0f, 0.0f};
  // Color.
  cv::Vec3f rgb_{-1.0f, -1.0f, -1.0f};
//  cv::Vec3f hsv_{-1.0f, -1.0f, -1.0f};
  bool hasColor_ = false;
};

// Throughout this file the first model is called P, the second is Q
// and the measure we adopt for the quality of a given transformation is the
// Largest Common Pointset (LCP), normalized to the size of P thus lies in
// [0,1] (See the paper for details).

const float kLargeNumber = 1e9;

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

class MatchShared4PCSImpl;



class MatchSuper4PCSImpl;

// Class for the computation of the 4PCS algorithm.
class MatchSuper4PCS {
 public:
  explicit MatchSuper4PCS(const Match4PCSOptions& options);

  ~MatchSuper4PCS();

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
  std::unique_ptr<MatchSuper4PCSImpl> pimpl_;
};

} // namespace match_4pcs



#endif //_SHARED_4PCS_H_
