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
// Authors: Dror Aiger, Yoni Weill, Nicolas Mellado
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

#include "super4pcs/utils/disablewarnings.h"

#include <Eigen/Core>

#include <vector>
#include <iostream>
#include <fstream>
#include <array>
#include <random>

namespace GlobalRegistration {

/// The basic 3D point structure. A point potentially contains also directional
/// information and color.
class Point3D {
 public:
  using Scalar = float; //_Scalar;
  using VectorType = Eigen::Matrix<Scalar, 3, 1>;

  inline Point3D(Scalar x, Scalar y, Scalar z) : pos_({ x, y, z}) {}
  inline Point3D(const Point3D& other):
      pos_(other.pos_),
      normal_(other.normal_),
      rgb_(other.rgb_) {}
  template<typename Scalar>
  explicit inline Point3D(const Eigen::Matrix<Scalar, 3, 1>& other):
      pos_({ other(0), other(1), other(2) }){
  }

  inline Point3D() {}
  inline VectorType& pos() { return pos_ ; }
  inline const VectorType& pos() const { return pos_ ; }
  inline const VectorType& rgb() const { return rgb_; }

  inline const VectorType& normal() const { return normal_; }
  inline void set_rgb(const VectorType& rgb) {
      rgb_ = rgb;
  }
  inline void set_normal(const VectorType& normal) {
      normal_ = normal.normalized();
  }

  inline void normalize() {
    pos_.normalize();
  }
  inline bool hasColor() const { return rgb_.squaredNorm() > Scalar(0.001); }

  Scalar& x() { return pos_.coeffRef(0); }
  Scalar& y() { return pos_.coeffRef(1); }
  Scalar& z() { return pos_.coeffRef(2); }

  Scalar x() const { return pos_.coeff(0); }
  Scalar y() const { return pos_.coeff(1); }
  Scalar z() const { return pos_.coeff(2); }



 private:
  /// Normal.
  VectorType pos_{0.0f, 0.0f, 0.0f};
  /// Normal.
  VectorType normal_{0.0f, 0.0f, 0.0f};
  /// Color.
  VectorType rgb_{-1.0f, -1.0f, -1.0f};
};



/// Holds a base from P. The base contains 4 points (indices) from the set P.
struct Quadrilateral {
    std::array <int, 4> vertices;
    inline Quadrilateral(int vertex0, int vertex1, int vertex2, int vertex3) {
        vertices = { vertex0, vertex1, vertex2, vertex3 };
    }

    inline bool operator< (const Quadrilateral& rhs) const {
        return    vertices[0] != rhs[0] ? vertices[0] < rhs[0]
                : vertices[1] != rhs[1] ? vertices[1] < rhs[1]
                : vertices[2] != rhs[2] ? vertices[2] < rhs[2]
                : vertices[3] < rhs[3];
    }

    inline bool operator== (const Quadrilateral& rhs) const {
        return  vertices[0] == rhs[0] &&
                vertices[1] == rhs[1] &&
                vertices[2] == rhs[2] &&
                vertices[3] == rhs[3];
    }

    int  operator[](int idx) const { return vertices[idx]; }
    int& operator[](int idx)       { return vertices[idx]; }
};

inline std::ofstream& operator<<(std::ofstream& ofs, const Quadrilateral& q){
    ofs << "[" << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << "]";
    return ofs;
}

// ----- 4PCS Options -----
/// delta and overlap_estimation are the application parameters. All other
/// parameters are more likely to keep fixed and they can be set via the setters.
struct Match4PCSOptions {
  using Scalar = typename Point3D::Scalar;
  Match4PCSOptions() {}

  /// The delta for the LCP (see the paper).
  Scalar delta = 5.0;

  /// Maximum normal difference.
  Scalar max_normal_difference = -1;
  /// Maximum translation distance. Set negative to ignore
  Scalar max_translation_distance = -1;
  /// Maximum rotation angle. Set negative to ignore
  Scalar max_angle = -1;
  /// Maximum color RGB distance between corresponding vertices. Set negative to ignore
  Scalar max_color_distance = -1;
  /// The number of points in the sample. We sample this number of points
  /// uniformly from P and Q.
  size_t sample_size = 200;
  /// Maximum time we allow the computation to take. This makes the algorithm
  /// an ANY TIME algorithm that can be stopped at any time, producing the best
  /// solution so far.
  int max_time_seconds = 60;
  /// use a constant default seed by default
  unsigned int randomSeed = std::mt19937::default_seed;

  inline bool configureOverlap(Scalar overlap_, Scalar terminate_threshold_ = Scalar(1)) {
      if(terminate_threshold_ < overlap_) return false;
      overlap_estimation = overlap_;
      terminate_threshold = terminate_threshold_;
      return true;
  }
  inline Scalar getTerminateThreshold() const { return terminate_threshold; }
  inline Scalar getOverlapEstimation()  const { return overlap_estimation; }

private:
  /// Threshold on the value of the target function (LCP, see the paper).
  /// It is used to terminate the process once we reached this value.
  Scalar terminate_threshold = 1.0;
  /// Estimated overlap between P and Q. This is the fraction of points in P that
  /// may have corresponding point in Q. It's being used to estimate the number
  /// of RANSAC iterations needed to guarantee small failure probability.
  Scalar overlap_estimation = 0.2;
};

} /// namespace Super4PCS



#endif //_SHARED_4PCS_H_
