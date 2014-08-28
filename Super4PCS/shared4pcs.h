#ifndef _SHARED_4PCS_H_
#define _SHARED_4PCS_H_

#include <vector>
#include <memory>

#include <opencv2/core/core.hpp>


namespace match_4pcs {
// The basic 3D point structure. A point potentially contains also directional
// information and color.
class Point3D : public cv::Point3f {
 public:
  inline Point3D(double x, double y, double z) : cv::Point3f(x, y, z) {}
  inline Point3D(const cv::Point3f& other): cv::Point3f(other) {}
  inline Point3D() : cv::Point3f(0.0f, 0.0f, 0.0f) {}
  inline const cv::Vec3f& rgb() const { return rgb_; }
  inline const cv::Point3f& normal() const { return normal_; }
  inline void set_rgb(const cv::Vec3f& rgb) { rgb_ = rgb; }
  inline void set_normal(const cv::Point3d& normal) { normal_ = normal; }
  inline void set_normal(const cv::Point3f& normal) { normal_ = normal; }
  inline void normalize() { 
    double n = cv::norm(*this);
    x /= n;
    y /= n;
    z /= n;
  }

 private:
  // Normal.
  cv::Point3f normal_{0.0f, 0.0f, 0.0f};
  // Color.
  cv::Vec3f rgb_{-1.0f, -1.0f, -1.0f};
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
  Match4PCSOptions() {};

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
} // namespace match_4pcs

#endif //_SHARED_4PCS_H_
