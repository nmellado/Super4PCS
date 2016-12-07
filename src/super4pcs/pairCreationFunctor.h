#ifndef PAIRCREATIONFUNCTOR_H
#define PAIRCREATIONFUNCTOR_H

#include <iostream>
#include <vector>
#include "shared4pcs.h"

template <typename _Scalar>
struct PairCreationFunctor{

private:
  using Point3D = match_4pcs::Point3D;

public:
  using Scalar      = _Scalar;
  using PairsVector = std::vector<std::pair<int, int>>;

  // Processing data
  Scalar norm_threshold;
  double pair_normals_angle;
  double pair_distance;
  double pair_distance_epsilon;

  // Shared data
  match_4pcs::Match4PCSOptions options_;
  const std::vector<Point3D>& Q_;

  PairsVector* pairs;

  std::vector<unsigned int> ids;


  // Internal data
  typedef Eigen::Matrix<Scalar, 3, 1> Point;
  typedef Super4PCS::Accelerators::PairExtraction::HyperSphere
  < typename PairCreationFunctor::Point, 3, Scalar> Primitive;

  std::vector< /*Eigen::Map<*/typename PairCreationFunctor::Point/*>*/ > points;
  std::vector< Primitive > primitives;

private:
  cv::Point3f segment1;
  std::vector<Point3D> base_3D_;
  int base_point1_, base_point2_;

  typename PairCreationFunctor::Point _gcenter;
  Scalar _ratio;
  static const typename PairCreationFunctor::Point half;

public:
  inline PairCreationFunctor(
    match_4pcs::Match4PCSOptions options,
    const std::vector<Point3D>& Q)
    :options_(options), Q_(Q),
     pairs(NULL), _ratio(1.f)
    { }

private:
  inline Point worldToUnit(
    const Eigen::MatrixBase<typename PairCreationFunctor::Point> &p) const {
    static const Point half = Point::Ones() * Scalar(0.5f);
    return (p-_gcenter) / _ratio + half;
  }


public:
  inline Point unitToWorld(
    const Eigen::MatrixBase<typename PairCreationFunctor::Point> &p) const {
    static const Point half = Point::Ones() * Scalar(0.5f);
    return (p - half) * _ratio + _gcenter;
  }


  inline Scalar unitToWorld( Scalar d) const {
    return d * _ratio;
  }


  inline Point getPointInWorldCoord(int i) const {
    return unitToWorld(points[i]);
  }


  inline void synch3DContent(){
    points.clear();
    primitives.clear();

    Super4PCS::AABB3D<Scalar> bbox;

    unsigned int nSamples = Q_.size();

    points.reserve(nSamples);
    primitives.reserve(nSamples);

    // Compute bounding box on fine data to be SURE to have all points in the
    // unit bounding box
    for (unsigned int i = 0; i < nSamples; ++i) {
        PairCreationFunctor::Point q ( Q_[i].x,
                                       Q_[i].y,
                                       Q_[i].z );
      points.push_back(q);
      bbox.extendTo(q);
    }

    _gcenter = bbox.center();
    // add a delta to avoid to have elements with coordinate = 1
    _ratio = std::max(bbox.depth() + 0.001,
             std::max(bbox.width() + 0.001,
                      bbox.height()+ 0.001));

    // update point cloud (worldToUnit use the ratio and gravity center
    // previously computed)
    // Generate primitives
    for (unsigned int i = 0; i < nSamples; ++i) {
      points[i] = worldToUnit(points[i]);

      primitives.push_back(Primitive(points[i], Scalar(1.)));
      ids.push_back(i);
    }

    std::cout << "Work with " << points.size() << " points" << std::endl;
  }

  inline void setRadius(Scalar radius) {
    const Scalar nRadius = radius/_ratio;
    for(typename std::vector< Primitive >::iterator it = primitives.begin();
        it != primitives.end(); ++it)
      (*it).radius() = nRadius;
  }

  inline Scalar getNormalizedEpsilon(Scalar eps){
    return eps/_ratio;
  }

  inline void setBase( int base_point1, int base_point2,
                       std::vector<Point3D>& base_3D){
    base_3D_     = base_3D;
    base_point1_ = base_point1;
    base_point2_ = base_point2;

    segment1 = base_3D_[base_point2_] - base_3D_[base_point1_];
    segment1 *= 1.0 / cv::norm(segment1);
  }


  inline void beginPrimitiveCollect(int /*primId*/){ }
  inline void endPrimitiveCollect(int /*primId*/){ }


  inline void process(int i, int j){
    if (i>j){
      const Point3D& p = Q_[j];
      const Point3D& q = Q_[i];

#ifndef MULTISCALE
      const float distance = cv::norm(q - p);
      if (std::abs(distance - pair_distance) > pair_distance_epsilon) return;
#endif
      const bool use_normals = norm(q.normal()) > 0 && norm(p.normal()) > 0;
      bool normals_good = true;
      if (use_normals) {
        const Scalar first_normal_angle = cv::norm(q.normal() - p.normal());
        const Scalar second_normal_angle = cv::norm(q.normal() + p.normal());
        // Take the smaller normal distance.
        const Scalar first_norm_distance =
            std::min(std::abs(first_normal_angle - pair_normals_angle),
                     std::abs(second_normal_angle - pair_normals_angle));
        // Verify appropriate angle between normals and distance.
        normals_good = first_norm_distance < norm_threshold;
      }
      if (!normals_good) return;
      cv::Point3f segment2 = q - p;
      segment2 *= 1.0 / cv::norm(segment2);
      // Verify restriction on the rotation angle, translation and colors.
      const bool use_rgb = (p.rgb()[0] >= 0 && q.rgb()[0] >= 0 &&
                            base_3D_[base_point1_].rgb()[0] >= 0 &&
                            base_3D_[base_point2_].rgb()[0] >= 0);
      const bool rgb_good =
          use_rgb ? cv::norm(p.rgb() - base_3D_[base_point1_].rgb()) <
                            options_.max_color_distance &&
                        cv::norm(q.rgb() - base_3D_[base_point2_].rgb()) <
                            options_.max_color_distance
                  : true;
      const bool dist_good = cv::norm(p - base_3D_[base_point1_]) <
                                 options_.max_translation_distance &&
                             cv::norm(q - base_3D_[base_point2_]) <
                                 options_.max_translation_distance;

      if (std::acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0 &&
          dist_good && rgb_good) {
        // Add ordered pair.
        pairs->push_back(std::pair<int, int>(j, i));
      }
      // The same for the second order.
      segment2 = p - q;
      segment2 *= 1.0 / cv::norm(segment2);
      if (std::acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0 &&
          dist_good && rgb_good) {
        // Add ordered pair.
        pairs->push_back(std::pair<int, int>(i, j));
      }
    }
  }
};

#endif // PAIRCREATIONFUNCTOR_H
