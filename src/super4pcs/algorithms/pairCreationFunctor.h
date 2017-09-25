#ifndef _SUPER4PCS_ALGO_PAIRCREATIONFUNCTOR_H
#define _SUPER4PCS_ALGO_PAIRCREATIONFUNCTOR_H

#include <iostream>
#include <vector>
#include "super4pcs/shared4pcs.h"

#include "super4pcs/accelerators/bbox.h"

#include "super4pcs/accelerators/pairExtraction/bruteForceFunctor.h"
#include "super4pcs/accelerators/pairExtraction/intersectionFunctor.h"
#include "super4pcs/accelerators/pairExtraction/intersectionPrimitive.h"

namespace GlobalRegistration {

template <typename _Scalar>
struct PairCreationFunctor{

public:
  using Scalar      = _Scalar;
  using PairsVector = std::vector<std::pair<int, int>>;
  using VectorType  = typename Point3D::VectorType;

  // Processing data
  Scalar norm_threshold;
  double pair_normals_angle;
  double pair_distance;
  double pair_distance_epsilon;

  // Shared data
  Match4PCSOptions options_;
  const std::vector<Point3D>& Q_;

  PairsVector* pairs;

  std::vector<unsigned int> ids;


  // Internal data
  typedef Eigen::Matrix<Scalar, 3, 1> Point;
  typedef Accelerators::PairExtraction::HyperSphere
  < typename PairCreationFunctor::Point, 3, Scalar> Primitive;

  std::vector< /*Eigen::Map<*/typename PairCreationFunctor::Point/*>*/ > points;
  std::vector< Primitive > primitives;

private:
  VectorType segment1;
  std::vector<Point3D> base_3D_;
  int base_point1_, base_point2_;

  typename PairCreationFunctor::Point _gcenter;
  Scalar _ratio;
  static const typename PairCreationFunctor::Point half;

public:
  inline PairCreationFunctor(
    Match4PCSOptions options,
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

    GlobalRegistration::AABB3D<Scalar> bbox;

    unsigned int nSamples = Q_.size();

    points.reserve(nSamples);
    primitives.reserve(nSamples);

    // Compute bounding box on fine data to be SURE to have all points in the
    // unit bounding box
    for (unsigned int i = 0; i < nSamples; ++i) {
        const VectorType &q = Q_[i].pos();
      points.push_back(q);
      bbox.extend(q);
    }

    _gcenter = bbox.center();
    // add a delta to avoid to have elements with coordinate = 1
    _ratio = bbox.diagonal().maxCoeff() + 0.001;

    // update point cloud (worldToUnit use the ratio and gravity center
    // previously computed)
    // Generate primitives
    for (unsigned int i = 0; i < nSamples; ++i) {
      points[i] = worldToUnit(points[i]);

      primitives.emplace_back(points[i], Scalar(1.));
      ids.push_back(i);
    }
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
                       const std::vector<Point3D>& base_3D){
    base_3D_     = base_3D;
    base_point1_ = base_point1;
    base_point2_ = base_point2;

    segment1 = (base_3D_[base_point2_].pos() -
                base_3D_[base_point1_].pos()).normalized();
  }


  inline void beginPrimitiveCollect(int /*primId*/){ }
  inline void endPrimitiveCollect(int /*primId*/){ }


  //! FIXME Pair filtering is the same than 4pcs. Need refactoring
  inline void process(int i, int j){
    if (i>j){
      const Point3D& p = Q_[j];
      const Point3D& q = Q_[i];

      // Compute the distance and two normal angles to ensure working with
      // wrong orientation. We want to verify that the angle between the
      // normals is close to the angle between normals in the base. This can be
      // checked independent of the full rotation angles which are not yet
      // defined by segment matching alone..
      const Scalar distance = (q.pos() - p.pos()).norm();
#ifndef MULTISCALE
      if (std::abs(distance - pair_distance) > pair_distance_epsilon) return;
#endif

      if ( options_.max_normal_difference > 0 &&
           q.normal().squaredNorm() > 0 &&
           p.normal().squaredNorm() > 0) {
          const Scalar norm_threshold =
                  0.5 * options_.max_normal_difference * M_PI / 180.0;
          const double first_normal_angle = (q.normal() - p.normal()).norm();
          const double second_normal_angle = (q.normal() + p.normal()).norm();
          // Take the smaller normal distance.
          const Scalar first_norm_distance =
                  std::min(std::abs(first_normal_angle - pair_normals_angle),
                           std::abs(second_normal_angle - pair_normals_angle));
          // Verify appropriate angle between normals and distance.

          if (first_norm_distance > norm_threshold) return;
      }
      // Verify restriction on the rotation angle, translation and colors.
      if (options_.max_color_distance > 0) {
          const bool use_rgb = (p.rgb()[0] >= 0 && q.rgb()[0] >= 0 &&
                  base_3D_[base_point1_].rgb()[0] >= 0 &&
                  base_3D_[base_point2_].rgb()[0] >= 0);
          bool color_good = (p.rgb() - base_3D_[base_point1_].rgb()).norm() <
                  options_.max_color_distance &&
                  (q.rgb() - base_3D_[base_point2_].rgb()).norm() <
                  options_.max_color_distance;

          if (use_rgb && ! color_good) return;
      }

      if (options_.max_translation_distance > 0) {
          const bool dist_good = (p.pos() - base_3D_[base_point1_].pos()).norm() <
                  options_.max_translation_distance &&
                  (q.pos() - base_3D_[base_point2_].pos()).norm() <
                  options_.max_translation_distance;
          if (! dist_good) return;
      }

      // need cleaning here
      if (options_.max_angle > 0){
          VectorType segment2 = (q.pos() - p.pos()).normalized();
          if (std::acos(segment1.dot(segment2)) <= options_.max_angle * M_PI / 180.0) {
              pairs->emplace_back(j, i);
          }

          if (std::acos(segment1.dot(- segment2)) <= options_.max_angle * M_PI / 180.0) {
              // Add ordered pair.
              pairs->emplace_back(i, j);
          }
      }else {
          pairs->emplace_back(j, i);
          pairs->emplace_back(i, j);
      }
    }
  }
};

} // namespace Super4PCS

#endif // PAIRCREATIONFUNCTOR_H
