// Copyright 2014 Nicolas Mellado
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
// Authors: Nicolas Mellado, Dror Aiger
//
// An implementation of the Super 4-points Congruent Sets (Super 4PCS) 
// algorithm presented in:
//
// Super 4PCS: Fast Global Pointcloud Registration via Smart Indexing
// Nicolas Mellado, Dror Aiger, Niloy J. Mitra
// Symposium on Geometry Processing 2014.
//
// Data acquisition in large-scale scenes regularly involves accumulating 
// information across multiple scans. A common approach is to locally align scan 
// pairs using Iterative Closest Point (ICP) algorithm (or its variants), but 
// requires static scenes and small motion between scan pairs. This prevents 
// accumulating data across multiple scan sessions and/or different acquisition
// modalities (e.g., stereo, depth scans). Alternatively, one can use a global
// registration algorithm allowing scans to be in arbitrary initial poses. The 
// state-of-the-art global registration algorithm, 4PCS, however has a quadratic
// time complexity in the number of data points. This vastly limits its 
// applicability to acquisition of large environments. We present Super 4PCS for
// global pointcloud registration that is optimal, i.e., runs in linear time (in 
// the number of data points) and is also output sensitive in the complexity of 
// the alignment problem based on the (unknown) overlap across scan pairs. 
// Technically, we map the algorithm as an ‘instance problem’ and solve it 
// efficiently using a smart indexing data organization. The algorithm is 
// simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in
// significant speedup over alternative approaches and allows unstructured 
// efficient acquisition of scenes at scales previously not possible. Complete 
// source code and datasets are available for research use at 
// http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.

#include "4pcs.h"

#include "Eigen/Core"
#include "Eigen/Geometry"                 // MatrixBase.homogeneous()
#include "Eigen/SVD"                      // Transform.computeRotationScaling()

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>

#include "accelerators/intersection.h"
#include "accelerators/primitives.h"
#include "accelerators/normalset.h"
#include "accelerators/normalHealSet.h"
#include "accelerators/kdtree.h"
#include "accelerators/bbox.h"

#include <hash_map>
#include <fstream>

#define sqr(x) ((x)*(x))
#define norm2(p) (sqr(p.x)+sqr(p.y)+sqr(p.z))

#ifndef MAX
#define MAX(a,b) a<b?a:b
#endif

//#define TEST_GLOBAL_TIMINGS

typedef std::vector<std::pair<int, int>> PairsVector;
typedef double Scalar;

#ifdef TEST_GLOBAL_TIMINGS
#include <chrono> //timers


class Timer {
public:
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::nanoseconds timestep;

    explicit inline Timer(bool run = false)
    {
        if (run) reset();
    }
    void reset()
    {
        _start = clock::now();
    }
    inline timestep elapsed() const
    {
        return std::chrono::duration_cast<timestep>(clock::now() - _start);
    }
    template <typename T, typename Traits>
    friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const Timer& timer)
    {
        return out << timer.elapsed().count();
    }
private:
    clock::time_point _start;
};

Scalar totalTime;
Scalar kdTreeTime;
Scalar verifyTime;

#endif

namespace match_4pcs {

using namespace std;

namespace {

class HashTable {
 private:
  const uint64 MAGIC1 = 100000007;
  const uint64 MAGIC2 = 161803409;
  const uint64 MAGIC3 = 423606823;
  const uint64 NO_DATA = 0xffffffffu;
  float voxel_;
  float scale_;
  vector<vector<int>> voxels_;
  vector<uint64> data_;

 public:
  HashTable(int maxpoints, float voxel) : voxel_(voxel), scale_(1.0f / voxel) {
    uint64 n = maxpoints;
    voxels_.resize(n);
    data_.resize(n, NO_DATA);
  }
  uint64& operator[](const Point3D& p) {
    vector<int> c(3);
    c[0] = static_cast<int>(floor(p.x * scale_));
    c[1] = static_cast<int>(floor(p.y * scale_));
    c[2] = static_cast<int>(floor(p.z * scale_));
    uint64 key = (MAGIC1 * c[0] + MAGIC2 * c[1] + MAGIC3 * c[2]) % data_.size();
    while (1) {
      if (data_[key] == NO_DATA) {
        voxels_[key] = c;
        break;
      } else if (voxels_[key] == c) {
        break;
      }
      key++;
      if (key == data_.size()) key = 0;
    }
    return data_[key];
  }
};

const int kNumberOfDiameterTrials = 1000;

// Local static helpers that do not depend on state.

inline float Square(float x) { return x * x; }

void DistUniformSampling(const std::vector<Point3D>& set, float delta,
                         std::vector<Point3D>* sample) {
  int num_input = set.size();
  sample->clear();
  HashTable hash(num_input, delta);
  for (int i = 0; i < num_input; i++) {
    uint64& ind = hash[set[i]];
    if (ind >= num_input) {
      sample->push_back(set[i]);
      ind = sample->size();
    }
  }
}

inline float PointsDistance(const Point3D& p, const Point3D& q) {
  return cv::norm(p - q);
}

// Compute the closest points between two 3D line segments and obtain the two
// invariants corresponding to the closet points. This is the "intersection"
// point that determines the invariants. Since the 4 points are not exactly
// planar, we use the center of the line segment connecting the two closest
// points as the "intersection".
float distSegmentToSegment(const cv::Point3f& p1, const cv::Point3f& p2,
                           const cv::Point3f& q1, const cv::Point3f& q2,
                           double* invariant1, double* invariant2) {
  if (invariant1 == 0 || invariant2 == 0) return kLargeNumber;
  const float kSmallNumber = 0.0001;
  cv::Point3f u = p2 - p1;
  cv::Point3f v = q2 - q1;
  cv::Point3f w = p1 - q1;
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
  cv::Point3f distance = w + (*invariant1 * u) - (*invariant2 * v);
  return cv::norm(distance);
}

}  // namespace

// Holds a base from P. The base contains 4 points (indices) from the set P.
struct Quadrilateral {
  int vertices[4];
  Quadrilateral(int vertex0, int vertex1, int vertex2, int vertex3) {
    vertices[0] = vertex0;
    vertices[1] = vertex1;
    vertices[2] = vertex2;
    vertices[3] = vertex3;
//    cout << vertex0 << " "
//         << vertex1 << " "
//         << vertex2 << " "
//         << vertex3 << " "
//         << (float(vertex0) + float(vertex1) + float(vertex2) + float(vertex3))/4.f
//         << endl;
  }
};


struct PairCreationFunctor{

public:
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
  typedef Super4PCS::HyperSphere< PairCreationFunctor::Point, 3, Scalar> Primitive;

  std::vector< /*Eigen::Map<*/PairCreationFunctor::Point/*>*/ > points;
  std::vector< Primitive > primitives;

private:
  cv::Point3f segment1;
  std::vector<Point3D> base_3D_;
  int base_point1_, base_point2_;

  PairCreationFunctor::Point _gcenter;
  Scalar _ratio;

public:
  inline PairCreationFunctor(
    Match4PCSOptions options,
    const std::vector<Point3D>& Q)
    :options_(options), Q_(Q),
     pairs(NULL), _ratio(1.f)
    { }


  inline PairCreationFunctor::Point worldToUnit(
    const PairCreationFunctor::Point &p) const {

    static const PairCreationFunctor::Point half =
              PairCreationFunctor::Point::Ones() * Scalar(0.5f);
    return (p-_gcenter) / _ratio + half;
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
    for (int i = 0; i < nSamples; ++i) {
        PairCreationFunctor::Point q ( Q_[i].x,
                                       Q_[i].y,
                                       Q_[i].z );
      points.push_back(q);
      bbox.extendTo(q);
    }

    _gcenter = bbox.center();
    // add a delta to avoid to have elements with coordinate = 1
    _ratio = MAX(bbox.depth() + 0.001,
             MAX(bbox.width() + 0.001,
                 bbox.height()+ 0.001));

    // update point cloud (worldToUnit use the ratio and gravity center
    // previously computed)
    // Generate primitives
    for (int i = 0; i < nSamples; ++i) {
      points[i] = worldToUnit(points[i]);

      primitives.push_back(Primitive(points[i], Scalar(1.)));
      ids.push_back(i);
    }

    cout << "Work with " << points.size() << " points" << endl;
  }

  inline void setRadius(Scalar radius) {
    const Scalar nRadius = radius/_ratio;
    for(std::vector< Primitive >::iterator it = primitives.begin();
        it != primitives.end(); it++)
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

      const float distance = cv::norm(q - p);
      if (fabs(distance - pair_distance) > pair_distance_epsilon) return;
      const bool use_normals = norm(q.normal()) > 0 && norm(p.normal()) > 0;
      bool normals_good = true;
      if (use_normals) {
        const double first_normal_angle = cv::norm(q.normal() - p.normal());
        const double second_normal_angle = cv::norm(q.normal() + p.normal());
        // Take the smaller normal distance.
        const float first_norm_distance =
            min(fabs(first_normal_angle - pair_normals_angle),
                fabs(second_normal_angle - pair_normals_angle));
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
};

class MatchSuper4PCSImpl {
 public:
  explicit MatchSuper4PCSImpl(const Match4PCSOptions& options)
      : number_of_trials_(0),
        max_base_diameter_(-1),
        P_mean_distance_(1.0),
        best_LCP_(0.0F),
        options_(options),
        pcfunctor_(options_, sampled_Q_3D_) {
    base_3D_.resize(4);
  }

  ~MatchSuper4PCSImpl() {
    Clear();
  }

  void Clear() {
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

  // Number of trials. Every trial picks random base from P.
  int number_of_trials_;
  // Maximum base diameter. It is computed automatically from the diameter of
  // P and the estimated overlap and used to limit the distance between the
  // points in the base in P so that the probability to have all points in
  // the base as inliers is increased.
  float max_base_diameter_;
  // The diameter of P.
  float P_diameter_;
  // KdTree used to compute the LCP
  Super4PCS::KdTree<Scalar> kd_tree_;
  // Mean distance between points and their nearest neighbor in the set P.
  // Used to normalize the "delta" which is given in terms of this distance.
  float P_mean_distance_;
  // The transformation matrix by wich we transform Q to P
  Eigen::Matrix<Scalar, 4, 4> transform_;
  // Quad centroids in first and second clouds
  // They are used temporarily and makes the transformations more robust to
  // noise. At the end, the direct transformation applied as a 4x4 matrix on
  // every points in Q is computed and returned.
  Eigen::Matrix<Scalar, 3, 1> qcentroid1_, qcentroid2_;
  // The points in the base (indices to P). It is being updated in every
  // RANSAC iteration.
  int base_[4];
  // The current congruent 4 points from Q. Every RANSAC iteration the
  // algorithm examines a set of such congruent 4-points from Q and retains
  // the best from them (the one that realizes the best LCP).
  int current_congruent_[4];
  // Sampled P (3D coordinates).
  std::vector<Point3D> sampled_P_3D_;
  // Sampled Q (3D coordinates).
  std::vector<Point3D> sampled_Q_3D_;
  // The 3D points of the base.
  std::vector<Point3D> base_3D_;
  // The copy of the input Q. We transform Q to match P and returned the
  // transformed version.
  std::vector<Point3D> Q_copy_;
  // The centroid of P.
  cv::Point3f centroid_P_;
  // The centroid of Q.
  cv::Point3f centroid_Q_;
  // The best LCP (Largest Common Point) fraction so far.
  float best_LCP_;
  // Current trial.
  int current_trial_;
  // Parameters.
  Match4PCSOptions options_;

  PairCreationFunctor pcfunctor_;

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
  void
  BruteForcePairs(double pair_distance, double pair_normals_angle,
                       double pair_distance_epsilon, int base_point1,
                       int base_point2,
                       PairsVector* pairs);

  // For each randomly picked base, verifies the computed transformation by
  // computing the number of points that this transformation brings near points
  // in Q. Returns the current LCP. R is the rotation matrix, (tx,ty,tz) is
  // the translation vector and (cx,cy,cz) is the center of transformation.template <class MatrixDerived>
  Scalar Verify(const Eigen::Matrix<Scalar, 4, 4>& mat);

  // Computes the mean distance between point in Q and its nearest neighbor.
  double MeanDistance();

  // Selects a quadrilateral from P and returns the corresponding invariants
  // and point indices. Returns true if a quadrilateral has been found, false
  // otherwise.
  bool SelectQuadrilateral(double* invariant1, double* invariant2, int* base1,
                           int* base2, int* base3, int* base4);

  // Select random triangle in P such that its diameter is close to
  // max_base_diameter_. This enables to increase the probability of having
  // all three points in the inlier set. Return true on success, false if such a
  // triangle cannot be found.
  bool SelectRandomTriangle(int* base1, int* base2, int* base3);

  // Takes quadrilateral as a base, computes robust intersection point
  // (approximate as the lines might not intersect) and returns the invariants
  // corresponding to the two selected lines. The method also updates the order
  // of the base base_3D_.
  bool TryQuadrilateral(double* invariant1, double* invariant2, int &base1, int &base2, int &base3, int &base4);

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
  // @param [in] Q_pointsFirst Point coordinates for the pairs first id
  // @param [in] Q_pointsSecond Point coordinates for the pairs second id
  // @param [out] quadrilaterals The set of congruent quadrilateral. In fact,
  // it's a super set from which we extract the real congruent set.
  bool FindCongruentQuadrilaterals(double invariant1, double invariant2,
                                   double distance_threshold1,
                                   double distance_threshold2,
                                   const PairsVector& P_pairs,
                                   const PairsVector& Q_pairs,
                                   const std::vector<Point3D>& Q_points,
                                   std::vector<Quadrilateral>* quadrilaterals);

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
  // @param [in] Q_pointsFirst Point coordinates for the pairs first id
  // @param [in] Q_pointsSecond Point coordinates for the pairs second id
  // @param [out] quadrilaterals The set of congruent quadrilateral. In fact,
  // it's a super set from which we extract the real congruent set.
  bool FindCongruentQuadrilateralsFast(double invariant1, double invariant2,
                                       double distance_threshold1,
                                       double distance_threshold2,
                                       const PairsVector& P_pairs,
                                       const PairsVector& Q_pairs,
                                       const std::vector<Point3D>& Q_points,
                                       std::vector<Quadrilateral>* quadrilaterals);

  // Computes the best rigid transformation between three corresponding pairs.
  // The transformation is characterized by rotation matrix, translation vector
  // and a center about which we rotate. The set of pairs is potentially being
  // updated by the best permutation of the second set. Returns the RMS of the
  // fit. The method is being called with 4 points but it applies the fit for
  // only 3 after the best permutation is selected in the second set (see
  // bellow). This is done because the solution for planar points is much
  // simpler.
  // The method is the closed-form solution by Horn:
  // people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf‎
  bool ComputeRigidTransformation(const vector< pair<Point3D, Point3D> >& pairs,
                                  const Eigen::Matrix<Scalar, 3, 1>& centroid1,
                                  const Eigen::Matrix<Scalar, 3, 1>& centroid2,
                                  Scalar max_angle,
                                  Eigen::Matrix<Scalar, 4, 4> &transform,
                                  Scalar& rms_,
                                  bool computeScale );

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
bool MatchSuper4PCSImpl::FindCongruentQuadrilateralsFast(
    double invariant1, double invariant2, double distance_threshold1,
    double distance_threshold2, const std::vector<std::pair<int, int>>& P_pairs,
    const std::vector<std::pair<int, int>>& Q_pairs,
    const std::vector<Point3D>& Q_points,
    std::vector<Quadrilateral>* quadrilaterals) {

  // Compute the angle formed by the two vectors of the basis
  Point3D b1 = base_3D_[1] - base_3D_[0];  b1.normalize();
  Point3D b2 = base_3D_[3] - base_3D_[2];  b2.normalize();
  double alpha = /*std::abs(*/b1.dot(b2)/*)*/;

  // 1. Datastructure construction
  typedef PairCreationFunctor::Point Point;
  const Scalar eps = pcfunctor_.getNormalizedEpsilon(
    std::min(distance_threshold1, distance_threshold2));
  typedef Super4PCS::IndexedNormalHealSet IndexedNormalSet3D;
  
  // Use the following definition to get ride of Healpix
//  typedef  IndexedNormalSet
//                  < Point,   //! \brief Point type used internally
//                    3,       //! \brief Nb dimension
//                    7,       //! \brief Nb cells/dim normal
//                    Scalar>  //! \brief Scalar type
//  IndexedNormalSet3D;


  IndexedNormalSet3D nset (eps, 2.);

  for (int i = 0; i < P_pairs.size(); ++i) {
    const Point& p1 = pcfunctor_.points[P_pairs[i].first];
    const Point& p2 = pcfunctor_.points[P_pairs[i].second];
    Point  n  = (p2 - p1).normalized();

//    cout << "new entry: " << endl
//         << p1.transpose() << "(" << P_pairs[i].first  << ")" << endl
//         << p2.transpose() << "(" << P_pairs[i].second << ")" << endl
//         << (p1+ invariant1       * (p2 - p1)).transpose() << endl
//         << n.transpose() << endl;

    nset.addElement( p1+ invariant1 * (p2 - p1),  n, i);
  }


  std::set< std::pair<unsigned int, unsigned int > > comb;

  unsigned int j = 0;
  std::vector<unsigned int> nei;
  // 2. Query time
  for (int i = 0; i < Q_pairs.size(); ++i) {
    const Point& p1 = pcfunctor_.points[Q_pairs[i].first];
    const Point& p2 = pcfunctor_.points[Q_pairs[i].second];

    const Point3D& pq1 = Q_points[Q_pairs[i].first];
    const Point3D& pq2 = Q_points[Q_pairs[i].second];

    nei.clear();

    Point   query  =  p1 + invariant2 * ( p2 - p1 );
    Point3D queryQ = pq1 + invariant2 * (pq2 - pq1);

    Point queryn = (p2 - p1).normalized();

//    cout << "query: " << endl
//         << p1.transpose() << "(" << Q_pairs[i].first  << ")" << endl
//         << p2.transpose() << "(" << Q_pairs[i].second << ")" << endl
//         << query.transpose() << endl
//         << queryn.transpose() << endl;

    nset.getNeighbors( query,  queryn, alpha, nei, false);


    Point3D invPoint;
    //const float distance_threshold2s = distance_threshold2 * distance_threshold2;
    for (unsigned int k = 0; k != nei.size(); k++){
      int id = nei[k];

      const Point3D& pp1 = Q_points[P_pairs[id].first];
      const Point3D& pp2 = Q_points[P_pairs[id].second];

      invPoint = pp1 + (pp2 - pp1) * invariant1;

       // use also distance_threshold2 for inv 1 and 2 in 4PCS
      if (cv::norm(queryQ-invPoint) <= distance_threshold2){
          comb.insert(std::pair<unsigned int, unsigned int>(id, i));
      }
    }
  }

  for (std::set< std::pair<unsigned int, unsigned int > >::const_iterator it =
             comb.cbegin();
       it != comb.cend(); it++){
    const unsigned int & id = (*it).first;
    const unsigned int & i  = (*it).second;

    quadrilaterals->push_back(
                Quadrilateral(P_pairs[id].first, P_pairs[id].second,
                              Q_pairs[i].first,  Q_pairs[i].second));
  }

  return quadrilaterals->size() != 0;
}

// Finds congruent candidates in the set Q, given the invariants and threshold
// distances.
bool MatchSuper4PCSImpl::FindCongruentQuadrilaterals(
    double /*invariant1*/,
    double /*invariant2*/,
    double /*distance_threshold1*/,
    double /*distance_threshold2*/,
    const std::vector<std::pair<int, int>>& /*P_pairs*/,
    const std::vector<std::pair<int, int>>& /*Q_pairs*/,
    const std::vector<Point3D>& /*Q_points*/,
    std::vector<Quadrilateral>* /*quadrilaterals*/) {

    cerr << "What the hell are you doing here ?? " << __FILE__ << ":" << __LINE__ << endl;
    exit(-1);

//  if (quadrilaterals == NULL) return false;

//  int number_of_points = 2 * P_pairs.size();

//  // We need a temporary ANN tree to store the new points corresponding to
//  // invariants in the P_pairs and then query them (for range search) for all
//  // the new points corresponding to the invariants in Q_pairs.
//  ANNpointArray data_points = annAllocPts(number_of_points, 3);
//  ANNpoint query_point = annAllocPt(3);
//  ANNidxArray near_neighbor_index = new ANNidx[number_of_points];
//  ANNdistArray distances = new ANNdist[number_of_points];

//  quadrilaterals->clear();

//  // Build the ANN tree using the invariants on P_pairs.
//  for (int i = 0; i < P_pairs.size(); ++i) {
//    const Point3D& p1 = Q_points[P_pairs[i].first];
//    const Point3D& p2 = Q_points[P_pairs[i].second];

//    // Init KdTree buffer
//    data_points[i * 2][0] = p1.x + invariant1 * (p2.x - p1.x);
//    data_points[i * 2][1] = p1.y + invariant1 * (p2.y - p1.y);
//    data_points[i * 2][2] = p1.z + invariant1 * (p2.z - p1.z);
//    data_points[i * 2 + 1][0] = p1.x + invariant2 * (p2.x - p1.x);
//    data_points[i * 2 + 1][1] = p1.y + invariant2 * (p2.y - p1.y);
//    data_points[i * 2 + 1][2] = p1.z + invariant2 * (p2.z - p1.z);
//  }

//  ANNkd_tree* tree = new ANNkd_tree(data_points, number_of_points, 3);

//  // Compute the angle formed by the two vectors of the basis
//  Point3D b1 = base_3D_[0] - base_3D_[1];  b1.normalize();
//  Point3D b2 = base_3D_[2] - base_3D_[3];  b2.normalize();
//  double alpha = b1.dot(b2);

//  // Query the ANN for all the points corresponding to the invariants in Q_pair.
//  for (int i = 0; i < Q_pairs.size(); ++i) {
//    const Point3D& p1 = Q_points[Q_pairs[i].first];
//    const Point3D& p2 = Q_points[Q_pairs[i].second];

//    query_point[0] = p1.x + invariant1 * (p2.x - p1.x);
//    query_point[1] = p1.y + invariant1 * (p2.y - p1.y);
//    query_point[2] = p1.z + invariant1 * (p2.z - p1.z);
//    tree->annkFRSearch(query_point, distance_threshold2, number_of_points,
//                       near_neighbor_index, distances, 0);

//    // This is a new candidate of a quadrilateral.
//    for (int j = 0; j < number_of_points; ++j) {
//      if (distances[j] != ANN_DIST_INF) {
//        int id = near_neighbor_index[j] / 2;

//          quadrilaterals->push_back(
//              Quadrilateral(P_pairs[id].first, P_pairs[id].second,
//                            Q_pairs[i].first, Q_pairs[i].second));
//      } else
//        break;
//    }



//    // We test the other order as our pairs are not ordered.
//    query_point[0] = p1.x + invariant2 * (p2.x - p1.x);
//    query_point[1] = p1.y + invariant2 * (p2.y - p1.y);
//    query_point[2] = p1.z + invariant2 * (p2.z - p1.z);
//    tree->annkFRSearch(query_point, distance_threshold2, number_of_points,
//                       near_neighbor_index, distances, 0);

//    for (int j = 0; j < number_of_points; ++j) {
//      if (distances[j] != ANN_DIST_INF) {
//        int id = near_neighbor_index[j] / 2;
//          quadrilaterals->push_back(
//              Quadrilateral(P_pairs[id].first, P_pairs[id].second,
//                            Q_pairs[i].first, Q_pairs[i].second));
//      } else
//        break;
//    }
//  }

//  annDeallocPt(query_point);
//  annDeallocPts(data_points);
//  delete[] near_neighbor_index;
//  delete[] distances;
//  delete tree;

//  return quadrilaterals->size() != 0;
}

// Try the current base in P and obtain the best pairing, i.e. the one that
// gives the smaller distance between the two closest points. The invariants
// corresponding the the base pairing are computed.
bool MatchSuper4PCSImpl::TryQuadrilateral(double* invariant1, double* invariant2,
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
  vector<Point3D> tmp = base_3D_;
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

// Selects a random triangle in the set P (then we add another point to keep the
// base as planar as possible). We apply a simple heuristic that works in most
// practical cases. The idea is to accept maximum distance, computed by the
// estimated overlap, multiplied by the diameter of P, and try to have
// a triangle with all three edges close to this distance. Wide triangles helps
// to make the transformation robust while too large triangles makes the
// probability of having all points in the inliers small so we try to trade-off.
bool MatchSuper4PCSImpl::SelectRandomTriangle(int* base1, int* base2, int* base3) {
  if (base1 == NULL || base2 == NULL || base3 == NULL) return false;
  
  int number_of_points = sampled_P_3D_.size();
  *base1 = *base2 = *base3 = -1;

  // Pick the first point at random.
  int first_point = rand() % number_of_points;

  // Try fixed number of times retaining the best other two.
  float best_wide = 0.0;
  for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
    // Pick and compute
    int second_point = rand() % number_of_points;
    int third_point = rand() % number_of_points;
    cv::Point3f u = sampled_P_3D_[second_point] - sampled_P_3D_[first_point];
    cv::Point3f w = sampled_P_3D_[third_point] - sampled_P_3D_[first_point];
    // We try to have wide triangles but still not too large.
    float how_wide = cv::norm(u.cross(w));
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

// Selects a good base from P and computes its invariants. Returns false if
// a good planar base cannot can be found.
bool MatchSuper4PCSImpl::SelectQuadrilateral(double* invariant1, double* invariant2,
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
      for (int i = 0; i < sampled_P_3D_.size(); ++i) {
        double d1 = PointsDistance(sampled_P_3D_[i], sampled_P_3D_[*base1]);
        double d2 = PointsDistance(sampled_P_3D_[i], sampled_P_3D_[*base2]);
        double d3 = PointsDistance(sampled_P_3D_[i], sampled_P_3D_[*base3]);
        float too_small = max_base_diameter_ * kBaseTooSmall;
        if (d1 >= too_small && d2 >= too_small && d3 >= too_small) {
          // Not too close to any of the first 3.
          double distance =
              fabs(A * sampled_P_3D_[i].x + B * sampled_P_3D_[i].y +
                   C * sampled_P_3D_[i].z - 1.0);
          // Search for the most planar.
          if (distance < best_distance) {
            best_distance = distance;
            *base4 = i;
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

// Computes the mean distance between points in Q and their nearest neighbor.
// We need this for normalization of the user delta (See the paper) to the
// "scale" of the set.
double MatchSuper4PCSImpl::MeanDistance() {
  const float kDiameterFraction = 0.2;
  Super4PCS::KdTree<Scalar>::VectorType query_point;

  int number_of_samples = 0;
  Scalar distance = 0.0;

  for (int i = 0; i < sampled_P_3D_.size(); ++i) {
    query_point << sampled_P_3D_[i].x, sampled_P_3D_[i].y, sampled_P_3D_[i].z;

    typename Super4PCS::KdTree<Scalar>::Index resId =
    kd_tree_.doQueryRestrictedClosest(query_point, P_diameter_ * kDiameterFraction, i);

    if (resId != Super4PCS::KdTree<Scalar>::invalidIndex()) {
      distance += cv::norm(sampled_P_3D_[i] - sampled_P_3D_[resId]);
      number_of_samples++;
    }
  }

  return distance / number_of_samples;
}

// Verify a given transformation by computing the number of points in P at
// distance at most (normalized) delta from some point in Q. In the paper
// we describe randomized verification. We apply deterministic one here with
// early termination. It was found to be fast in practice.
Scalar MatchSuper4PCSImpl::Verify(const Eigen::Matrix<Scalar, 4, 4>& mat) {

#ifdef TEST_GLOBAL_TIMINGS
    Timer t_verify (true);
#endif

  // We allow factor 2 scaling in the normalization.
  float epsilon = options_.delta;
  int good_points = 0;
  int number_of_points = sampled_Q_3D_.size();
  int terminate_value = best_LCP_ * number_of_points;

  Eigen::Matrix<Scalar, 4, 1> query_point;

  Scalar sq_eps = epsilon*epsilon;

  for (int i = 0; i < number_of_points; ++i) {
    Point3D p = sampled_Q_3D_[i];
    //Transform(rotation, center, translate, &p);

    query_point << p.x,
                   p.y,
                   p.z,
                   Scalar(1);
    query_point = (mat * query_point).eval();


    // Use the kdtree to get the nearest neighbor
#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif
    typename Super4PCS::KdTree<Scalar>::Index resId =
    kd_tree_.doQueryRestrictedClosest(query_point.head<3>(), sq_eps);

#ifdef TEST_GLOBAL_TIMINGS
    kdTreeTime += Scalar(t.elapsed().count()) / Scalar(1000000.);
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
  verifyTime += Scalar(t_verify.elapsed().count()) / Scalar(1000000.);
#endif
  return static_cast<float>(good_points) / number_of_points;
}

// Constructs two sets of pairs in Q, each corresponds to one pair in the base
// in P, by having the same distance (up to some tolerantz) and optionally the
// same angle between normals and same color.
void
MatchSuper4PCSImpl::BruteForcePairs(double pair_distance,
                                    double pair_normals_angle,
                                    double pair_distance_epsilon,
                                    int base_point1, int base_point2,
                                    PairsVector* pairs) {

  pcfunctor_.pairs = pairs;

  pairs->clear();
  pairs->reserve(2 * pcfunctor_.Q_.size());

  pcfunctor_.pair_distance         = pair_distance;
  pcfunctor_.pair_distance_epsilon = pair_distance_epsilon;
  pcfunctor_.pair_normals_angle    = pair_normals_angle;
  pcfunctor_.norm_threshold =
      0.5 * options_.max_normal_difference * M_PI / 180.0;

  pcfunctor_.setRadius(pair_distance);
  pcfunctor_.setBase(base_point1, base_point2, base_3D_);
  
  Super4PCS::IntersectionFunctor
  <PairCreationFunctor::Primitive,
   PairCreationFunctor::Point, 3, Scalar> interFunctor;

  Scalar eps = pcfunctor_.getNormalizedEpsilon(pair_distance_epsilon);

  interFunctor.process(pcfunctor_.primitives,
                       pcfunctor_.points,
                       eps,
                       50,
                       pcfunctor_);
}

// Pick one base, finds congruent 4-points in Q, verifies for all
// transformations, and retains the best transformation and LCP. This is
// a complete RANSAC iteration.
bool MatchSuper4PCSImpl::TryOneBase() {
  vector<pair<Point3D, Point3D>> congruent_points(4);
  double invariant1, invariant2;
  int base_id1, base_id2, base_id3, base_id4;
  float distance_factor = 2.0;

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

  if (!SelectQuadrilateral(&invariant1, &invariant2, &base_id1, &base_id2,
                           &base_id3, &base_id4)) {
    return false;
  }
#endif

  // Computes distance between pairs.
  double distance1 = PointsDistance(base_3D_[0], base_3D_[1]);
  double distance2 = PointsDistance(base_3D_[2], base_3D_[3]);

  vector<pair<int, int>> pairs1, pairs2;
  vector<Quadrilateral> congruent_quads;

  // Compute normal angles.
  double normal_angle1 = cv::norm(base_3D_[0].normal() - base_3D_[1].normal());
  double normal_angle2 = cv::norm(base_3D_[2].normal() - base_3D_[3].normal());

  BruteForcePairs(distance1, normal_angle1, distance_factor * options_.delta, 0,
                  1, &pairs1);
  BruteForcePairs(distance2, normal_angle2, distance_factor * options_.delta, 2,
                  3, &pairs2);

  if (pairs1.size() == 0 || pairs2.size() == 0) {
    return false;
  }

  float factor;
  if (sampled_Q_3D_.size() <= 200) factor = 3.5;
  else if (sampled_Q_3D_.size() >=500) factor = 1.5;
  else factor = 3.5 - (sampled_Q_3D_.size() - 200) * (3.5 - 1.5) / (500.0 - 200.0);
  if (!FindCongruentQuadrilateralsFast(invariant1, invariant2,
                                   distance_factor * factor * options_.delta,
                                   distance_factor * factor * options_.delta,
                                   pairs1,
                                   pairs2,
                                   sampled_Q_3D_,
                                   &congruent_quads)) {
    return false;
  }

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
  congruent_points.resize(4);
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
                               false);             // state: compute scale ratio ?

    if (ok) {

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

struct eqstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1, s2) == 0;
  }
};

// Initialize all internal data structures and data members.
void MatchSuper4PCSImpl::Initialize(const std::vector<Point3D>& P,
                               const std::vector<Point3D>& Q) {

  const float kSmallError = 0.00001;
  const int kMinNumberOfTrials = 4;
  const float kDiameterFraction = 0.3;

  centroid_P_.x = 0;
  centroid_P_.y = 0;
  centroid_P_.z = 0;
  centroid_Q_.x = 0;
  centroid_Q_.y = 0;
  centroid_Q_.z = 0;

  sampled_P_3D_.clear();
  sampled_Q_3D_.clear();

#ifdef TEST_GLOBAL_TIMINGS
    kdTreeTime = 0;
    totalTime  = 0;
    verifyTime = 0;
#endif


  // prepare P
  if (P.size() > options_.sample_size){
      std::vector<Point3D> uniform_P;
      DistUniformSampling(P, options_.delta, &uniform_P);

      int sample_fraction_P = 1;  // We prefer not to sample P but any number can be
                                  // placed here.

      // Sample the sets P and Q uniformly.
      for (int i = 0; i < uniform_P.size(); ++i) {
        if (rand() % sample_fraction_P == 0) {
          sampled_P_3D_.push_back(uniform_P[i]);
        }
      }
  }
  else
  {
      cout << "(P) More samples requested than available: use whole cloud" << endl;
      sampled_P_3D_ = P;
  }



  // prepare Q
  if (Q.size() > options_.sample_size){
      std::vector<Point3D> uniform_Q;
      DistUniformSampling(Q, options_.delta, &uniform_Q);
      int sample_fraction_Q =
          max(1, static_cast<int>(uniform_Q.size() / options_.sample_size));

      for (int i = 0; i < uniform_Q.size(); ++i) {
        if (rand() % sample_fraction_Q == 0) {
          sampled_Q_3D_.push_back(uniform_Q[i]);
        }
      }
  }
  else
  {
      cout << "(Q) More samples requested than available: use whole cloud" << endl;
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

  pcfunctor_.synch3DContent();

  // Compute the diameter of P approximately (randomly). This is far from being
  // Guaranteed close to the diameter but gives good results for most common
  // objects if they are densely sampled.
  P_diameter_ = 0.0;
  for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
    int at = rand() % sampled_Q_3D_.size();
    int bt = rand() % sampled_Q_3D_.size();
    cv::Point3f u(sampled_Q_3D_[bt].x - sampled_Q_3D_[at].x,
                  sampled_Q_3D_[bt].y - sampled_Q_3D_[at].y,
                  sampled_Q_3D_[bt].z - sampled_Q_3D_[at].z);
    double l = cv::norm(u);
    if (l > P_diameter_) {
      P_diameter_ = l;
    }
  }

  // Mean distance and a bit more... We increase the estimation to allow for
  // noise, wrong estimation and non-uniform sampling.
  // can be used to automatically select delta
  //P_mean_distance_ = MeanDistance();

  // Normalize the delta (See the paper) and the maximum base distance.
  // delta = P_mean_distance_ * delta;
  max_base_diameter_ = P_diameter_;  // * estimated_overlap_;

  // RANSAC probability and number of needed trials.
  double first_estimation =
      log(kSmallError) / log(1.0 - pow(options_.overlap_estimation,
                                       static_cast<float>(kMinNumberOfTrials)));
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
}

bool MatchSuper4PCSImpl::ComputeRigidTransformation(const vector< pair<Point3D, Point3D> >& pairs,
                                                    const Eigen::Matrix<Scalar, 3, 1>& centroid1,
                                                    const Eigen::Matrix<Scalar, 3, 1>& centroid2,
                                                    Scalar max_angle,
                                                    Eigen::Matrix<Scalar, 4, 4> &transform,
                                                    Scalar& rms_,
                                                    bool computeScale ) {

  rms_ = kLargeNumber;

  if (pairs.size() == 0 || pairs.size() % 2 != 0)
      return false;


  Scalar kSmallNumber = 1e-6;
  cv::Mat rotation = cv::Mat::eye(3, 3, CV_64F);

  // We only use the first 3 pairs. This simplifies the process considerably
  // because it is the planar case.

  const cv::Point3f& p0 = pairs[0].first;
  const cv::Point3f& p1 = pairs[1].first;
  const cv::Point3f& p2 = pairs[2].first;
  const cv::Point3f& q0 = pairs[0].second;
  const cv::Point3f& q1 = pairs[1].second;
  const cv::Point3f& q2 = pairs[2].second;

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

  rotation = rotate_p.t() * rotate_q;

  // Discard singular solutions. The rotation should be orthogonal.
  cv::Mat unit = rotation * rotation.t();
  if (std::abs(unit.at<double>(0, 0) - 1.0) > kSmallNumber ||
      std::abs(unit.at<double>(1, 1) - 1.0) > kSmallNumber ||
      std::abs(unit.at<double>(2, 2) - 1.0) > kSmallNumber){
      return false;
  }

  // Discard too large solutions (todo: lazy evaluation during boolean computation
  float theta_x = std::abs(atan2(rotation.at<double>(2, 1), rotation.at<double>(2, 2)));
  float theta_y = std::abs(atan2(-rotation.at<double>(2, 0),
                             sqrt(Square(rotation.at<double>(2, 1)) +
                                  Square(rotation.at<double>(2, 2)))));
  float theta_z = std::abs(atan2(rotation.at<double>(1, 0), rotation.at<double>(0, 0)));
  if (theta_x > max_angle ||
      theta_y > max_angle ||
      theta_z > max_angle)
      return false;


  // Compute rms and return it.
  rms_ = Scalar(0.0);
  {
      cv::Mat first(3, 1, CV_64F), transformed;
      for (int i = 0; i < 3; ++i) {
          first.at<double>(0, 0) = pairs[i].second.x - centroid2(0);
          first.at<double>(1, 0) = pairs[i].second.y - centroid2(1);
          first.at<double>(2, 0) = pairs[i].second.z - centroid2(2);
          transformed = rotation * first;
          rms_ += sqrt(Square(transformed.at<double>(0, 0) -
                              (pairs[i].first.x - centroid1(0))) +
                       Square(transformed.at<double>(1, 0) -
                              (pairs[i].first.y - centroid1(1))) +
                       Square(transformed.at<double>(2, 0) -
                              (pairs[i].first.z - centroid1(2))));
      }
  }

  rms_ /= Scalar(pairs.size());

  Eigen::Transform<Scalar, 3, Eigen::Affine> etrans (Eigen::Transform<Scalar, 3, Eigen::Affine>::Identity());

  // Compute scale factor if needed
  if (computeScale){
      std::cerr << __FILE__ << ":" << __LINE__ << " Scale computation not supported yet" << std::endl;
  }

  // compute rotation and translation
  {
      Eigen::Matrix<Scalar, 3, 3> rot;
      cv::cv2eigen(rotation, rot);

      etrans.translate(centroid1);  // translation between quads
      etrans.rotate(rot);           // rotate to align frames
      etrans.translate(-centroid2); // move to congruent quad frame

      transform = etrans.matrix();
  }

  return true;
}


// Performs N RANSAC iterations and compute the best transformation. Also,
// transforms the set Q by this optimal transformation.
bool MatchSuper4PCSImpl::Perform_N_steps(int n, cv::Mat* transformation,
                                    std::vector<Point3D>* Q) {
  if (transformation == NULL || Q == NULL) return false;

#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif

  float last_best_LCP = best_LCP_;
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
    // ok means that we already have the desired LCP.
    if (ok || i > number_of_trials_ || fraction > 0.99) break;
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
        transform_.col(3) = (qcentroid1_ + centroid_P - ( rot * (qcentroid2_ + centroid_Q))).homogeneous();
    }

    cv::eigen2cv( transform_, *transformation );

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
#ifdef TEST_GLOBAL_TIMINGS
    totalTime += Scalar(t.elapsed().count()) / Scalar(1000000.);
#endif

  return ok || current_trial_ >= number_of_trials_;
}

// The main 4PCS function. Computes the best rigid transformation and transfoms
// Q toward P by this transformation.
float MatchSuper4PCSImpl::ComputeTransformation(const std::vector<Point3D>& P,
                                           std::vector<Point3D>* Q,
                                           cv::Mat* transformation) {

  if (Q == nullptr || transformation == nullptr) return kLargeNumber;
  Initialize(P, *Q);

  *transformation = cv::Mat(4, 4, CV_64F, cv::Scalar(0.0));
  for (int i = 0; i < 4; ++i) transformation->at<double>(i, i) = 1.0;
  Perform_N_steps(number_of_trials_, transformation, Q);

#ifdef TEST_GLOBAL_TIMINGS
  cout << "----------- Timings (msec) -------------"           << endl;
  cout << " Total computation time  : " << totalTime           << endl;
  cout << " Total verify time       : " << verifyTime          << endl;
  cout << "    Kdtree query         : " << kdTreeTime          << endl;
#endif

  return best_LCP_;
}

MatchSuper4PCS::MatchSuper4PCS(const Match4PCSOptions& options)
    : pimpl_{new MatchSuper4PCSImpl{options}} {}

MatchSuper4PCS::~MatchSuper4PCS() {}

float
MatchSuper4PCS::ComputeTransformation(const std::vector<Point3D>& P,
                                       std::vector<Point3D>* Q,
                                       cv::Mat* transformation) {
  return pimpl_->ComputeTransformation(P, Q, transformation);
}
}

