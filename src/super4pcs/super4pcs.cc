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
// Technically, we map the algorithm as an 'instance problem' and solve it
// efficiently using a smart indexing data organization. The algorithm is
// simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in
// significant speedup over alternative approaches and allows unstructured
// efficient acquisition of scenes at scales previously not possible. Complete
// source code and datasets are available for research use at
// http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.

#include "4pcs.h"
#include "match4pcsBase.h"

#include "Eigen/Core"
#include "Eigen/Geometry"                 // MatrixBase.homogeneous()
#include "Eigen/SVD"                      // Transform.computeRotationScaling()

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>

#include "accelerators/pairExtraction/bruteForceFunctor.h"
#include "accelerators/pairExtraction/intersectionFunctor.h"
#include "accelerators/pairExtraction/intersectionPrimitive.h"
#include "accelerators/normalset.h"
#include "accelerators/normalHealSet.h"
#include "accelerators/bbox.h"

#include "pairCreationFunctor.h"

#include <fstream>
#include <array>
#include <time.h>

#define sqr(x) ((x)*(x))
#define norm2(p) (sqr(p.x)+sqr(p.y)+sqr(p.z))

#ifdef TEST_GLOBAL_TIMINGS
#   include "utils/timer.h"
#endif

//#define MULTISCALE

namespace match_4pcs {

using namespace std;


class MatchSuper4PCSImpl  : public Super4PCS::Match4PCSBase {
public:
    using Base = Super4PCS::Match4PCSBase;
    using Scalar = Base::Scalar;

    typedef PairCreationFunctor<Scalar>::PairsVector PairsVector;

#ifdef TEST_GLOBAL_TIMINGS

    Scalar totalTime;
    Scalar kdTreeTime;
    Scalar verifyTime;

    using Timer = Super4PCS::Utils::Timer;

#endif
 public:
  explicit MatchSuper4PCSImpl(const Match4PCSOptions& options)
      : Base(options),
        pcfunctor_(options_, sampled_Q_3D_) { }

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

  // Read access to the sampled clouds used for the registration
  inline const std::vector<Point3D>& getFirstSampled() const {
      return sampled_P_3D_;
  }

  // Read access to the sampled clouds used for the registration
  inline const std::vector<Point3D>& getSecondSampled() const {
      return sampled_Q_3D_;
  }

 private:
  // Private data contains parameters and internal variables that are computed
  // and change during the match computation. All parameters have default
  // values.

  // Internal data members.

  PairCreationFunctor<Scalar> pcfunctor_;

protected:
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
  ExtractPairs(Scalar pair_distance, Scalar pair_normals_angle,
                       Scalar pair_distance_epsilon, int base_point1,
                       int base_point2,
                       PairsVector* pairs) override;

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
  bool FindCongruentQuadrilaterals(Scalar invariant1, Scalar invariant2,
                                   Scalar distance_threshold1,
                                   Scalar distance_threshold2,
                                   const PairsVector& P_pairs,
                                   const PairsVector& Q_pairs,
                                   std::vector<Super4PCS::Quadrilateral>* quadrilaterals) override;

private:
  // Initializes the data structures and needed values before the match
  // computation.
  // @param [in] point_P First input set.
  // @param [in] point_Q Second input set.
  // expected to be in the inliers.
  void Initialize(const std::vector<Point3D>& P, const std::vector<Point3D>& Q);

};

// Finds congruent candidates in the set Q, given the invariants and threshold
// distances.
bool MatchSuper4PCSImpl::FindCongruentQuadrilaterals(
    Scalar invariant1, Scalar invariant2, Scalar distance_threshold1,
    Scalar distance_threshold2, const std::vector<std::pair<int, int>>& P_pairs,
    const std::vector<std::pair<int, int>>& Q_pairs,
    std::vector<Super4PCS::Quadrilateral>* quadrilaterals) {


  if (quadrilaterals == NULL) return false;
  quadrilaterals->clear();

  int number_of_points = 2 * P_pairs.size();

  // Compute the angle formed by the two vectors of the basis
  Point3D b1 = base_3D_[1] - base_3D_[0];  b1.normalize();
  Point3D b2 = base_3D_[3] - base_3D_[2];  b2.normalize();
  Scalar alpha = /*std::abs(*/b1.dot(b2)/*)*/;

  // 1. Datastructure construction
  typedef PairCreationFunctor<Scalar>::Point Point;
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

  for (unsigned int i = 0; i < P_pairs.size(); ++i) {
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
  for (unsigned int i = 0; i < Q_pairs.size(); ++i) {
    const Point& p1 = pcfunctor_.points[Q_pairs[i].first];
    const Point& p2 = pcfunctor_.points[Q_pairs[i].second];

    const Point3D& pq1 = sampled_Q_3D_[Q_pairs[i].first];
    const Point3D& pq2 = sampled_Q_3D_[Q_pairs[i].second];

    nei.clear();

    Point   query  =  p1 + invariant2 * ( p2 - p1 );
    Point3D queryQ = pq1 + invariant2 * (pq2 - pq1);

    Point queryn = (p2 - p1).normalized();

//    cout << "query: " << endl
//         << p1.transpose() << "(" << Q_pairs[i].first  << ")" << endl
//         << p2.transpose() << "(" << Q_pairs[i].second << ")" << endl
//         << query.transpose() << endl
//         << queryn.transpose() << endl;

    nset.getNeighbors( query,  queryn, alpha, nei);


    Point3D invPoint;
    //const float distance_threshold2s = distance_threshold2 * distance_threshold2;
    for (unsigned int k = 0; k != nei.size(); k++){
      int id = nei[k];

      const Point3D& pp1 = sampled_Q_3D_[P_pairs[id].first];
      const Point3D& pp2 = sampled_Q_3D_[P_pairs[id].second];

      invPoint = pp1 + (pp2 - pp1) * invariant1;

       // use also distance_threshold2 for inv 1 and 2 in 4PCS
      if (cv::norm(queryQ-invPoint) <= distance_threshold2){
          comb.insert(std::make_pair(id, i));
      }
    }
  }

  for (std::set< std::pair<unsigned int, unsigned int > >::const_iterator it =
             comb.cbegin();
       it != comb.cend(); it++){
    const unsigned int & id = (*it).first;
    const unsigned int & i  = (*it).second;

    quadrilaterals->push_back(
                Super4PCS::Quadrilateral(P_pairs[id].first, P_pairs[id].second,
                              Q_pairs[i].first,  Q_pairs[i].second));
  }

  return quadrilaterals->size() != 0;
}


// Constructs two sets of pairs in Q, each corresponds to one pair in the base
// in P, by having the same distance (up to some tolerantz) and optionally the
// same angle between normals and same color.
void
MatchSuper4PCSImpl::ExtractPairs(Scalar pair_distance,
                                    Scalar pair_normals_angle,
                                    Scalar pair_distance_epsilon,
                                    int base_point1, int base_point2,
                                    PairsVector* pairs) {

  using namespace Super4PCS::Accelerators::PairExtraction;

  pcfunctor_.pairs = pairs;

  pairs->clear();
  pairs->reserve(2 * pcfunctor_.points.size());

  pcfunctor_.pair_distance         = pair_distance;
  pcfunctor_.pair_distance_epsilon = pair_distance_epsilon;
  pcfunctor_.pair_normals_angle    = pair_normals_angle;
  pcfunctor_.norm_threshold =
      0.5 * options_.max_normal_difference * M_PI / 180.0;

  pcfunctor_.setRadius(pair_distance);
  pcfunctor_.setBase(base_point1, base_point2, base_3D_);


#ifdef MULTISCALE
  BruteForceFunctor
  <PairCreationFunctor<Scalar>::Point, 3, Scalar> interFunctor;
#else
  IntersectionFunctor
          <PairCreationFunctor<Scalar>::Primitive,
          PairCreationFunctor<Scalar>::Point, 3, Scalar> interFunctor;
#endif

  Scalar eps = pcfunctor_.getNormalizedEpsilon(pair_distance_epsilon);

  interFunctor.process(pcfunctor_.primitives,
                       pcfunctor_.points,
                       eps,
                       50,
                       pcfunctor_);
}




// Initialize all internal data structures and data members.
void MatchSuper4PCSImpl::Initialize(const std::vector<Point3D>& P,
                               const std::vector<Point3D>& Q) {
    Base::init(P,Q);

#ifdef TEST_GLOBAL_TIMINGS
    kdTreeTime = 0;
    totalTime  = 0;
    verifyTime = 0;
#endif

  pcfunctor_.synch3DContent();

  best_LCP_ = Verify(transform_);
  printf("Initial LCP: %f\n", best_LCP_);

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

const std::vector<Point3D>&
MatchSuper4PCS::getFirstSampled() const{
  return pimpl_->getFirstSampled();
}

const std::vector<Point3D>&
MatchSuper4PCS::getSecondSampled() const{
  return pimpl_->getSecondSampled();
}

}

