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

#include "super4pcs/algorithms/super4pcs.h"
#include "super4pcs/accelerators/bbox.h"

#ifdef SUPER4PCS_USE_CHEALPIX
#include "super4pcs/accelerators/normalHealSet.h"
#else
#include "super4pcs/accelerators/normalset.h"
#endif

#include <fstream>
#include <array>
#include <time.h>


//#define MULTISCALE

namespace GlobalRegistration {


MatchSuper4PCS::MatchSuper4PCS(const Match4PCSOptions& options,
                               const Utils::Logger &logger)
    : Base(options
           , logger
#ifdef SUPER4PCS_USE_OPENMP
           , 1
#endif
           ),
      pcfunctor_(options_, sampled_Q_3D_) { }

MatchSuper4PCS::~MatchSuper4PCS() { }

// Finds congruent candidates in the set Q, given the invariants and threshold
// distances.
bool
MatchSuper4PCS::FindCongruentQuadrilaterals(
        Scalar invariant1,
        Scalar invariant2,
        Scalar /*distance_threshold1*/,
        Scalar distance_threshold2,
        const std::vector<std::pair<int, int>>& P_pairs,
        const std::vector<std::pair<int, int>>& Q_pairs,
        std::vector<Quadrilateral>* quadrilaterals) const {

    typedef PairCreationFunctor<Scalar>::Point Point;

#ifdef SUPER4PCS_USE_CHEALPIX
    typedef GlobalRegistration::IndexedNormalHealSet IndexedNormalSet3D;
#else
    typedef  GlobalRegistration::IndexedNormalSet
                    < Point,   //! \brief Point type used internally
                      3,       //! \brief Nb dimension
                      7,       //! \brief Nb cells/dim normal
                      Scalar>  //! \brief Scalar type
    IndexedNormalSet3D;
#endif


  if (quadrilaterals == nullptr) return false;

  quadrilaterals->clear();

  // Compute the angle formed by the two vectors of the basis
  const Scalar alpha =
          (base_3D_[1].pos() - base_3D_[0].pos()).normalized().dot(
          (base_3D_[3].pos() - base_3D_[2].pos()).normalized());

  // 1. Datastructure construction
  const Scalar eps = pcfunctor_.getNormalizedEpsilon(distance_threshold2);

  IndexedNormalSet3D nset (eps);

  for (size_t i = 0; i <  P_pairs.size(); ++i) {
    const Point& p1 = pcfunctor_.points[P_pairs[i].first];
    const Point& p2 = pcfunctor_.points[P_pairs[i].second];
    const Point  n  = (p2 - p1).normalized();

    nset.addElement((p1+ Point::Scalar(invariant1) * (p2 - p1)).eval(), n, i);
  }


  std::set< std::pair<unsigned int, unsigned int > > comb;

  unsigned int j = 0;
  std::vector<unsigned int> nei;
  // 2. Query time
  for (unsigned int i = 0; i < Q_pairs.size(); ++i) {
    const Point& p1 = pcfunctor_.points[Q_pairs[i].first];
    const Point& p2 = pcfunctor_.points[Q_pairs[i].second];

    const VectorType& pq1 = sampled_Q_3D_[Q_pairs[i].first].pos();
    const VectorType& pq2 = sampled_Q_3D_[Q_pairs[i].second].pos();

    nei.clear();

    const Point      query  =  p1 + invariant2 * ( p2 - p1 );
    const VectorType queryQ = pq1 + invariant2 * (pq2 - pq1);

    const Point queryn = (p2 - p1).normalized();

    nset.getNeighbors( query, queryn, alpha, nei);


    VectorType invPoint;
    //const Scalar distance_threshold2s = distance_threshold2 * distance_threshold2;
    for (unsigned int k = 0; k != nei.size(); k++){
      const int id = nei[k];

      const VectorType& pp1 = sampled_Q_3D_[P_pairs[id].first].pos();
      const VectorType& pp2 = sampled_Q_3D_[P_pairs[id].second].pos();

      invPoint = pp1 + (pp2 - pp1) * invariant1;

       // use also distance_threshold2 for inv 1 and 2 in 4PCS
      if ((queryQ-invPoint).squaredNorm() <= distance_threshold2){
          comb.emplace(id, i);
      }
    }
  }

  for (std::set< std::pair<unsigned int, unsigned int > >::const_iterator it =
             comb.cbegin();
       it != comb.cend(); it++){
    const unsigned int & id = (*it).first;
    const unsigned int & i  = (*it).second;

    quadrilaterals->emplace_back(P_pairs[id].first, P_pairs[id].second,
                                 Q_pairs[i].first,  Q_pairs[i].second);
  }

  return quadrilaterals->size() != 0;
}


// Constructs two sets of pairs in Q, each corresponds to one pair in the base
// in P, by having the same distance (up to some tolerantz) and optionally the
// same angle between normals and same color.
void
MatchSuper4PCS::ExtractPairs(Scalar pair_distance,
                             Scalar pair_normals_angle,
                             Scalar pair_distance_epsilon,
                             int base_point1,
                             int base_point2,
                             PairsVector* pairs) const {

  using namespace GlobalRegistration::Accelerators::PairExtraction;

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
  <PairCreationFunctor<Scalar>::Primitive, PairCreationFunctor<Scalar>::Point, 3, Scalar> interFunctor;
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
void
MatchSuper4PCS::Initialize(const std::vector<Point3D>& /*P*/,
                           const std::vector<Point3D>& /*Q*/) {
  pcfunctor_.synch3DContent();
}


} // namespace Super4PCS
