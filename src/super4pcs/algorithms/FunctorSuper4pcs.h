//
// Created by shiro on 24/04/18.
//

#ifndef SUPER4PCS_FUNCTORSUPER4PCS_H
#define SUPER4PCS_FUNCTORSUPER4PCS_H


#include <vector>
#include "../shared4pcs.h"
#include "../algorithms/pairCreationFunctor.h"
#include "../accelerators/bbox.h"

#ifdef SUPER4PCS_USE_CHEALPIX
#include "super4pcs/accelerators/normalHealSet.h"
#else
#include "../accelerators/normalset.h"
#include "../accelerators/utils.h"

#endif

#include <fstream>
#include <array>
#include <time.h>

namespace GlobalRegistration {
    struct MatchSuper4PCS {
    public :
        using TypeBase = std::vector<Point3D>;
        using Scalar      = typename Point3D::Scalar;
        using PairsVector = std::vector< std::pair<int, int> >;
        using VectorType  = typename Point3D::VectorType;
        using OptionType  = Match4PCSOptions;


    private :
        Match4PCSOptions myOptions_;
        std::vector<Point3D> mySampled_Q_3D_;
        TypeBase myBase_3D_;

        /// Private data contains parameters and internal variables that are computed
        /// and change during the match computation. All parameters have default
        /// values.

        /// Internal data members.

        mutable PairCreationFunctor<Scalar> pcfunctor_;

    public :
        // Initialize all internal data structures and data members.
        inline void Initialize(const std::vector<Point3D>& /*P*/,
                                   const std::vector<Point3D>& /*Q*/) {
            pcfunctor_.synch3DContent();
        }

        inline void setOptions (Match4PCSOptions options) {
            myOptions_ = options;
        }

        inline void setSampled_Q_3D (std::vector<Point3D> sampled_Q_3D_) {
            mySampled_Q_3D_ = sampled_Q_3D_;
        }

        inline void setBase_3D (TypeBase base_3D_) {
            myBase_3D_ = base_3D_;
        }

        // Constructs two sets of pairs in Q, each corresponds to one pair in the base
        // in P, by having the same distance (up to some tolerantz) and optionally the
        // same angle between normals and same color.
        inline void ExtractPairs(Scalar pair_distance,
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
                    0.5 * myOptions_.max_normal_difference * M_PI / 180.0;

            pcfunctor_.setRadius(pair_distance);
            pcfunctor_.setBase(base_point1, base_point2, myBase_3D_);


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

// Finds congruent candidates in the set Q, given the invariants and threshold
// distances.
        inline bool FindCongruentQuadrilaterals(
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
                    (myBase_3D_[1].pos() - myBase_3D_[0].pos()).normalized().dot(
                            (myBase_3D_[3].pos() - myBase_3D_[2].pos()).normalized());

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

                const VectorType& pq1 = mySampled_Q_3D_[Q_pairs[i].first].pos();
                const VectorType& pq2 = mySampled_Q_3D_[Q_pairs[i].second].pos();

                nei.clear();

                const Point      query  =  p1 + invariant2 * ( p2 - p1 );
                const VectorType queryQ = pq1 + invariant2 * (pq2 - pq1);

                const Point queryn = (p2 - p1).normalized();

                nset.getNeighbors( query, queryn, alpha, nei);


                VectorType invPoint;
                //const Scalar distance_threshold2s = distance_threshold2 * distance_threshold2;
                for (unsigned int k = 0; k != nei.size(); k++){
                    const int id = nei[k];

                    const VectorType& pp1 = mySampled_Q_3D_[P_pairs[id].first].pos();
                    const VectorType& pp2 = mySampled_Q_3D_[P_pairs[id].second].pos();

                    invPoint = pp1 + (pp2 - pp1) * invariant1;

                    // use also distance_threshold2 for inv 1 and 2 in 4PCS
                    if ((queryQ-invPoint).squaredNorm() <= distance_threshold2){
                        comb.emplace(id, i);
                    }
                }
            }

            for (std::set< std::pair<unsigned int, unsigned int>>::const_iterator it = comb.cbegin();
                    it != comb.cend(); it++) {
                const unsigned int & id = (*it).first;
                const unsigned int & i  = (*it).second;

                quadrilaterals->emplace_back(P_pairs[id].first, P_pairs[id].second,
                                             Q_pairs[i].first,  Q_pairs[i].second);
            }

            return quadrilaterals->size() != 0;
        }

    };
}

#endif //SUPER4PCS_FUNCTORSUPER4PCS_H
