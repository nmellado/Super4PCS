// Copyright 2017 Nicolas Mellado
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
// This file is part of the implementation of the 4-points Congruent Sets (4PCS)
// algorithm presented in:
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

#ifndef _SUPER4PCS_ALGO_MATCH_4PCS_BASE_
#define _SUPER4PCS_ALGO_MATCH_4PCS_BASE_

#include <vector>

#ifdef SUPER4PCS_USE_OPENMP
#include <omp.h>
#endif

#include "super4pcs/shared4pcs.h"
#include "super4pcs/sampling.h"
#include "super4pcs/accelerators/kdtree.h"
#include "super4pcs/utils/logger.h"

#ifdef TEST_GLOBAL_TIMINGS
#   include "super4pcs/utils/timer.h"
#endif

namespace GlobalRegistration{


class Match4PCSBase {

public:
    using PairsVector =  std::vector< std::pair<int, int> >;
    using Scalar = typename Point3D::Scalar;
    using VectorType = typename Point3D::VectorType;
    using MatrixType = Eigen::Matrix<Scalar, 4, 4>;
    using LogLevel = Utils::LogLevel;
    struct DummyTransformVisitor {
        inline void operator() (float, float, Eigen::Ref<Match4PCSBase::MatrixType>) const {}
        constexpr bool needsGlobalTransformation() const { return false; }
    };
    using DefaultSampler = Sampling::UniformDistSampler;

    static constexpr int kNumberOfDiameterTrials = 1000;
    static constexpr Scalar kLargeNumber = 1e9;
    static constexpr Scalar distance_factor = 2.0;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual ~Match4PCSBase();

    /// Read access to the sampled clouds used for the registration
    inline const std::vector<Point3D>& getFirstSampled() const {
        return sampled_P_3D_;
    }

    /// Read access to the sampled clouds used for the registration
    inline const std::vector<Point3D>& getSecondSampled() const {
        return sampled_Q_3D_;
    }

    /// Computes an approximation of the best LCP (directional) from Q to P
    /// and the rigid transformation that realizes it. The input sets may or may
    /// not contain normal information for any point.
    /// @param [in] P The first input set.
    /// @param [in] Q The second input set.
    /// as a fraction of the size of P ([0..1]).
    /// @param [out] transformation Rigid transformation matrix (4x4) that brings
    /// Q to the (approximate) optimal LCP. Initial value is considered as a guess
    /// @return the computed LCP measure.
    /// The method updates the coordinates of the second set, Q, applying
    /// the found transformation.
    template <typename Sampler = DefaultSampler,
              typename Visitor = DummyTransformVisitor>
    Scalar
    ComputeTransformation(const std::vector<Point3D>& P,
                          std::vector<Point3D>* Q,
                          Eigen::Ref<MatrixType> transformation,
                          const Sampler& sampler = Sampler(),
                          const Visitor& v = Visitor());


protected:
    /// Number of trials. Every trial picks random base from P.
    int number_of_trials_;
    /// Maximum base diameter. It is computed automatically from the diameter of
    /// P and the estimated overlap and used to limit the distance between the
    /// points in the base in P so that the probability to have all points in
    /// the base as inliers is increased.
    Scalar max_base_diameter_;
    /// The diameter of P.
    Scalar P_diameter_;
    /// Mean distance between points and their nearest neighbor in the set P.
    /// Used to normalize the "delta" which is given in terms of this distance.
    Scalar P_mean_distance_;
    /// The transformation matrix by wich we transform Q to P
    Eigen::Matrix<Scalar, 4, 4> transform_;
    /// Quad centroids in first and second clouds
    /// They are used temporarily and makes the transformations more robust to
    /// noise. At the end, the direct transformation applied as a 4x4 matrix on
    /// every points in Q is computed and returned.
    Eigen::Matrix<Scalar, 3, 1> qcentroid1_, qcentroid2_;
    /// The points in the base (indices to P). It is being updated in every
    /// RANSAC iteration.
    int base_[4];
    /// The current congruent 4 points from Q. Every RANSAC iteration the
    /// algorithm examines a set of such congruent 4-points from Q and retains
    /// the best from them (the one that realizes the best LCP).
    int current_congruent_[4];
    /// Sampled P (3D coordinates).
    std::vector<Point3D> sampled_P_3D_;
    /// Sampled Q (3D coordinates).
    std::vector<Point3D> sampled_Q_3D_;
    /// The 3D points of the base.
    std::vector<Point3D> base_3D_;
    /// The copy of the input Q. We transform Q to match P and returned the
    /// transformed version.
    std::vector<Point3D> Q_copy_;
    /// The centroid of P.
    VectorType centroid_P_;
    /// The centroid of Q.
    VectorType centroid_Q_;
    /// The best LCP (Largest Common Point) fraction so far.
    Scalar best_LCP_;
    /// Current trial.
    int current_trial_;
    /// KdTree used to compute the LCP
    KdTree<Scalar> kd_tree_;
    /// Parameters.
    const Match4PCSOptions options_;

    std::mt19937 randomGenerator_;

    const Utils::Logger &logger_;

#ifdef SUPER4PCS_USE_OPENMP
    /// number of threads used to verify the congruent set
    const int omp_nthread_congruent_;
#endif

#ifdef TEST_GLOBAL_TIMINGS

    Scalar totalTime;
    Scalar kdTreeTime;
    Scalar verifyTime;

    using Timer = GlobalRegistration::Utils::Timer;

#endif

protected:

    Match4PCSBase(const Match4PCSOptions& options
                  , const Utils::Logger &logger
#ifdef SUPER4PCS_USE_OPENMP
                  , const int omp_nthread_congruent = omp_get_max_threads()
#endif
        );

    template <Utils::LogLevel level, typename...Args>
    inline void Log(Args...args) const { logger_.Log<level>(args...); }


    /// Computes the mean distance between points in Q and their nearest neighbor.
    /// We need this for normalization of the user delta (See the paper) to the
    /// "scale" of the set.
    Scalar MeanDistance();


    /// Selects a random triangle in the set P (then we add another point to keep the
    /// base as planar as possible). We apply a simple heuristic that works in most
    /// practical cases. The idea is to accept maximum distance, computed by the
    /// estimated overlap, multiplied by the diameter of P, and try to have
    /// a triangle with all three edges close to this distance. Wide triangles helps
    /// to make the transformation robust while too large triangles makes the
    /// probability of having all points in the inliers small so we try to trade-off.
    bool SelectRandomTriangle(int& base1, int& base2, int& base3);

    /// Takes quadrilateral as a base, computes robust intersection point
    /// (approximate as the lines might not intersect) and returns the invariants
    /// corresponding to the two selected lines. The method also updates the order
    /// of the base base_3D_.
    bool TryQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                          int &base1, int &base2, int &base3, int &base4);


    /// Computes the best rigid transformation between three corresponding pairs.
    /// The transformation is characterized by rotation matrix, translation vector
    /// and a center about which we rotate. The set of pairs is potentially being
    /// updated by the best permutation of the second set. Returns the RMS of the
    /// fit. The method is being called with 4 points but it applies the fit for
    /// only 3 after the best permutation is selected in the second set (see
    /// bellow). This is done because the solution for planar points is much
    /// simpler.
    /// The method is the closed-form solution by Horn:
    /// people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
    bool ComputeRigidTransformation(const std::array<Point3D, 4>& ref,
                                    const std::array<Point3D, 4>& candidate,
                                    const Eigen::Matrix<Scalar, 3, 1>& centroid1,
                                    Eigen::Matrix<Scalar, 3, 1> centroid2,
                                    Scalar max_angle,
                                    Eigen::Ref<MatrixType> transform,
                                    Scalar& rms_,
                                    bool computeScale ) const;

    /// For each randomly picked base, verifies the computed transformation by
    /// computing the number of points that this transformation brings near points
    /// in Q. Returns the current LCP. R is the rotation matrix, (tx,ty,tz) is
    /// the translation vector and (cx,cy,cz) is the center of transformation.template <class MatrixDerived>
    Scalar Verify(const Eigen::Ref<const MatrixType> & mat) const;

    /// Performs n RANSAC iterations, each one of them containing base selection,
    /// finding congruent sets and verification. Returns true if the process can be
    /// terminated (the target LCP was obtained or the maximum number of trials has
    /// been reached), false otherwise.
    template <typename Visitor>
    bool Perform_N_steps(int n,
                         Eigen::Ref<MatrixType> transformation,
                         std::vector<Point3D>* Q,
                         const Visitor& v);

    /// Tries one base and finds the best transformation for this base.
    /// Returns true if the achieved LCP is greater than terminate_threshold_,
    /// else otherwise.
    template <typename Visitor>
    bool TryOneBase(const Visitor &v);

    /// Initializes the data structures and needed values before the match
    /// computation.
    /// @param [in] point_P First input set.
    /// @param [in] point_Q Second input set.
    /// expected to be in the inliers.
    /// This method is called once the internal state of the Base class as been
    /// set.
    virtual void
    Initialize(const std::vector<Point3D>& P,
               const std::vector<Point3D>& Q) = 0;

    template <typename Sampler>
    void init(const std::vector<Point3D>& P,
              const std::vector<Point3D>& Q,
              const Sampler& sampler);

    /// Selects a quadrilateral from P and returns the corresponding invariants
    /// and point indices. Returns true if a quadrilateral has been found, false
    /// otherwise.
    bool SelectQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                             int& base1, int& base2, int& base3, int& base4);

    const std::vector<Point3D>& base3D() const { return base_3D_; }

    /// Constructs pairs of points in Q, corresponding to a single pair in the
    /// in basein P.
    /// @param [in] pair_distance The distance between the pairs in P that we have
    /// to match in the pairs we select from Q.
    /// @param [in] pair_normal_distance The angle between the normals of the pair
    /// in P.
    /// @param [in] pair_distance_epsilon Tolerance on the pair distance. We allow
    /// candidate pair in Q to have distance of
    /// pair_distance+-pair_distance_epsilon.
    /// @param [in] base_point1 The index of the first point in P.
    /// @param [in] base_point2 The index of the second point in P.
    /// @param [out] pairs A set of pairs in Q that match the pair in P with
    /// respect to distance and normals, up to the given tolerance.
    virtual void
    ExtractPairs( Scalar pair_distance,
                  Scalar pair_normals_angle,
                  Scalar pair_distance_epsilon, int base_point1,
                  int base_point2,
                  PairsVector* pairs) const = 0;

    /// Finds congruent candidates in the set Q, given the invariants and threshold
    /// distances. Returns true if a non empty set can be found, false otherwise.
    /// @param invariant1 [in] The first invariant corresponding to the set P_pairs
    /// of pairs, previously extracted from Q.
    /// @param invariant2 [in] The second invariant corresponding to the set
    /// Q_pairs of pairs, previously extracted from Q.
    /// @param [in] distance_threshold1 The distance for verification.
    /// @param [in] distance_threshold2 The distance for matching middle points due
    /// to the invariants (See the paper for e1, e2).
    /// @param [in] P_pairs The first set of pairs.
    /// @param [in] Q_pairs The second set of pairs.
    /// @param [out] quadrilaterals The set of congruent quadrilateral. In fact,
    /// it's a super set from which we extract the real congruent set.
    virtual bool
    FindCongruentQuadrilaterals(Scalar invariant1, Scalar invariant2,
                                Scalar distance_threshold1,
                                Scalar distance_threshold2,
                                const PairsVector& P_pairs,
                                const PairsVector& Q_pairs,
                                std::vector<Quadrilateral>* quadrilaterals) const = 0;

    /// Loop over the set of congruent 4-points and test compatiliby with the
    /// input base.
    /// \param [out] Nb Number of quads corresponding to valid configurations
    template <typename Visitor>
    bool TryCongruentSet(int base_id1,
                         int base_id2,
                         int base_id3,
                         int base_id4,
                         const std::vector<Quadrilateral> &congruent_quads,
                         const Visitor &v,
                         size_t &nbCongruent);
private:
    void initKdTree();

}; /// class Match4PCSBase
} /// namespace Super4PCS

#include "super4pcs/algorithms/match4pcsBase.hpp"

#endif
