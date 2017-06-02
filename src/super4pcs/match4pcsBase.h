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

#ifndef _MATCH_4PCS_BASE_
#define _MATCH_4PCS_BASE_

#include <vector>
#include <opencv2/highgui/highgui.hpp>

#include "shared4pcs.h"
#include "sampling.h"
#include "accelerators/kdtree.h"

namespace Super4PCS{

// Holds a base from P. The base contains 4 points (indices) from the set P.
struct Quadrilateral {
    std::array <int, 4> vertices;
    inline Quadrilateral(int vertex0, int vertex1, int vertex2, int vertex3) {
        vertices = { vertex0, vertex1, vertex2, vertex3 };
    }
};


class Match4PCSBase {

public:
    using Point3D = match_4pcs::Point3D;
    using PairsVector =  std::vector< std::pair<int, int> >;
    using Scalar = double;

    static constexpr int kNumberOfDiameterTrials = 1000;


protected:
    // Number of trials. Every trial picks random base from P.
    int number_of_trials_;
    // Maximum base diameter. It is computed automatically from the diameter of
    // P and the estimated overlap and used to limit the distance between the
    // points in the base in P so that the probability to have all points in
    // the base as inliers is increased.
    float max_base_diameter_;
    // The diameter of P.
    float P_diameter_;
    // Mean distance between points and their nearest neighbor in the set P.
    // Used to normalize the "delta" which is given in terms of this distance.
    float P_mean_distance_;
    // The centroid about which we rotate a congruent set in Q to match the base
    // in P. It is used temporarily and makes the transformations more robust to
    // noise. At the end, the direct transformation applied as a 4x4 matrix on
    // every points in Q is computed and returned.
    cv::Point3f centroid_;
    // The translation vector by which we move Q to match P.
    cv::Point3f translate_;
    // The rotation matrix by which we rotate Q to match P.
    cv::Mat rotate_;
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
    // KdTree used to compute the LCP
    Super4PCS::KdTree<Scalar> kd_tree_;
    // Parameters.
    match_4pcs::Match4PCSOptions options_;


    Match4PCSBase(const match_4pcs::Match4PCSOptions& options);

protected:
    // Computes the mean distance between points in Q and their nearest neighbor.
    // We need this for normalization of the user delta (See the paper) to the
    // "scale" of the set.
    Scalar MeanDistance();

    void init(const std::vector<Point3D>& P,
                     const std::vector<Point3D>& Q);


    // Selects a random triangle in the set P (then we add another point to keep the
    // base as planar as possible). We apply a simple heuristic that works in most
    // practical cases. The idea is to accept maximum distance, computed by the
    // estimated overlap, multiplied by the diameter of P, and try to have
    // a triangle with all three edges close to this distance. Wide triangles helps
    // to make the transformation robust while too large triangles makes the
    // probability of having all points in the inliers small so we try to trade-off.
    bool SelectRandomTriangle(int* base1, int* base2, int* base3);

private:
    void initKdTree();

}; // class Match4PCSBase
} // namespace Super4PCS

#endif
