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
// Authors: Dror Aiger, Yoni Weill
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

#ifndef _4PCS_H_
#define _4PCS_H_

#include "shared4pcs.h"

namespace match_4pcs {
class Match4PCSImpl;

// Class for the computation of the 4PCS algorithm.
class Match4PCS {
 public:
  explicit Match4PCS(const Match4PCSOptions& options);

  ~Match4PCS();

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
  std::unique_ptr<Match4PCSImpl> pimpl_;
};
} // namespace match_4pcs

#endif  // _4PCS_H_
