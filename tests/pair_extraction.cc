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
// Authors: Nicolas Mellado
//
// This test check the validity of the pair extraction subroutine on random data
// of different dimensions (2,3 and 4).
//
//
// This test is part of the implementation of the Super 4-points Congruent Sets
// (Super 4PCS) algorithm presented in:
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

#include "Eigen/Dense"

#include <fstream>
#include <iostream>
#include <string>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>
#include "accelerators/pairExtraction/bruteForceFunctor.h"
#include "accelerators/pairExtraction/intersectionFunctor.h"
#include "accelerators/pairExtraction/intersectionPrimitive.h"
#include "utils/timer.h"
#include "bbox.h"

#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <utility> // pair

#include "testing.h"

//#include "normalset.h"
//#include "normalHealSet.h"

using namespace match_4pcs;

struct PairCreationFunctor{
  typedef std::pair<unsigned int, unsigned int>ResPair;
  std::vector< ResPair >pairs;

  std::vector<unsigned int> ids;

  inline void beginPrimitiveCollect(int primId){
  }

  inline void endPrimitiveCollect(int primId){
  }
  inline void process(int primId, int pointId){
    //if(pointId >= 10)
      if (primId>pointId)
        pairs.push_back(ResPair(pointId, primId));
  }
};


/*!
 * \brief Generate a set of random points and spheres, and test the pair extraction

   Two tests are operated here:
    - Check the validity of the sphere to point intersection test
    - Check the rendering

   Note here that rendering timings are not optimal because we have a volume
   uniformly sampled and not a surface.

 */
template<typename Scalar, typename Point, typename Primitive, typename Functor>
void testFunction( Scalar r, Scalar epsilon,
                   unsigned int nbPoints,
                   unsigned int minNodeSize){

  // Init required structures
  Super4PCS::Utils::Timer t;
  std::vector< std::pair<unsigned int, unsigned int> > p2;
  p2.reserve(nbPoints*nbPoints);

  PairCreationFunctor functor;
  functor.ids.clear();
  for(unsigned int i = 0; i < nbPoints; i++)
    functor.ids.push_back(i);
  functor.pairs.reserve(nbPoints*nbPoints);

  //std::cout << "**************************" << std::endl;
  //std::cout << "Epsilon = " << epsilon << std::endl;

  // Init Random Positions
  std::vector<Point> points;
  std::vector<Primitive> primitives;

  //std::cout << "Points: " << std::endl;
  Point half (Point::Ones()/2.f);
  for(unsigned int i = 0; i != nbPoints; i++){
    // generate points on a sphere
    Point p (0.5f*Point::Random().normalized() + half);
    points.push_back(p);
    primitives.push_back(Primitive(p, r));
    //std::cout << p.transpose() << std::endl;
  }

  // Test test intersection procedure
  // Here we compare the brute force pair extraction and the sphere to point
  // intersection procedure
  {
    Primitive& sphere = primitives.front();
    for(unsigned int i = 0; i != nbPoints; i++){
      const Point&p = points[i];
      VERIFY( sphere.intersectPoint(p, epsilon) ==
              SQR((p - sphere.center()).norm()- sphere.radius()) < SQR(epsilon));
    }
  }


  // Test Rendering process
  {
    Functor IF;

    // Extract pairs using rendering process
    t.reset();
    IF.process(primitives, points, epsilon, minNodeSize, functor);
    const auto IFtimestep = t.elapsed();


    // Extract pairs using brute force
    t.reset();
    for(unsigned int i = 0; i != nbPoints; i++)
      for(unsigned int j = i+1; j < nbPoints; j++)
         if (primitives[j].intersectPoint(points[i], epsilon))
          p2.push_back(std::make_pair(i,j));
    const auto BFtimestep = t.elapsed();

    std::cout << "Timers (" << (IFtimestep.count() < BFtimestep.count()
                                ? "PASSED" : "NOT PASSED")
              << "): \t Functor: " << IFtimestep.count()/1000
              << "\t BruteForce: " << BFtimestep.count()/1000 << std::endl;

    // Check we get the same set size
    std::cout << "Size check (" << (functor.pairs.size() == p2.size()
                                ? "PASSED" : "NOT PASSED")
              << "): \t Functor: " << functor.pairs.size()
              << " \t BruteForce: " << p2.size() << std::endl;

    // sort to ensure containers consistency
    std::sort(functor.pairs.begin(), functor.pairs.end());
    std::sort(p2.begin(), p2.end());


//    std::cout << "Functor: " << std::endl;
//    for (const auto&p0 : functor.pairs)
//        std::cout << "\t" << p0.first << " - " << p0.second << std::endl;
//    std::cout << "Brute Force: " << std::endl;
//    for (const auto&p0 : p2)
//        std::cout << "\t" << p0.first << " - " << p0.second << std::endl;


    VERIFY( functor.pairs.size() == p2.size() );
    VERIFY( std::equal(functor.pairs.begin(), functor.pairs.end(), p2.begin()));

  }
}


template<typename Scalar,
         int Dim,
         template <typename,typename,int,typename> class _Functor>
void callSubTests()
{
    using namespace Super4PCS::Accelerators::PairExtraction;

    typedef  Eigen::Matrix<Scalar, Dim, 1> EigenPoint;
    typedef  HyperSphere< EigenPoint, Dim, Scalar > Sphere;
    typedef _Functor<Sphere, EigenPoint, Dim, Scalar> Functor;

    Scalar   r = 0.5; // radius of the spheres
    Scalar eps = GetRoundedEpsilonValue(0.125/16.); // epsilon value
    unsigned int nbPoint = 5000;  // size of Q point cloud
    int minNodeSize = 100;

    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Scalar,
                                    EigenPoint,
                                    Sphere,
                                    Functor>(r, eps, nbPoint, minNodeSize) ));
    }
}

int main(int argc, const char **argv) {
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    using std::cout;
    using std::endl;
    using namespace Super4PCS::Accelerators::PairExtraction;

    cout << "Extract pairs in 2 dimensions (BRUTE FORCE)..." << endl;
    callSubTests<float, 2, BruteForceFunctor>();
    callSubTests<double, 2, BruteForceFunctor>();
    callSubTests<long double, 2, BruteForceFunctor>();
    cout << "Ok..." << endl;

    cout << "Extract pairs in 2 dimensions (RENDERING)..." << endl;
    callSubTests<float, 2, IntersectionFunctor>();
    callSubTests<double, 2, IntersectionFunctor>();
    callSubTests<long double, 2, IntersectionFunctor>();
    cout << "Ok..." << endl;

    cout << "Extract pairs in 3 dimensions (BRUTE FORCE)..." << endl;
    callSubTests<float, 3, BruteForceFunctor>();
    callSubTests<double, 3, BruteForceFunctor>();
    callSubTests<long double, 3, BruteForceFunctor>();
    cout << "Ok..." << endl;

    cout << "Extract pairs in 3 dimensions (RENDERING)..." << endl;
    callSubTests<float, 3, IntersectionFunctor>();
    callSubTests<double, 3, IntersectionFunctor>();
    callSubTests<long double, 3, IntersectionFunctor>();
    cout << "Ok..." << endl;

    cout << "Extract pairs in 4 dimensions (BRUTE FORCE)..." << endl;
    callSubTests<float, 4, BruteForceFunctor>();
    callSubTests<double, 4, BruteForceFunctor>();
    callSubTests<long double, 4, BruteForceFunctor>();
    cout << "Ok..." << endl;

    cout << "Extract pairs in 4 dimensions (RENDERING)..." << endl;
    callSubTests<float, 4, IntersectionFunctor>();
    callSubTests<double, 4, IntersectionFunctor>();
    callSubTests<long double, 4, IntersectionFunctor>();
    cout << "Ok..." << endl;

    return EXIT_SUCCESS;
}
