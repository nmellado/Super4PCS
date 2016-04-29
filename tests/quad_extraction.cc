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
// This test check the validity of the quad extraction subroutine on random data
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
#include "accelerators/normalHealSet.h"
#include "accelerators/normalset.h"
#include "accelerators/pairExtraction/bruteForceFunctor.h"
#include "bbox.h"

#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <utility> // pair

#include "testing.h"

//#include "normalset.h"
//#include "normalHealSet.h"

using namespace match_4pcs;

namespace Utilities{

    static double Timer_getTime(void)
    {
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);
        return (double)tv.tv_sec + 1.e-6 * (double)tv.tv_usec;
    }

    typedef double Timer;

    static inline void startTimer(Timer& timer) { timer = Timer_getTime(); }
    static inline void stopTimer(Timer& timer) { timer = Timer_getTime() - timer; }
}

/*!
 * \brief Generate a set of random points and spheres, and test the pair extraction
 
   Two tests are operated here:
    - Check the validity of the sphere to point intersection test
    - Check the rendering 
    
   Note here that rendering timings are not optimal because we have a volume 
   uniformly sampled and not a surface.
   
 */
template<typename Scalar,
         int Dim,
         typename Point,
         typename NormalSet>
void testFunction( Scalar epsilon,
                   unsigned int nbPoints){

  using PairFunctor = Super4PCS::Accelerators::PairExtraction::
    BruteForceFunctor<Point, Dim, Scalar>;

  Utilities::Timer t; 

  auto getRandomPoint = [](Scalar factor = Scalar(1)){
      static const Point half (Point::Ones()/Scalar(2));
      return Point (Scalar(0.5 * factor)*Point::Random() + half);
  };



  // Init random basis
  struct Basis{
      Point a, b, c, d;
      Scalar inv1, inv2;
  } ;
  Basis basis;

  // get random positions in [0:0.5]^n
  Scalar factor (0.2);
  basis.a = getRandomPoint(factor);
  basis.b = getRandomPoint(factor);
  basis.inv1 = factor*Scalar(rand() % 10000) / Scalar(10000);
  basis.inv2 = factor*Scalar(rand() % 10000) / Scalar(10000);
  Point tmp = getRandomPoint(factor);
  Point pivot = basis.a + (basis.b-basis.a) * basis.inv1;
  basis.c = pivot - tmp * basis.inv2;
  basis.d = pivot + tmp * (Scalar(1) - basis.inv2);

  std::cout << "Basis: "
            << "\n  " << basis.a.transpose()
            << "\n  " << basis.b.transpose()
            << "\n  " << basis.c.transpose()
            << "\n  " << basis.d.transpose() << std::endl;

  NormalSet _set(epsilon);

  struct IndexedPoint{
      Point p, n;
      int id;
  };


  // generate random point cloud
  std::vector< Point > points;
  points.reserve(nbPoints);

  for(unsigned int i = 0; i != nbPoints; i++){
      points.push_back(Point::Random());
  }

  struct ExtractFunctor{
      Basis basis;
      std::vector< std::pair<int,int> >pairs1;
      std::vector< std::pair<int,int> >pairs2;
      std::vector< Point > &points;
      Scalar epsilon;
      inline ExtractFunctor(std::vector< Point > &p) : points(p) {}

      inline
      void process(int i, int j){
          if (i>j){
              if (std::abs(  (points[j] - points[i]).norm()
                           - (basis.b   - basis.a).norm())
                      < epsilon)
                  pairs1.push_back(std::pair<int,int>(i,j));
              if (std::abs(  (points[j] - points[i]).norm()
                           - (basis.d   - basis.c).norm())
                      < epsilon)
                  pairs2.push_back(std::pair<int,int>(i,j));
          }
      }
  };
  ExtractFunctor extractfunctor (points);
  extractfunctor.basis  = basis;
  extractfunctor.epsilon = epsilon;

  //extract pairs
  PairFunctor pfunctor;
  pfunctor.process(points, points, epsilon, 0, extractfunctor);
/*      IndexedPoint p;

      p.p  = getRandomPoint();
      p.n  = getRandomPoint().normalized();
      p.id = i;

      points.push_back(p);
      _set.addElement(p.p,p.n, p.id)*/;
 //

  //////////////////////////////////////////////////////////////////////////////
  /// quadrangles extraction

  // Compute the angle formed by the two vectors of the basis
  Point3D b1 ((basis.b - basis.a).eval());  b1.normalize();
  Point3D b2 ((basis.d - basis.c).eval());  b2.normalize();
  double alpha = b1.dot(b2);


  struct Quadrangle { std::array<int, 4> ids; };
  std::vector< Quadrangle > gt_quads;
  for(unsigned int i_a = 0; i_a < nbPoints; ++i_a){
//      const Point& p_a = points.
//      for(unsigned int i_b = i_a+1; i_b < nbPoints; ++i_b){
//          for(unsigned int i_c = 0; i_c < nbPoints; ++i_c){
//              for(unsigned int i_d = i_c+1; i_d < nbPoints; ++i_d){

//              }
//          }
//      }
  }

  
//  // Init random Spheres
//  std::vector<Primitive> primitives;
//  for(unsigned int i = 0; i != nbPrimitives; i++){
//    Point p (0.5f*Point::Random() + half);

//    primitives.push_back(Primitive(p, r));
//  }
  
//  // Test test intersection procedure
//  // Here we compare the brute force pair extraction and the sphere to point
//  // intersection procedure
//  {
//    Primitive& sphere = primitives.front();
//    for(unsigned int i = 0; i != nbPoints; i++){
//      const Point&p = points[i];
//      VERIFY( sphere.intersectPoint(p, epsilon) ==
//              SQR((p - sphere.center()).norm()- sphere.radius()) < SQR(epsilon));
//    }
//  }
  
    
//  // Test Rendering process
//  {
//    Functor IF;
    
//    // Extract pairs using rendering process
//    Utilities::startTimer(t);
//    IF.process(primitives, points, epsilon, 20, functor);
//    Utilities::stopTimer(t);
           
     
//    // Extract pairs using brute force
//    Utilities::startTimer(t);
//    for(unsigned int i = 0; i != nbPoints; i++)
//      for(unsigned int j = i+1; j < primitives.size(); j++)
//         if (primitives[j].intersectPoint(points[i], epsilon))
//          p2.push_back(std::pair<unsigned int, unsigned int>(i,j));
//    Utilities::stopTimer(t);
    
//    // Check we get the same set size
//    VERIFY( functor.pairs.size() == p2.size());
//  }
}


template<typename Scalar, 
         int Dim>
void callSubTests()
{
    typedef  Eigen::Matrix<Scalar, Dim, 1> EigenPoint;

    Scalar eps = 0.125/8.; // epsilon value
    
    unsigned int nbPoint = 1000;  // size of Q point cloud

    using NormalSet0 = Super4PCS::IndexedNormalHealSet;
    using NormalSet1 = Super4PCS::IndexedNormalSet <EigenPoint, Dim, 6, Scalar>;
    
    for(int i = 0; i < g_repeat; ++i)
    {
//        CALL_SUBTEST(( testFunction<Scalar,
//                                    EigenPoint,
//                                    NormalSet0>(eps, nbPoint) ));
        CALL_SUBTEST(( testFunction<Scalar,
                                    Dim,
                                    EigenPoint,
                                    NormalSet1>(eps, nbPoint) ));
    }
}

int main(int argc, const char **argv) {
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    using std::cout;
    using std::endl;

//    cout << "Extract pairs in 2 dimensions..." << endl;
//    callSubTests<float, 2>();
//    callSubTests<double, 2>();
//    callSubTests<long double, 2>();
//    cout << "Ok..." << endl;

    cout << "Extract pairs in 3 dimensions..." << endl;
    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
    cout << "Ok..." << endl;

//    cout << "Extract pairs in 4 dimensions..." << endl;
//    callSubTests<float, 4>();
//    callSubTests<double, 4>();
//    callSubTests<long double, 4>();
//    cout << "Ok..." << endl;

    return EXIT_SUCCESS;
}
