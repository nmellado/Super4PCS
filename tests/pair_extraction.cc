#include "4pcs.h"

#include "Eigen/Dense"

#include <fstream>
#include <iostream>
#include <string>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>
#include "intersection.h"
#include "primitives.h"
#include "bbox.h"

#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>

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

struct PairCreationFunctor{
  typedef pair<unsigned int, unsigned int>ResPair;
  vector< ResPair >pairs;
  
  std::vector<unsigned int> ids;
  
  inline void beginPrimitiveCollect(int primId){
  }
  
  inline void endPrimitiveCollect(int primId){
  }
  inline void process(int primId, int pointId){
    if(pointId >= 10)
      pairs.push_back(ResPair(pointId, primId));
  }
};


#define dim 3
typedef double Scalar;
typedef Eigen::Matrix<Scalar, dim, 1> EigenPoint;
typedef HyperSphere< EigenPoint, dim, Scalar > Sphere ;


void
testIntersection(std::vector<EigenPoint>&points, 
                 std::vector<Sphere> &primitives,
                 Scalar epsilon){
                 
  unsigned int nbPoints = points.size();
    
  cout << "Compute for " << nbPoints << " points ...." << flush;
  IntersectionFunctor<Sphere, EigenPoint, dim, Scalar> IF;
  
  Utilities::Timer t; 
  PairCreationFunctor functor;
  
  Utilities::startTimer(t);
  //std::vector< std::pair<unsigned int, unsigned int> > pairs = 
  IF.process(primitives, points, epsilon, 0, functor);
  Utilities::stopTimer(t);
  
  
  cout << "DONE " << t << endl;

  cout << functor.pairs.size() << " pairs extracted: " 
       << " for epsilon = " << epsilon << endl;
       
  Utilities::startTimer(t);
       
  vector< pair<unsigned int, unsigned int> > p2;
  p2.reserve(nbPoints*primitives.size());
  
  cout << "Compute euclidean...." << flush;
  for(unsigned int i = 0; i != nbPoints; i++)
    for(unsigned int j = 0; j != primitives.size(); j++)
       if (primitives[j].intersectPoint(points[i], epsilon))
        p2.push_back(std::pair<unsigned int, unsigned int>(i,j));
  Utilities::stopTimer(t);
        
  cout << "DONE " << t << endl;  
  cout << p2.size() << " pairs extracted using complete search" << endl;
  
#ifdef DEBUG_PLOT_SEGMENT
// The gnuplot command is:
// set xrange [0:1]
// set yrange [0:1]
// plot "./detected.plot" u 1:2 w points t "detected", "./euclidean.plot" u 1:2 w points t "euclidean"

#if dim != 2
#error Gnuplot files generation is valid only in 2D
#endif

  ofstream detected, outliers, euclidean;
  
  detected.open ("detected.plot",  ios::out | ios::trunc);
  euclidean.open("euclidean.plot", ios::out | ios::trunc);  
  
  if (detected.is_open() && euclidean.is_open()) {

    // print all points selected using euclidean distance
//    for(std::vector< std::pair<unsigned int, unsigned int> >::const_iterator it =
//        p2.begin(); it != p2.end(); it++)
//      euclidean << points[(*it).first].transpose() << endl;
      
      
    for(unsigned int i = 0; i != nbPoints; i++){      
      for(unsigned int j = 0; j != primitives.size(); j++){

        if (primitives[j].intersectPoint(points[i], epsilon)){
          euclidean << points[i].transpose() << endl;
          break;
        }
      } 
    }
      
    // print all points selected using the structure
    for(std::vector< std::pair<unsigned int, unsigned int> >::const_iterator it =
        functor.pairs.begin(); it != functor.pairs.end(); it++)
      detected << points[(*it).first].transpose() << endl;
     
    
    detected.close ();
    euclidean.close();
  }
#endif
}



#if dim == 3
void
testWithInput(const vector<Point3D> &set, float epsilon, int nbPrim){
  std::vector<EigenPoint> points;
  std::vector<Sphere> primitives; 
  
  points.reserve(set.size());
  primitives.reserve(set.size());
  
  BoundingBox3D<Point3D> bbox;
  
  // init Eigen data and primitives
  for (int i = 0; i < set.size(); ++i) {
    points.push_back( EigenPoint( set[i].x, set[i].y, set[i].z ));
    bbox.extendTo(set[i]);
  }    
  
  Point3D gcenter3D = bbox.getCenter();
  EigenPoint gcenter (gcenter3D.x, gcenter3D.y, gcenter3D.z);
  float ratio = MAX(bbox.depth(), MAX(bbox.width(), bbox.height()));
  EigenPoint half = EigenPoint::Ones() * 0.5f;
  
  bbox.clear();
  epsilon /= ratio;
  
  for (int i = 0; i < points.size(); ++i) {
    points[i] -= gcenter;
    points[i] /= ratio;
    points[i] += half;
    
    if (nbPrim != 0){
      if ( i < nbPrim)
        primitives.push_back(Sphere(points[i], 0.0078125));
    }else
      primitives.push_back(Sphere(points[i], 0.0078125));
  }
  
  testIntersection(points, primitives, epsilon);
}
#endif

void
testWithRandom(){
    // Init input point
  unsigned int nbPoints = 10000;
  
  std::vector<EigenPoint> points;
  std::vector<Sphere> primitives; 
  
  EigenPoint half (EigenPoint::Ones()/2.f);
  
  cout << "Input points: " << endl;
  for(unsigned int i = 0; i != nbPoints; i++){
    EigenPoint p (0.5f*EigenPoint::Random() + half);
    points.push_back(p);   
    
    if(i < 30)
    primitives.push_back(Sphere( p, 0.2));
  } 
    
    
  // Compute intersection
  float epsilon = 0.125f/4.f;//0.125f/2.f;
  
  testIntersection(points, primitives, epsilon);
}


//void
//testNormalIndex() {
//    // Init input point
//  unsigned int nbPoints = 100000;

//  Scalar epsilon = 0.125f/2.f;
//  IndexedNormalSet<EigenPoint, dim, 7, Scalar> nset(epsilon);
//  IndexedNormalHealSet nhset(epsilon, 4);
//  
//  std::vector<EigenPoint> points;
//  std::vector<EigenPoint> normals;
//  EigenPoint normal(1,0,0);
//  
//  Utilities::Timer t;
//  
//  EigenPoint half (EigenPoint::Ones()/2.f);
//  
//  float tEuclidean = 0,
//        tAngular   = 0;
//  // filling
//  for(unsigned int i = 0; i != nbPoints; i++){
//    EigenPoint p (0.5f*EigenPoint::Random() + half);
//    EigenPoint normal (EigenPoint::Random());
//    normal.normalize();
//    points.push_back(p);  
//    normals.push_back(normal); 
//      
//    Utilities::startTimer(t);
//    nset.addElement(p, normal, i);
//    Utilities::stopTimer(t);
//    tEuclidean += t;
//    
//    Utilities::startTimer(t);
//    nhset.addElement(p, normal, i);
//    Utilities::stopTimer(t);
//    tAngular += t;
//  }
//  
//  cout << "Fill Euclidean: " << tEuclidean << endl;
//  cout << "Fill Anguular:  " << tAngular << endl;
//  
//  // Query
//  unsigned int nbTry = 100000;
//  
//  tEuclidean = 0;
//  tAngular   = 0;
//  
//  Scalar eEuclidean = 0,
//         eAngular   = 0;
//  Scalar dEuclidean = 0,
//         dAngular   = 0;
//  unsigned int nbNeiEuclidean = 0;
//  unsigned int nbNeiAngular   = 0;
//  for(unsigned int i = 0; i != nbTry; i++){
//  std::vector<unsigned int> nei;    
//  EigenPoint query  = (0.5f*EigenPoint::Random() + half);
//  EigenPoint queryn = normal;  
//  
////  cout << "epsilon = " << epsilon << endl;
//  
////  cout << "Start gathering neighborhood" << endl;
//  Utilities::startTimer(t);
//  nset.getNeighbors( query, 
//                     queryn,
//                     0.2,
//                     nei);
//  Utilities::stopTimer(t);
//  tEuclidean += t;
//  
//  for (unsigned int k = 0; k != nei.size(); k++){   
//    int id = nei[k]; 
//    dEuclidean += (points[id] - query).norm();
//    eEuclidean += std::abs(normals[id].dot(queryn));
//  }   
//  nbNeiEuclidean += nei.size();
//  
//  nei.clear();
//  Utilities::startTimer(t);
//  nhset.getNeighbors( query, 
//                     queryn,
//                     0.2,
//                     nei);
//  Utilities::stopTimer(t);
//  tAngular += t;
//  
//  for (unsigned int k = 0; k != nei.size(); k++){   
//    int id = nei[k]; 
//    dAngular += (points[id] - query).norm();
//    eAngular += std::abs(normals[id].dot(queryn));
//  }  
//  nbNeiAngular += nei.size();
//  
//  
////  cout << "Done" << endl;
////              
////  for (unsigned int k = 0; k != nei.size(); k++){   
////    int id = nei[k]; 
////    cout << (points[id] - query).norm() << " " << normals[id].dot(queryn) << endl;
////  }   
//  }
//  
//  cout << "Query Euclidean: " << tEuclidean << endl;
//  cout << "Query Angular:  " << tAngular << endl;
//  
//  cout << "Error Euclidean distances: " << dEuclidean / Scalar(nbNeiEuclidean) << endl;
//  cout << "Error Angular distances:   " << dAngular / Scalar(nbNeiAngular) << endl;
//  
//  cout << "Error Euclidean dot: " << eEuclidean / Scalar(nbNeiEuclidean) << endl;
//  cout << "Error Angular dot:   " << eAngular / Scalar(nbNeiAngular) << endl;
//  
//  cout << "Nb Points Euclidean: " << nbNeiEuclidean << endl;
//  cout << "Nb Points Angular:   " << nbNeiAngular << endl;
//  //for (unsigned int i = 0; i != nbPoints; i++){
//  //  if((points[i] - query).norm() <= epsilon)
//  //    cout << (points[i] - query).norm() << endl;  
//  //}
//}


int main(int argc, char **argv) {
  return true;
}
