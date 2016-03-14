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


#ifndef _INDEXED_NORMAL_HEAL_SET_H_
#define _INDEXED_NORMAL_HEAL_SET_H_

#include "utils.h"
#include "chealpix.h"
#include <vector>

namespace Super4PCS{

/*!
  Work only in 3D, based on healpix.
  Scalar are constrained to be double
  Points are constrained to be Eigen vec3d
 */
class IndexedNormalHealSet{
public:
    typedef Eigen::Vector3d Point;
    typedef Eigen::Vector3i Index3D;
  typedef std::vector<std::vector<unsigned int>> ChealMap;

#ifdef DEBUG
#define VALIDATE_INDICES true
#pragma message "Index validation is enabled in debug mode"
#else
#define VALIDATE_INDICES false
#endif

  //! \brief State of the index validation, disabled when compiled in release mode
  enum{ INDEX_VALIDATION_ENABLED = VALIDATE_INDICES };

#undef VALIDATE_INDICES
  
private:
  double _epsilon;
  int _resolution;
  int _egSize;    //! <\brief Size of the euclidean grid for each dimension
  long _ngLength;  //! ,\brief Length of the normal map from healpix
  std::vector<ChealMap*> _grid;

  //! \brief Return the coordinates corresponding to index coord
  inline int indexCoordinates( const Index3D& coord) const
  { return Utils::UnrollIndexLoop<INDEX_VALIDATION_ENABLED>( coord,  2,  _egSize ); }

  //! \brief Return the index corresponding to position p
  inline int indexPos   ( const Point& p) const{
    return indexCoordinates( coordinatesPos(p) );
  }
  
  //! \brief Return the index corresponding to normal n  \warning Bounds are not tested
  inline int indexNormal( const Point& n) const {
    long id;
    vec2pix_ring(_resolution, n.data(), &id);
    return id;
  }

  //! \brief Return the coordinates corresponding to position p
  inline Index3D coordinatesPos   ( const Point& p) const
  {
    return (p/_epsilon)
            .unaryExpr(std::ptr_fun<Point::Scalar,Point::Scalar>(std::floor))
            .cast<typename Index3D::Scalar>();
  }
  
  
public:
  inline IndexedNormalHealSet(double epsilon, int resolution = 4)
  : _epsilon(epsilon), _resolution(resolution) {
    // We need to check if epsilon is a power of two and correct it if needed
    const int gridDepth = -std::log2(epsilon);
    _egSize = std::pow(2,gridDepth);
    _grid = std::vector<ChealMap*> (std::pow(_egSize, 3), NULL);
    
    _ngLength = nside2npix(resolution);
    
    _epsilon = 1.f/_egSize;
  }
  
  virtual ~IndexedNormalHealSet();
  
  //! \brief Add a new couple pos/normal, and its associated id
  bool addElement(const Point& pos, 
                  const Point& normal, 
                  unsigned int id);

  template <typename otherPoint>
  inline bool addElement(const otherPoint& pos,
                         const otherPoint& normal,
                         unsigned int id){
      return addElement(pos.template cast<double>(), normal.template cast<double>(), id);
  }
  
  //! \return NULL if the grid does not exist or p is out of bound
  inline ChealMap* getMap(const Point& p) { 
    const int pId = indexPos(p);
    if (pId == -1) return NULL;
    return _grid[pId]; 
  }

  //! \return a vector of maps containing points that can be close to p (according to input epsilon)
  inline std::vector<ChealMap*> getEpsilonMaps(const Point& p){
     const Index3D pId3 = coordinatesPos(p);
     const int pId = indexCoordinates(pId3);
     std::vector<ChealMap*> result;

     if (pId == -1) return result;


     // Here we extract the 8-neighorhood of the map in pId
     // We use three umbricated for loops, need to check how loop unrolling
     // is performed.
     // The following is really ugly, metaprog would be better
     // \FIXME
     const int lastId = _egSize-1;
     int imin = pId3[0] == 0 ? 0 : -1;
     int jmin = pId3[1] == 0 ? 0 : -1;
     int kmin = pId3[2] == 0 ? 0 : -1;
     int imax = pId3[0] == lastId ? 1 : 2;
     int jmax = pId3[1] == lastId ? 1 : 2;
     int kmax = pId3[2] == lastId ? 1 : 2;

     for (int i = imin; i < imax; ++i) // x axis
         for (int j = jmin; j < jmax; ++j) // y axis
             for (int k = kmin; k < kmax; ++k) // z axis
             {
                 const int id = indexCoordinates(pId3 + Index3D(i,j,k));
                 ChealMap* g = _grid[id];
                 if (g != nullptr)
                     result.push_back(g);
             }
     return result;
  }

  
  //! Get closest points in euclidean space
  void getNeighbors( const Point& p, 
                     std::vector<unsigned int>&nei);
  //! Get closest points in euclidean an normal space
  void getNeighbors( const Point& p, 
                     const Point& n,
                     std::vector<unsigned int>&nei);
  //! Get closest poitns in euclidean an normal space with angular deviation
  void getNeighbors( const Point& p,
                     const Point& n,
                     double alpha,
                     std::vector<unsigned int>&nei);

}; // class IndexedNormalHealSet
} // namespace Super4PCS

#endif // _INDEXED_NORMAL_SET_H_
