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

namespace Super4PCS{

/*!
  Work only in 3D, based on healpix.
  Scalar are constrained to be double
  Points are constrained to be Eigen vec3d
 */
class IndexedNormalHealSet{
public:
  typedef Eigen::Vector3d Point;
  typedef std::vector<std::vector<unsigned int>> ChealMap;
  
  
private:
  double _epsilon;
  int _resolution;
  int _egSize;    //! <\brief Size of the euclidean grid for each dimension
  long _ngLength;  //! ,\brief Length of the normal map from healpix
  std::vector<ChealMap*> _grid;
  
  // Get the index corresponding to position p \warning Bounds are not tested
  inline int indexPos   ( const Point& p) const{
    return UnrollIndexLoop( coordinatesPos(p),  2/*dim-1*/,  _egSize );  
  }
  
  // Get the index corresponding to normal n   \warning Bounds are not tested
  inline int indexNormal( const Point& n) const {
    long id;
    vec2pix_ring(_resolution, n.data(), &id);
    return id;
  }
  
  // Get the coordinates corresponding to position p \warning Bounds are not tested
  inline Point coordinatesPos   ( const Point& p) const
  { return p/_epsilon;  }
  // Get the index corresponding to normal n   \warning Bounds are not tested
  
  // Get the coordinates corresponding to position p \warning Bounds are not tested
  inline int indexCoordinatesPos   ( const Point& pCoord) const{  
    return UnrollIndexLoop( pCoord,  2/*dim-1*/,  _egSize );  
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
  
  //! \return NULL if the grid does not exist or p is out of bound
  inline ChealMap* getMap(const Point& p) { 
    const int pId = indexPos(p);
    if (pId == -1) return NULL;
    return _grid[pId]; 
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
                     std::vector<unsigned int>&nei,
                     bool tryReverse = false);

}; // class IndexedNormalHealSet
} // namespace Super4PCS

#endif // _INDEXED_NORMAL_SET_H_
