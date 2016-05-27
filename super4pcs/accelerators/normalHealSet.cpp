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
// Technically, we map the algorithm as an 'instance problem' and solve it 
// efficiently using a smart indexing data organization. The algorithm is 
// simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in
// significant speedup over alternative approaches and allows unstructured 
// efficient acquisition of scenes at scales previously not possible. Complete 
// source code and datasets are available for research use at 
// http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.


#include <math.h>
#include <set>
#include <Eigen/Geometry>
#include <iostream>

#include "normalHealSet.h"

using namespace std;

namespace Super4PCS{

IndexedNormalHealSet::~IndexedNormalHealSet(){
  for(unsigned int i = 0; i != _grid.size(); i++)
    delete(_grid[i]);
}

bool
IndexedNormalHealSet::addElement( 
  const Point& p,
  const Point& n,
  unsigned int id)
{
  const int pId = indexPos(p);
  if (pId == -1) return false;
  
  const int nId = indexNormal(n); 
  if (nId == -1) return false;
  
  if (_grid[pId] == NULL) _grid[pId] = new ChealMap(_ngLength);
  (_grid[pId])->at(nId).push_back(id);
  
  return true;
}


void 
IndexedNormalHealSet::getNeighbors( 
  const Point& p, 
  std::vector<unsigned int>&nei)
{
  using ChealMapIterator = ChealMap::const_iterator;
  ChealMap* grid = getMap(p);
  if ( grid == NULL ) return;
  
  for(ChealMapIterator it = grid->cbegin();
      it != grid->cend(); it++){
    const std::vector<unsigned int>& lnei = *it;
    nei.insert( nei.end(), lnei.begin(), lnei.end() );
  }
}


void 
IndexedNormalHealSet::getNeighbors( 
  const Point& p, 
  const Point& n,
  std::vector<unsigned int>&nei)
{
  ChealMap* grid = getMap(p);
  if ( grid == NULL ) return;
  
  const std::vector<unsigned int>& lnei = grid->at(indexNormal(n));
  nei.insert( nei.end(), lnei.begin(), lnei.end() );
}


void 
IndexedNormalHealSet::getNeighbors( 
  const Point& p, 
  const Point& n,
  double cosAlpha,
  std::vector<unsigned int>&nei)
{
  //ChealMap* grid = getMap(p);
  std::vector<ChealMap*> grids = getEpsilonMaps(p);
  if ( grids.empty() ) return;
  
  const double alpha          = std::acos(cosAlpha);
  //const double perimeter      = double(2) * M_PI * std::atan(alpha);
  const unsigned int nbSample = std::pow(2,_resolution+1);
  const double angleStep      = double(2) * M_PI / double(nbSample);

  
  const double sinAlpha       = std::sin(alpha);
  
  Eigen::Quaternion<double> q;
  q.setFromTwoVectors(Point(0.,0.,1.), n);
  
  // store a pair with
  // first  = grid id in grids
  // second = normal id in grids[first]
  typedef std::pair<unsigned int,unsigned int> PairId;
  std::set< PairId > colored;
  const int nbgrid = grids.size();
  
  // Do the rendering independently of the content
  for(unsigned int a = 0; a != nbSample; a++){
    double theta    = double(a) * angleStep;
    const Point dir = ( q * Point(sinAlpha*std::cos(theta), 
                              sinAlpha*std::sin(theta), 
                              cosAlpha ) ).normalized();
    int id = indexNormal( dir );      

    for (int i = 0; i != nbgrid; ++i){
        if(grids[i]->at(id).size() != 0){
          colored.insert(PairId(i,id));
        }
    }

  }
  
  for( std::set<PairId>::const_iterator it = colored.cbegin();
       it != colored.cend(); it++){
    const std::vector<unsigned int>& lnei = grids[it->first]->at(it->second);
    nei.insert( nei.end(), lnei.begin(), lnei.end() );
  }
}

} // namespace Super4PCS
