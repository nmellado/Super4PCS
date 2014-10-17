#ifndef BBOX_H
#define BBOX_H


#include <iostream>
#include <limits>
using namespace std;


template <class Point3D>
class BoundingBox3D
{
 public:

  BoundingBox3D () :
          _min(Point3D(  std::numeric_limits<float>::max()/2,
                      std::numeric_limits<float>::max()/2,
                      std::numeric_limits<float>::max()/2)),
          _max(Point3D(- std::numeric_limits<float>::max()/2,
                    - std::numeric_limits<float>::max()/2,
                    - std::numeric_limits<float>::max()/2)) { }
  BoundingBox3D ( const Point3D& min, const Point3D& max ) : _min(min), _max(max) { }
  BoundingBox3D ( const BoundingBox3D& bb ) : _min(bb._min), _max(bb._max) { }
  ~BoundingBox3D () { }

  inline BoundingBox3D& operator= ( const BoundingBox3D& bb ) {
    _min = bb._min;
    _max = bb._max;
    return (*this);
  }

  inline void extendTo ( const float x, const float y, const float z ) {
    if (x > _max.x) _max.x = x;
    if (x < _min.x) _min.x = x;
    if (y > _max.y) _max.y = y;
    if (y < _min.y) _min.y = y;
    if (z > _max.z) _max.z = z;
    if (z < _min.z) _min.z = z;
  }

  inline void extendTo ( const Point3D& vec ) {
    extendTo(vec.x,vec.y,vec.z);
  }

  inline void extendTo ( const BoundingBox3D& bb ) {
    extendTo(bb._min);
    extendTo(bb._max);
  }


  inline float width () {
    return (_max.x - _min.x);
  }

  inline float height () {
    return (_max.y - _min.y);
  }

  inline float depth () {
    return (_max.z - _min.z);
  }

  inline Point3D getCenter () const { return Point3D ((_min.x + _max.x)/2.,
                                                      (_min.y + _max.y)/2.,
                                                      (_min.z + _max.z)/2.); }

  inline Point3D& min () { return _min; }
  inline Point3D min () const { return _min; }

  inline Point3D& max () { return _max; }
  inline Point3D max () const { return _max; }

  inline void clear () {
      _min = Point3D(   std::numeric_limits<float>::max()/2,
                        std::numeric_limits<float>::max()/2,
                        std::numeric_limits<float>::max()/2);
      _max = Point3D( - std::numeric_limits<float>::max()/2,
                      - std::numeric_limits<float>::max()/2,
                      - std::numeric_limits<float>::max()/2);
  }
 private:
  Point3D _min;
  Point3D _max;
};

#endif // BBOX_H




