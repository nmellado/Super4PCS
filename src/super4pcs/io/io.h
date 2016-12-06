#ifndef IO_H
#define IO_H

#include "shared4pcs.h"
#include "utils/disablewarnings.h"

#include <fstream>
#include <iostream>
#include <string>
#include <opencv2/highgui/highgui.hpp>

#include "Eigen/Core"

#ifndef _MSC_VER
#include <sys/time.h>
#include <unistd.h>
#endif
#include <stdlib.h>


struct tripple {
  int a;
  int b;
  int c;
  int n1;
  int n2;
  int n3;
  int t1;
  int t2;
  int t3;
  tripple() {}
  tripple(int _a, int _b, int _c) : a(_a), b(_b), c(_c) {}
};


using namespace std;
using namespace match_4pcs;

class IOManager{
public:
  enum MATRIX_MODE {
      POLYWORKS //! <\brief Matrix file to be loaded and applied to polyworks layers
  };

public:
  // Obj read/write simple functions.
  bool ReadObject(const char *name, vector<Point3D> &v, vector<cv::Point2f> &tex_coords,
                  vector<cv::Point3f> &normals, vector<tripple> &tris,
                  vector<std::string> &mtls);
  bool WriteObject(const char *name, const vector<Point3D> &v,
                   const vector<cv::Point2f> &tex_coords, const vector<cv::Point3f> &normals,
                   const vector<tripple> &tris, const vector<string> &mtls);

  bool WriteMatrix(const string& name, const cv::Mat& mat, MATRIX_MODE mode);
private:
  bool
  ReadPly(const char *name, vector<Point3D> &v, vector<cv::Point3f> &normals);

  /*!
   * \brief ReadPtx
   * \param name
   * \param v
   * \return
   *
   * \note Transformations declared in file are ignored
   *
   * Implementation inspired by
   *            http://github.com/adasta/pcl_io_extra/blob/master/src/ptx_io.cpp
   */
  bool
  ReadPtx(const char *name, vector<Point3D> &v);
  
  bool 
  ReadObj(const char *name, vector<Point3D> &v, vector<cv::Point2f> &tex_coords,
          vector<cv::Point3f> &normals, vector<tripple> &tris,
          vector<std::string> &mtls);
                  
  bool
  WritePly(string name, const vector<Point3D> &v, const vector<cv::Point3f> &normals);

  bool 
  WriteObj(string name, const vector<Point3D> &v,
           const vector<cv::Point2f> &tex_coords, const vector<cv::Point3f> &normals,
           const vector<tripple> &tris, const vector<string> &mtls);


  /*!
   * \brief formatPolyworksMatrix Format 4x4 matrice so it can be loaded by polyworks
   * \param mat
   * \param sstr
   * \return
   */
  std::ofstream &
  formatPolyworksMatrix(const cv::Mat& mat, std::ofstream &sstr);
}; // class IOMananger

#endif // IO_H
