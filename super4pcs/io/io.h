#ifndef IO_H
#define IO_H

#include "4pcs.h"

#include <fstream>
#include <iostream>
#include <string>
#include <opencv2/highgui/highgui.hpp>

#include "Eigen/Core"

#include <sys/time.h>
#include <unistd.h>
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
  // Obj read/write simple functions.
  bool ReadObject(char *name, vector<Point3D> &v, vector<cv::Point2f> &tex_coords,
                  vector<cv::Point3f> &normals, vector<tripple> &tris,
                  vector<std::string> &mtls);
  bool WriteObject(char *name, vector<Point3D> &v,
                   vector<cv::Point2f> &tex_coords, vector<cv::Point3f> &normals,
                   vector<tripple> &tris, vector<std::string> &mtls);
private:
  bool 
  ReadPly(char *name, vector<Point3D> &v, vector<cv::Point3f> &normals);
  
  bool 
  ReadObj(char *name, vector<Point3D> &v, vector<cv::Point2f> &tex_coords,
          vector<cv::Point3f> &normals, vector<tripple> &tris,
          vector<std::string> &mtls);
                  
  bool
  WritePly(string name, vector<Point3D> &v, vector<cv::Point3f> &normals);

  bool 
  WriteObj(string name, vector<Point3D> &v,
           vector<cv::Point2f> &tex_coords, vector<cv::Point3f> &normals,
           vector<tripple> &tris, vector<std::string> &mtls);
}; // class IOMananger

#endif // IO_H
