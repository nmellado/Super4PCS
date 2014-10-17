#include "4pcs.h"
#include "Eigen/Dense"

#include <fstream>
#include <iostream>
#include <string>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>

#include "io/io.h"

#define sqr(x) ((x) * (x))

using namespace std;
using namespace match_4pcs;

// First input.
std::string input1 = "input1.obj";

// Second input.
std::string input2 = "input2.obj";

// Output. The transformed second input.
std::string output = "output.obj";

// Delta (see the paper).
double delta = 5.0;

// Estimated overlap (see the paper).
double overlap = 0.2;

// Threshold of the computed overlap for termination. 1.0 means don't terminate
// before the end.
double thr = 1.0;

// Maximum norm of RGB values between corresponded points. 1e9 means don't use.
double max_color = 150;

// Number of sampled points in both files. The 4PCS allows a very aggressive
// sampling.
int n_points = 200;

// Maximum angle (degrees) between corresponded normals.
double norm_diff = 90.0;

// Maximum allowed computation time.
int max_time_seconds = 10;

bool use_super4pcs = true;

void getArgs(int argc, char **argv) {
  int i = 1;
  while (i < argc) {
    if (!strcmp(argv[i], "-i")) {
      input1 = std::string(argv[++i]);
      input2 = std::string(argv[++i]);
    } else if (!strcmp(argv[i], "-o")) {
      overlap = atof(argv[++i]);
    } else if (!strcmp(argv[i], "-d")) {
      delta = atof(argv[++i]);
    } else if (!strcmp(argv[i], "-c")) {
      max_color = atof(argv[++i]);
    } else if (!strcmp(argv[i], "-t")) {
      max_time_seconds = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-a")) {
      norm_diff = atof(argv[++i]);
    } else if (!strcmp(argv[i], "-n")) {
      n_points = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-r")) {
      output = argv[++i];
    } else if (!strcmp(argv[i], "-x")) {
      use_super4pcs = false;
    } else if (!strcmp(argv[i], "-h")) {
      fprintf(stderr, "\nUsage: %s -i input1 input2\n\n", argv[0]);
      fprintf(stderr, "\t[ -o overlap (%2.2f) ]\n", overlap);
      fprintf(stderr, "\t[ -d delta (%2.2f) ]\n", delta);
      fprintf(stderr, "\t[ -n n_points (%d) ]\n", n_points);
      fprintf(stderr, "\t[ -a norm_diff (%f) ]\n", norm_diff);
      fprintf(stderr, "\t[ -c max_color_diff (%f) ]\n", max_color);
      fprintf(stderr, "\t[ -t max_time_seconds (%d) ]\n", max_time_seconds);
      fprintf(stderr, "\t[ -r result_file_name (%s) ]\n", output.c_str());
      fprintf(stderr, "\t[ -x (use 4pcs: false by default) ]\n");
      exit(0);
    } else if (argv[i][0] == '-') {
      cerr << "Unknown flag\n";
      exit(-1);
    };
    i++;
  }
}


void CleanInvalidNormals( vector<Point3D> &v, 
                          vector<cv::Point3f> &normals){
  if (v.size() == normals.size()){
    vector<Point3D>::iterator itV = v.begin();
    vector<cv::Point3f>::iterator itN = normals.begin();
  
    float norm;
    unsigned int nb = 0;
    for( ; itV != v.end(); ){
      norm = cv::norm((*itV).normal());
      if (norm < 0.1){
        itN = normals.erase(itN);
        itV = v.erase(itV);
        nb++;
      }else{
        if (norm != 1.){
          (*itN).x /= norm;
          (*itN).y /= norm;
          (*itN).z /= norm;
        }
        itV++;
        itN++;
      }
    }
    
    if (nb != 0){
      cout << "Removed " << nb << " invalid points/normals" << endl; 
    }
  }
}


int main(int argc, char **argv) {

  vector<Point3D> set1, set2;
  vector<cv::Point2f> tex_coords1, tex_coords2;
  vector<cv::Point3f> normals1, normals2;
  vector<tripple> tris1, tris2;
  vector<std::string> mtls1, mtls2;

  getArgs(argc, argv);
  
  IOManager iomananger;

  // Read the inputs.
  if (!iomananger.ReadObject((char *)input1.c_str(), set1, tex_coords1, normals1, tris1,
                  mtls1)) {
    perror("Can't read input set1");
    exit(-1);
  }

  if (!iomananger.ReadObject((char *)input2.c_str(), set2, tex_coords2, normals2, tris2,
                  mtls2)) {
    perror("Can't read input set2");
    exit(-1);
  }
  
  // clean only when we have pset to avoid wrong face to point indexation
  if (tris1.size() == 0)
    CleanInvalidNormals(set1, normals1);
  if (tris2.size() == 0)
    CleanInvalidNormals(set2, normals2);

  // Our matcher.
  Match4PCSOptions options;  

  // Set parameters.
  cv::Mat mat = cv::Mat::eye(4, 4, CV_64F);
  options.overlap_estimation = overlap;
  options.sample_size = n_points;
  options.max_normal_difference = norm_diff;
  options.max_color_distance = max_color;
  options.max_time_seconds = max_time_seconds;
  options.delta = delta;
  // Match and return the score (estimated overlap or the LCP).  
  float score = 0;
  
  if(use_super4pcs){
    MatchSuper4PCS matcher(options);
    cout << "Use Super4PCS" << endl;
    score = matcher.ComputeTransformation(set1, &set2, &mat);
  }else{
    Match4PCS matcher(options);
    cout << "Use old 4PCS" << endl;
    score = matcher.ComputeTransformation(set1, &set2, &mat);  
  }

  cout << "Score: " << score << endl;
  cerr <<  score << endl;
  printf("(Homogeneous) Transformation from %s to %s:\n", input2.c_str(),
         input1.c_str());
  printf(
      "\n\n%25.3f %25.3f %25.3f %25.3f\n%25.3f %25.3f %25.3f %25.3f\n%25.3f "
      "%25.3f %25.3f %25.3f\n%25.3f %25.3f %25.3f %25.3f\n\n",
      mat.at<double>(0, 0), mat.at<double>(0, 1), mat.at<double>(0, 2),
      mat.at<double>(0, 3), mat.at<double>(1, 0), mat.at<double>(1, 1),
      mat.at<double>(1, 2), mat.at<double>(1, 3), mat.at<double>(2, 0),
      mat.at<double>(2, 1), mat.at<double>(2, 2), mat.at<double>(2, 3),
      mat.at<double>(3, 0), mat.at<double>(3, 1), mat.at<double>(3, 2),
      mat.at<double>(3, 3));

  // If the default images are the input then we need to test the result.
  if (input1 == "input1.obj" && input2 == "input2.obj") {  // test!
    float gt_mat[16] = {0.977,  -0.180,  -0.114, 91.641, 0.070, 0.778,
                        -0.624, 410.029, 0.201,  0.602,  0.773, 110.810,
                        0.000,  0.000,   0.000,  1.000};
    float norm_val = 0.0;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        norm_val += sqr(mat.at<double>(i, j) - gt_mat[i * 4 + j]);
      }
    }
    assert(norm_val > 9.46881e-10);
  }
  
  iomananger.WriteObject((char *)output.c_str(), 
                         set2, 
                         tex_coords2, 
                         normals2, 
                         tris2,
	                       mtls2);

  return 0;
}
