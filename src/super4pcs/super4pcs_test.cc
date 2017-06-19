#include "algorithms/4pcs.h"
#include "algorithms/super4pcs.h"
#include "Eigen/Dense"

#include <fstream>
#include <iostream>
#include <string>

#include "io/io.h"
#include "utils/geometry.h"

#define sqr(x) ((x) * (x))

using namespace std;
using namespace match_4pcs;

// First input.
std::string input1 = "input1.obj";

// Second input.
std::string input2 = "input2.obj";

// Output. The transformed second input.
std::string output = "";
// Default name for the '.obj' output file
std::string defaultObjOutput = "output.obj";
// Default name for the '.ply' output file
std::string defaultPlyOutput = "output.ply";

// Transformation matrice.
std::string outputMat = "";

// Sampled cloud 1
std::string outputSampled1 = "";

// Sampled cloud 2
std::string outputSampled2 = "";

// Delta (see the paper).
double delta = 5.0;

// Estimated overlap (see the paper).
double overlap = 0.2;

// Threshold of the computed overlap for termination. 1.0 means don't terminate
// before the end.
double thr = 1.0;

// Maximum norm of RGB values between corresponded points. 1e9 means don't use.
double max_color = -1;

// Number of sampled points in both files. The 4PCS allows a very aggressive
// sampling.
int n_points = 200;

// Maximum angle (degrees) between corresponded normals.
double norm_diff = -1;

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
    } else if (!strcmp(argv[i], "-m")) {
      outputMat = argv[++i];
    } else if (!strcmp(argv[i], "-x")) {
      use_super4pcs = false;
    } else if (!strcmp(argv[i], "--sampled1")) {
      outputSampled1 = argv[++i];
    } else if (!strcmp(argv[i], "--sampled2")) {
      outputSampled2 = argv[++i];
    } else if (!strcmp(argv[i], "-h")) {
      fprintf(stderr, "\nUsage: %s -i input1 input2\n\n", argv[0]);
      fprintf(stderr, "\t[ -o overlap (%2.2f) ]\n", overlap);
      fprintf(stderr, "\t[ -d delta (%2.2f) ]\n", delta);
      fprintf(stderr, "\t[ -n n_points (%d) ]\n", n_points);
      fprintf(stderr, "\t[ -a norm_diff (%f) ]\n", norm_diff);
      fprintf(stderr, "\t[ -c max_color_diff (%f) ]\n", max_color);
      fprintf(stderr, "\t[ -t max_time_seconds (%d) ]\n", max_time_seconds);
      fprintf(stderr, "\t[ -r result_file_name (%s) ]\n", output.c_str());
      fprintf(stderr, "\t[ -m output matrix file (%s) ]\n", outputMat.c_str());
      fprintf(stderr, "\t[ -x (use 4pcs: false by default) ]\n");
      fprintf(stderr, "\t[ --sampled1 (output sampled cloud 1 -- debug+super4pcs only) ]\n");
      fprintf(stderr, "\t[ --sampled2 (output sampled cloud 2 -- debug+super4pcs only) ]\n");
      exit(0);
    } else if (argv[i][0] == '-') {
      cerr << "Unknown flag\n";
      exit(-1);
    };
    i++;
  }

  // if no output file (geometry/matrix) is set, force 3d mesh
  if (output.empty() && outputMat.empty()) output = defaultObjOutput;

}


int main(int argc, char **argv) {
  using namespace Super4PCS;

  vector<Point3D> set1, set2;
  vector<Eigen::Matrix2f> tex_coords1, tex_coords2;
  vector<typename Point3D::VectorType> normals1, normals2;
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
    Utils::CleanInvalidNormals(set1, normals1);
  if (tris2.size() == 0)
    Utils::CleanInvalidNormals(set2, normals2);

  // Our matcher.
  Match4PCSOptions options;

  // Set parameters.
  Match4PCSBase::MatrixType mat;
  options.overlap_estimation = overlap;
  options.sample_size = n_points;
  options.max_normal_difference = norm_diff;
  options.max_color_distance = max_color;
  options.max_time_seconds = max_time_seconds;
  options.delta = delta;
  // Match and return the score (estimated overlap or the LCP).
  typename Point3D::Scalar score = 0;

  try {

      if (use_super4pcs) {
          MatchSuper4PCS matcher(options);
          cout << "Use Super4PCS" << endl;
          score = matcher.ComputeTransformation(set1, &set2, mat);

          if(! outputSampled1.empty() ){
              std::cout << "Exporting Sampled cloud 1 to "
                        << outputSampled1.c_str()
                        << "..." << std::flush;
              iomananger.WriteObject((char *)outputSampled1.c_str(),
                                     matcher.getFirstSampled(),
                                     vector<Eigen::Matrix2f>(),
                                     vector<typename Point3D::VectorType>(),
                                     vector<tripple>(),
                                     vector<string>());
              std::cout << "DONE" << std::endl;
          }
          if(! outputSampled2.empty() ){
              std::cout << "Exporting Sampled cloud 2 to "
                        << outputSampled2.c_str()
                        << "..." << std::flush;
              iomananger.WriteObject((char *)outputSampled2.c_str(),
                                     matcher.getSecondSampled(),
                                     vector<Eigen::Matrix2f>(),
                                     vector<typename Point3D::VectorType>(),
                                     vector<tripple>(),
                                     vector<string>());
              std::cout << "DONE" << std::endl;
          }
      }
      else {
          Match4PCS matcher(options);
          cout << "Use old 4PCS" << endl;
          score = matcher.ComputeTransformation(set1, &set2, mat);
      }

  }
  catch (const std::exception& e) {
      std::cout << "[Error]: " << e.what() << '\n';
      std::cout << "Aborting with code -2 ..." << std::endl;
      return -2;
  }
  catch (...) {
      std::cout << "[Unknown Error]: Aborting with code -3 ..." << std::endl;
      return -3;
  }

  cout << "Score: " << score << endl;
  cout <<"(Homogeneous) Transformation from " << input2.c_str()
       << " to "<< input1.c_str() << ":\n";

  cout << mat << std::endl;


  if(! outputMat.empty() ){
      std::cout << "Exporting Matrix to "
                << outputMat.c_str()
                << "..." << std::flush;
      iomananger.WriteMatrix(outputMat, mat.cast<double>(), IOManager::POLYWORKS);
      std::cout << "DONE" << std::endl;
  }

  if (! output.empty() ){

      std::cout << "Exporting Registered geometry to "
                << output.c_str()
                << "..." << std::flush;
      iomananger.WriteObject((char *)output.c_str(),
                             set2,
                             tex_coords2,
                             normals2,
                             tris2,
                             mtls2);
      std::cout << "DONE" << std::endl;
  }

  return 0;
}
