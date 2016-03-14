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
// This test runs super4PCS on multiple models and check the computed transformation
// matrix is the same than the one computed during a previous run.
// Dataset used: Armadillo scans, Stanford University Computer Graphics Laboratory
//               http://graphics.stanford.edu/data/3Dscanrep
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
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>

#include <boost/filesystem.hpp>

#include "io/io.h"

#include "testing.h"

#define sqr(x) ((x) * (x))

using namespace std;
using namespace match_4pcs;

typedef double Scalar;
enum {Dim = 3};
typedef Eigen::Transform<Scalar, Dim, Eigen::Affine> Transform;

const int nbSet = 2;

/*!
  Datasets and associated parameters
  In the current state we use a conservative configuration, better timing
  could be achieved by tuning a bit more the parameters.
  */

std::array<std::string, nbSet> confFiles = {
    "./datasets/bunny/data/bun.conf",
    "./datasets/armadillo/ArmadilloSide.conf"
};

std::array<Scalar, nbSet> deltas  = {
    0.001,
    0.001,
};

std::array<Scalar, nbSet> overlaps = {
    0.5,
    0.7,
};

std::array<Scalar, nbSet> n_points = {
    900,
    700,
};

// Maximum norm of RGB values between corresponded points. 1e9 means don't use.
double max_color = 1e9;

// Maximum angle (degrees) between corresponded normals.
double norm_diff = 360.0;

// Maximum allowed computation time.
int max_time_seconds = 1e9;

bool use_super4pcs = true;

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

/*!
  Read a configuration file from Standford 3D shape repository and
  output a set of filename and eigen transformations
  */
inline void
extractFilesAndTrFromStandfordConfFile(
        const std::string &confFilePath,
        std::vector<Transform>& transforms,
        std::vector<string>& files
        ){
    using namespace boost;
    using namespace std;

    VERIFY (filesystem::exists(confFilePath) && filesystem::is_regular_file(confFilePath));

    // extract the working directory for the configuration path
    const string workingDir = filesystem::path(confFilePath).parent_path().native();
    VERIFY (filesystem::exists(workingDir));

    // read the configuration file and call the matching process
    string line;
    ifstream confFile;
    confFile.open(confFilePath);
    VERIFY (confFile.is_open());

    while ( getline (confFile,line) )
    {
        istringstream iss (line);
        vector<string> tokens{istream_iterator<string>{iss},
                              istream_iterator<string>{}};

        // here we know that the tokens are:
        // [0]: keyword, must be bmesh
        // [1]: 3D object filename
        // [2-4]: target translation with previous object
        // [5-8]: target quaternion with previous object

        if (tokens.size() == 9){
            if (tokens[0].compare("bmesh") == 0){

                string inputfile = filesystem::path(confFilePath).parent_path().string()+string("/")+tokens[1];
                VERIFY(filesystem::exists(inputfile) && filesystem::is_regular_file(inputfile));

                // build the Eigen rotation matrix from the rotation and translation stored in the files
                Eigen::Matrix<Scalar, Dim, 1> tr (
                            std::atof(tokens[2].c_str()),
                            std::atof(tokens[3].c_str()),
                            std::atof(tokens[4].c_str()));

                Eigen::Quaternion<Scalar> quat(
                            std::atof(tokens[8].c_str()), // eigen starts by w
                            std::atof(tokens[5].c_str()),
                            std::atof(tokens[6].c_str()),
                            std::atof(tokens[7].c_str()));

                quat.normalize();

                Transform transform (Transform::Identity());
                transform.rotate(quat);
                transform.translate(-tr);

                transforms.push_back(transform);
                files.push_back(inputfile);


            }
        }
    }
    confFile.close();
}

void test_model(const vector<Transform> &transforms,
                const vector<string> &files,
                vector<Point3D> &mergedset,
                int i,
                int param_i){
    const string input1 = files.at(i-1);
    const string input2 = files.at(i);

    cout << "Matching " << input2.c_str() << endl;

    vector<Point3D> set1, set2;
    vector<cv::Point2f> tex_coords1, tex_coords2;
    vector<cv::Point3f> normals1, normals2;
    vector<tripple> tris1, tris2;
    vector<std::string> mtls1, mtls2;

    IOManager iomananger;
    VERIFY(iomananger.ReadObject((char *)input1.c_str(), set1, tex_coords1, normals1, tris1, mtls1));
    VERIFY(iomananger.ReadObject((char *)input2.c_str(), set2, tex_coords2, normals2, tris2, mtls2));

    // clean only when we have pset to avoid wrong face to point indexation
    if (tris1.size() == 0)
        CleanInvalidNormals(set1, normals1);
    if (tris2.size() == 0)
        CleanInvalidNormals(set2, normals2);

    // first transform the first mesh to its gt coordinates:
    // we compare only pairwise matching, so we don't want to
    // accumulate error during the matching process
    // Transforms Q by the new transformation.
    {
        cv::Mat transformation = cv::Mat::eye(4, 4, CV_64F);
        cv::eigen2cv(transforms[i-1].inverse().matrix(), transformation);
        for (int i = 0; i < set1.size(); ++i) {
            cv::Mat first(4, 1, CV_64F), transformed;
            first.at<double>(0, 0) = set1[i].x;
            first.at<double>(1, 0) = set1[i].y;
            first.at<double>(2, 0) = set1[i].z;
            first.at<double>(3, 0) = 1;
            transformed = transformation * first;
            set1[i].x = transformed.at<double>(0, 0);
            set1[i].y = transformed.at<double>(1, 0);
            set1[i].z = transformed.at<double>(2, 0);
        }
    }

    mergedset.insert(mergedset.end(), set1.begin(), set1.end());

    // Our matcher.
    Match4PCSOptions options;

    // Set parameters.
    cv::Mat mat = cv::Mat::eye(4, 4, CV_64F);
    options.overlap_estimation = overlaps[param_i];
    options.sample_size = n_points[param_i];
    options.max_normal_difference = norm_diff;
    options.max_color_distance = max_color;
    options.max_time_seconds = max_time_seconds;
    options.delta = deltas[param_i];

    Scalar score = 0.;

    if(use_super4pcs){
        MatchSuper4PCS matcher(options);
        cout << "./Super4PCS -i "
             << input1.c_str() << " "
             << input2.c_str()
             << " -o " << options.overlap_estimation
             << " -d " << options.delta
             << " -n " << options.sample_size
             << " -a " << options.max_normal_difference
             << " -c " << options.max_color_distance
             << " -t " << options.max_time_seconds
             << endl;
        score = matcher.ComputeTransformation(mergedset, &set2, &mat);
    }else{
        Match4PCS matcher(options);
        cout << "./Super4PCS -i "
             << input1.c_str() << " "
             << input2.c_str()
             << " -o " << options.overlap_estimation
             << " -d " << options.delta
             << " -n " << options.sample_size
             << " -a " << options.max_normal_difference
             << " -c " << options.max_color_distance
             << " -t " << options.max_time_seconds
             << " -x "
             << endl;
        score = matcher.ComputeTransformation(mergedset, &set2, &mat);
    }


#ifdef WRITE_OUTPUT_FILES
    stringstream iss;
    iss << input2;
    iss << "_aligned.ply";
    iomananger.WriteObject(iss.str().c_str(),
                           set2,
                           tex_coords2,
                           normals2,
                           tris2,
                           mtls2);
#endif

    // convert matrix to eigen matrix, and
    Eigen::Matrix<Scalar,4, 4> mat_eigen;
    cv::cv2eigen(mat,mat_eigen);

    Transform transformEst (mat_eigen);

    cout << "Reference: " << endl << transforms[i].matrix() << endl;
    cout << "Estimation: " << endl << transformEst.matrix() << endl;

    Eigen::Quaternion<Scalar>
            q    (transformEst.rotation()),
            qref (transforms[i].rotation());
    cout << " " << q.x()
         << " " << q.y()
         << " " << q.z()
         << " " << q.w()  << endl;
    cout << " " << qref.x()
         << " " << qref.y()
         << " " << qref.z()
         << " " << qref.w()  << endl;

    Scalar rotDiff = std::abs( ( q.vec().array().abs() - ( qref.vec().array().abs() ) ).abs().sum()) + std::abs(std::abs(q.w()) - std::abs(qref.w()));
    Scalar trDiff  = std::abs(transformEst.translation().dot(transforms[i].translation()));

    cout << "rotDiff = " << rotDiff << endl;
    cout << "trDiff = " << trDiff << endl;

    // use different tests, just to get more verbose output
    VERIFY(rotDiff <= 0.2);
    VERIFY(trDiff <= 0.1);
    VERIFY(rotDiff + trDiff <= 0.2);

#ifdef WRITE_OUTPUT_FILES
    stringstream iss2;
    iss2 << input1;
    iss2 << "_merged.ply";
    iomananger.WriteObject(iss2.str().c_str(),
                           mergedset,
                           vector<cv::Point2f>(),
                           vector<cv::Point3f>(),
                           vector<tripple>(),
                           vector<std::string>());
#endif
}

int main(int argc, const char **argv) {
    using std::string;

    const char* custom_argv [1] = {"matching"};

    if(!init_testing(1, custom_argv))
    {
        return EXIT_FAILURE;
    }

    if(argc < 2) // we accept only 1 argument: the test id
    {
        return EXIT_FAILURE;
    }


    int i = std::atoi(argv[1]);
    if (i > confFiles.size()) return EXIT_FAILURE;

    cout << "Starting to process test " << i << endl;

    vector<Transform> transforms;
    vector<string> files;
    extractFilesAndTrFromStandfordConfFile(confFiles.at(i), transforms, files);

    VERIFY(transforms.size() == files.size());
    const int nbTests = transforms.size()-1;


    // In this test we match each frame to the union of the previously matched
    // frames.
    // In practice we can't use the output of Super4PCS directly, it would require
    // a local ICP to avoid error accumulation. So we sum-up the GT transformations,
    // and check our Super4PCS is working well by comparing the estimated transformation
    // matrix and the GT.
    vector<Point3D> mergedset;


    for (int j = 1; j <= nbTests; ++j){
        CALL_SUBTEST(( test_model(transforms, files, mergedset, j, i) ));
    }

    return EXIT_SUCCESS;
}
