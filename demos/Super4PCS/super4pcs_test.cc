#include "super4pcs/algorithms/4pcs.h"
#include "super4pcs/algorithms/super4pcs.h"
#include "super4pcs/io/io.h"
#include "super4pcs/utils/geometry.h"

#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <string>

#include "../demo-utils.h"

#define sqr(x) ((x) * (x))

using namespace std;
using namespace GlobalRegistration;
using namespace GlobalRegistration::Demo;


static inline void printS4PCSParameterList(){
    fprintf(stderr, "\t[ -r result_file_name (%s) ]\n", output.c_str());
    fprintf(stderr, "\t[ -m output matrix file (%s) ]\n", outputMat.c_str());
    fprintf(stderr, "\t[ -x (use 4pcs: false by default) ]\n");
    fprintf(stderr, "\t[ --sampled1 (output sampled cloud 1 -- debug+super4pcs only) ]\n");
    fprintf(stderr, "\t[ --sampled2 (output sampled cloud 2 -- debug+super4pcs only) ]\n");
}

struct TransformVisitor {
    inline void operator()(
            float fraction,
            float best_LCP,
            Eigen::Ref<typename  Match4PCSBase<>::MatrixType> /*transformation*/) const {
      if (fraction >= 0)
        {
          printf("done: %d%c best: %f                  \r",
               static_cast<int>(fraction * 100), '%', best_LCP);
          fflush(stdout);
        }
    }
    constexpr bool needsGlobalTransformation() const { return false; }
};

int main(int argc, char **argv) {
  using namespace GlobalRegistration;

  vector<Point3D> set1, set2;
  vector<Eigen::Matrix2f> tex_coords1, tex_coords2;
  vector<typename Point3D::VectorType> normals1, normals2;
  vector<tripple> tris1, tris2;
  vector<std::string> mtls1, mtls2;

  // Match and return the score (estimated overlap or the LCP).
  typename Point3D::Scalar score = 0;

  constexpr Utils::LogLevel loglvl = Utils::Verbose;
  using SamplerType   = GlobalRegistration::Sampling::UniformDistSampler;
  using TrVisitorType = typename std::conditional <loglvl==Utils::NoLog,
                            typename Match4PCSBase<>::DummyTransformVisitor,
                            TransformVisitor>::type;
  SamplerType sampler;
  TrVisitorType visitor;
  Utils::Logger logger(loglvl);

  /// TODO Add proper error codes
  if(argc < 4){
      Demo::printUsage(argc, argv);
      exit(-2);
  }
  if(int c = Demo::getArgs(argc, argv) != 0)
  {
    Demo::printUsage(argc, argv);
    printS4PCSParameterList();
    exit(std::max(c,0));
  }

  // prepare matcher ressources
  Match4PCSOptions options; //TODO : MatchOptions
  using MatrixType = typename Match4PCSBase<>::MatrixType;
  MatrixType mat (MatrixType::Identity());

  if(! Demo::setOptionsFromArgs(options, logger))
  {
    exit(-3);
  }

  // load data
  IOManager iomananger;

  // Read the inputs.
  if (!iomananger.ReadObject((char *)input1.c_str(), set1, tex_coords1, normals1, tris1,
                  mtls1)) {
    logger.Log<Utils::ErrorReport>("Can't read input set1");
    exit(-1);
  }

  if (!iomananger.ReadObject((char *)input2.c_str(), set2, tex_coords2, normals2, tris2,
                  mtls2)) {
    logger.Log<Utils::ErrorReport>("Can't read input set2");
    exit(-1);
  }

  // clean only when we have pset to avoid wrong face to point indexation
  if (tris1.size() == 0)
    Utils::CleanInvalidNormals(set1, normals1);
  if (tris2.size() == 0)
    Utils::CleanInvalidNormals(set2, normals2);

  try {

      if (use_super4pcs) {
          MatchSuper4PCS matcher(options, logger); //TODO : MatchBase<FunctorSuper4pcs> matcher(options, logger); ?
          logger.Log<Utils::Verbose>( "Use Super4PCS" );
          score = matcher.ComputeTransformation(set1, &set2, mat, sampler, visitor);

          if(! outputSampled1.empty() ){
              logger.Log<Utils::Verbose>( "Exporting Sampled cloud 1 to ",
                                          outputSampled1.c_str(),
                                          " ..." );
              iomananger.WriteObject((char *)outputSampled1.c_str(),
                                     matcher.getFirstSampled(),
                                     vector<Eigen::Matrix2f>(),
                                     vector<typename Point3D::VectorType>(),
                                     vector<tripple>(),
                                     vector<string>());
              logger.Log<Utils::Verbose>( "Export DONE" );
          }
          if(! outputSampled2.empty() ){
              logger.Log<Utils::Verbose>( "Exporting Sampled cloud 2 to ",
                                          outputSampled2.c_str(),
                                          " ..." );
              iomananger.WriteObject((char *)outputSampled2.c_str(),
                                     matcher.getSecondSampled(),
                                     vector<Eigen::Matrix2f>(),
                                     vector<typename Point3D::VectorType>(),
                                     vector<tripple>(),
                                     vector<string>());
              logger.Log<Utils::Verbose>( "Export DONE" );
          }
      }
      else {
          Match4PCS matcher(options, logger); //TODO : MatchBase<Functor4pcs> matcher(options, logger); ?
          logger.Log<Utils::Verbose>( "Use old 4PCS" );
          score = matcher.ComputeTransformation(set1, &set2, mat, sampler, visitor);
      }

  }
  catch (const std::exception& e) {
      logger.Log<Utils::ErrorReport>( "[Error]: " , e.what() );
      logger.Log<Utils::ErrorReport>( "Aborting with code -2 ..." );
      return -2;
  }
  catch (...) {
      logger.Log<Utils::ErrorReport>( "[Unknown Error]: Aborting with code -3 ..." );
      return -3;
  }

  logger.Log<Utils::Verbose>( "Score: ", score );
  logger.Log<Utils::Verbose>( "(Homogeneous) Transformation from ",
                              input2.c_str(),
                              " to ",
                              input1.c_str(),
                              ": \n",
                              mat);


  if(! outputMat.empty() ){
      logger.Log<Utils::Verbose>( "Exporting Matrix to ",
                                  outputMat.c_str(),
                                  "..." );
      iomananger.WriteMatrix(outputMat, mat.cast<double>(), IOManager::POLYWORKS);
      logger.Log<Utils::Verbose>( "Export DONE" );
  }

  if (! output.empty() ){

      logger.Log<Utils::Verbose>( "Exporting Registered geometry to ",
                                  output.c_str(),
                                  "..." );
      iomananger.WriteObject((char *)output.c_str(),
                             set2,
                             tex_coords2,
                             normals2,
                             tris2,
                             mtls2);
      logger.Log<Utils::Verbose>( "Export DONE" );
  }

  return 0;
}
