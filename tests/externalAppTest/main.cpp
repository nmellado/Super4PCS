#include "super4pcs/algorithms/super4pcs.h"
#include "super4pcs/io/io.h"
#include "super4pcs/utils/geometry.h"

#include <Eigen/Dense>


int main(int argc, char **argv) {
  using namespace GlobalRegistration;
  using namespace std;

  vector<Point3D> set1, set2;
  vector<Eigen::Matrix2f> tex_coords1, tex_coords2;
  vector<typename Point3D::VectorType> normals1, normals2;
  vector<tripple> tris1, tris2;
  vector<std::string> mtls1, mtls2;

  IOManager iomanager;

  // dummy call, to test symbols accessibility
  iomanager.ReadObject("", set1, tex_coords1, normals1, tris1, mtls1);

  // check availability of the Utils functions
  if (tris1.size() == 0)
    Utils::CleanInvalidNormals(set1, normals1);

  // Our matcher.
  Match4PCSOptions options; //TODO : MatchOptions

  // Set parameters.
  Match4PCSBase<>::MatrixType mat;
  double overlap (1);
  options.configureOverlap(overlap);

  typename Point3D::Scalar score = 0;

  constexpr Utils::LogLevel loglvl = Utils::Verbose;
  using TrVisitorType = Match4PCSBase<>::DummyTransformVisitor;
  using SamplerType   = Match4PCSBase<>::DefaultSampler;
  Utils::Logger logger(loglvl);

  MatchSuper4PCS matcher(options, logger);
  score = matcher.ComputeTransformation<SamplerType,TrVisitorType>(set1, &set2, mat);

  logger.Log<Utils::Verbose>( "Score: ", score );

  iomanager.WriteMatrix("output.map", mat.cast<double>(), IOManager::POLYWORKS);

  return 0;
}

