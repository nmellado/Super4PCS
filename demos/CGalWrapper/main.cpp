#include <Eigen/Core>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/property_map.h>
#include <CGAL/super4pcs.h>

#include "../demo-utils.h"

#include <utility>
#include <vector>
#include <fstream>

// Types
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point  = Kernel::Point_3;
using Vector = Kernel::Vector_3;

typedef CGAL::cpp11::tuple<Point, Vector> PNCI;
typedef CGAL::Nth_of_tuple_property_map<0, PNCI> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNCI> Normal_map;

using namespace GlobalRegistration;

// Align a rigid object to a scene with clutter and occlusions
int
main (int argc, char **argv)
{
  // Get input object and scene
  std::vector<PNCI> ppoints, qpoints;
  if (argc < 4)
  {
    std::cerr << "\nUsage: " << argv[0] << " -i input1.ply input2.ply\n";
    Demo::printParameterList();
    return (-1);
  }

  // Load object and scene
  std::cout << "Loading point clouds..." << std::endl;
  auto loadFile = [] (std::vector<PNCI>& buf, char* fname) {
    std::ifstream in(fname);
    return ( in &&
             CGAL::read_ply_points_with_properties
             (in,
              std::back_inserter (buf),
              CGAL::make_ply_point_reader (Point_map()),
              CGAL::make_ply_normal_reader (Normal_map())
              ));
  };

  if (! loadFile(ppoints, argv[0])  || ! loadFile(qpoints, argv[1]) )
  {
    std::cerr << "Error: cannot read input file." << std::endl;
    return (-1);
  }

  // Load Super4pcs parameters
  Demo::getArgs(argc, argv);

  GlobalRegistration::Match4PCSOptions opt;
  Demo::setOptionsFromArgs(opt);

  CGAL::Aff_transformation_3<Kernel> tr =
      CGAL::align<Kernel>(ppoints.begin(), ppoints.end(),
                  qpoints.begin(), qpoints.end(),
                  Point_map(),
                  opt);


  return (0);
}
