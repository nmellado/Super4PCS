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
// Technically, we map the algorithm as an ‘instance problem’ and solve it
// efficiently using a smart indexing data organization. The algorithm is
// simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in
// significant speedup over alternative approaches and allows unstructured
// efficient acquisition of scenes at scales previously not possible. Complete
// source code and datasets are available for research use at
// http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.

// Part of this file has been adapted from the Eigen library.


#ifndef _SUPER4PCS_TESTING_H_
#define _SUPER4PCS_TESTING_H_

#include <iostream>
#include <vector>
#include <cerrno>
#include <cstdlib>
#include <sstream>
#include <ctime>

#include "super4pcs/shared4pcs.h"

#define DEFAULT_REPEAT 10

#define SUPER4PCS_PP_MAKE_STRING2(S) #S
#define SUPER4PCS_PP_MAKE_STRING(S) SUPER4PCS_PP_MAKE_STRING2(S)

namespace GlobalRegistration {
namespace Testing {

//! Matcher class providing public visilibility of internal routines
template <class BaseMatcher>
class TestMatcher : public BaseMatcher {
public:
    using Scalar                = typename BaseMatcher::Scalar;
    using MatrixType            = typename BaseMatcher::MatrixType;
    using PairsVector           = typename BaseMatcher::PairsVector;
    using DefaultSampler        = typename BaseMatcher::DefaultSampler;
    using DummyTransformVisitor = typename BaseMatcher::DummyTransformVisitor;

    using BaseMatcher::BaseMatcher;

    template < typename Sampler = DefaultSampler,
               typename Visitor = DummyTransformVisitor>
    inline Scalar
    ComputeTransformation(const std::vector<Point3D>& P,
                          std::vector<Point3D>* Q,
                          Eigen::Ref<MatrixType> transformation,
                          const Sampler& s = Sampler(),
                          const Visitor& v = Visitor()){
        return BaseMatcher::ComputeTransformation(P, Q,
                                                  transformation,
                                                  s,
                                                  v);
    }

    template < typename Sampler = DefaultSampler>
    inline void init(const std::vector<Point3D>& P,
                     const std::vector<Point3D>& Q,
                     const Sampler& sampler = Sampler())
    { BaseMatcher::init(P,Q, sampler); }

    inline bool SelectQuadrilateral(Scalar &inv1, Scalar &inv2,
                                    int& base1, int& base2,
                                    int& base3, int& base4)
    {
        return BaseMatcher::SelectQuadrilateral(inv1, inv2,
                                                base1, base2, base3, base4);
    }

    inline const std::vector<Point3D>& base3D() const
    { return BaseMatcher::base3D(); }

    virtual void
    ExtractPairs( Scalar pair_distance,
                  Scalar pair_normals_angle,
                  Scalar pair_distance_epsilon, int base_point1,
                  int base_point2,
                  PairsVector* pairs) const
    {
        BaseMatcher::ExtractPairs(pair_distance,
                                  pair_normals_angle,
                                  pair_distance_epsilon,
                                  base_point1,
                                  base_point2,
                                  pairs);
    }

    virtual bool
    FindCongruentQuadrilaterals(Scalar inv1, Scalar inv2,
                                Scalar distance_threshold1,
                                Scalar distance_threshold2,
                                const PairsVector& P_pairs,
                                const PairsVector& Q_pairs,
                                std::vector<Quadrilateral>* quadrilaterals) const
    {
        return BaseMatcher::FindCongruentQuadrilaterals(inv1, inv2,
                                                        distance_threshold1,
                                                        distance_threshold2,
                                                        P_pairs,
                                                        Q_pairs,
                                                        quadrilaterals);
    }

    inline bool TryCongruentSet(int base_id1,
                                int base_id2,
                                int base_id3,
                                int base_id4,
                                const std::vector<Quadrilateral> &congruent_quads,
                                size_t &nbCongruent){
        return BaseMatcher::TryCongruentSet(base_id1, base_id2, base_id3, base_id4,
                                            congruent_quads,
                                            nbCongruent);
    }
};


static inline void
generateSphereCloud (std::vector<Point3D>& cloud,
                     size_t len){
  cloud.resize(len);
  for(auto& p : cloud)
  {
    typename Point3D::VectorType eigenpos =
             Point3D::VectorType::Random().normalized();
    p = Point3D(eigenpos);
    p.normalize();
  }
}


// extract pairs using brute force
template <typename Scalar, typename PairsVector>
static inline void
extractPairs( Scalar pair_distance,
              Scalar pair_distance_epsilon,
              const std::vector<Point3D>& cloud,
              PairsVector& pairs){
    pairs.clear();
    pairs.reserve(2 * cloud.size());

    // extract pairs using full search
    // Go over all ordered pairs in Q.
    for (int j = 0; j < cloud.size(); ++j) {
        const auto& p = cloud[j].pos();
        for (int i = j + 1; i < cloud.size(); ++i) {
            const auto& q = cloud[i].pos();
            const Scalar distance = (q - p).norm();
            if (std::abs(distance - pair_distance) <= pair_distance_epsilon) {
                pairs.emplace_back(j, i);
                pairs.emplace_back(i, j);
            }
        }
    }
}




static std::vector<std::string> g_test_stack;
static int g_repeat;
static unsigned int g_seed;
static bool g_has_set_repeat, g_has_set_seed;

void verify_impl(bool condition, const char *testname, const char *file, int line, const char *condition_as_string)
{
  if (!condition)
  {
    std::cerr << "Test " << testname << " failed in " << file << " (" << line << ")"
      << std::endl << "    " << condition_as_string << std::endl;
    abort();
  }
}

#define VERIFY(a) GlobalRegistration::Testing::verify_impl(a, \
                  GlobalRegistration::Testing::g_test_stack.back().c_str(), __FILE__, \
                  __LINE__, SUPER4PCS_PP_MAKE_STRING(a))

#if defined (_MSC_VER) && ! defined (__INTEL_COMPILER)
#  define MYPRAGMA(X) __pragma(X)
#else
#  define MYPRAGMA(X) _Pragma(X)
#endif

#define CALL_SUBTEST(FUNC) do { \
    MYPRAGMA("omp critical") \
    { GlobalRegistration::Testing::g_test_stack.push_back(SUPER4PCS_PP_MAKE_STRING(FUNC)); }\
    FUNC; \
    MYPRAGMA("omp critical") \
    { GlobalRegistration::Testing::g_test_stack.pop_back(); } \
  } while (0)

inline void set_repeat_from_string(const char *str)
{
  errno = 0;
  g_repeat = int(strtoul(str, 0, 10));
  if(errno || g_repeat <= 0)
  {
    std::cout << "Invalid repeat value " << str << std::endl;
    exit(EXIT_FAILURE);
  }
  g_has_set_repeat = true;
}

inline void set_seed_from_string(const char *str)
{
  errno = 0;
  g_seed = int(strtoul(str, 0, 10));
  if(errno || g_seed == 0)
  {
    std::cout << "Invalid seed value " << str << std::endl;
    exit(EXIT_FAILURE);
  }
  g_has_set_seed = true;
}

static bool init_testing(int argc, const char *argv[])
{
  g_has_set_repeat = false;
  g_has_set_seed   = false;
  bool need_help   = false;

  for(int i = 1; i < argc; i++)
  {
    if(argv[i][0] == 'r')
    {
      if(g_has_set_repeat)
      {
        std::cout << "Argument " << argv[i] << " conflicting with a former argument" << std::endl;
        return 1;
      }
      set_repeat_from_string(argv[i]+1);
    }
    else if(argv[i][0] == 's')
    {
      if(g_has_set_seed)
      {
        std::cout << "Argument " << argv[i] << " conflicting with a former argument" << std::endl;
        return false;
      }
        set_seed_from_string(argv[i]+1);
    }
    else
    {
      need_help = true;
    }
  }

  if(need_help)
  {
    std::cout << "This test application takes the following optional arguments:" << std::endl;
    std::cout << "  rN     Repeat each test N times (default: " << DEFAULT_REPEAT << ")" << std::endl;
    std::cout << "  sN     Use N as seed for random numbers (default: based on current time)" << std::endl;
    std::cout << std::endl;
    std::cout << "If defined, the environment variables SUPER4PCS_REPEAT and SUPER4PCS_SEED" << std::endl;
    std::cout << "will be used as default values for these parameters." << std::endl;
    return false;
  }

  char *env_SUPER4PCS_REPEAT = getenv("SUPER4PCS_REPEAT");
  if(!g_has_set_repeat && env_SUPER4PCS_REPEAT)
    set_repeat_from_string(env_SUPER4PCS_REPEAT);
  char *env_SUPER4PCS_SEED = getenv("SUPER4PCS_SEED");
  if(!g_has_set_seed && env_SUPER4PCS_SEED)
    set_seed_from_string(env_SUPER4PCS_SEED);

  if(!g_has_set_seed) g_seed = (unsigned int) time(NULL);
  if(!g_has_set_repeat) g_repeat = DEFAULT_REPEAT;

  std::cout << "Initializing random number generator with seed " << g_seed << std::endl;
  std::stringstream ss;
  ss << "Seed: " << g_seed;
  g_test_stack.push_back(ss.str());
  srand(g_seed);
  std::cout << "Repeating each test " << g_repeat << " times" << std::endl;

  return true;
}

} // namespace Testing
} // namespace Super4PCS

#endif // _SUPER4PCS_TESTING_H_
