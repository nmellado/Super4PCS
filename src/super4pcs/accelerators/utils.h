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
// Technically, we map the algorithm as an 'instance problem' and solve it
// efficiently using a smart indexing data organization. The algorithm is
// simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in
// significant speedup over alternative approaches and allows unstructured
// efficient acquisition of scenes at scales previously not possible. Complete
// source code and datasets are available for research use at
// http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.


#ifndef _SUPER4PCS_ACCELERATORS_UTILS_H_
#define _SUPER4PCS_ACCELERATORS_UTILS_H_


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace GlobalRegistration{
namespace Utils{

//! \brief Compile time pow
template<typename baseT, typename expoT>
constexpr baseT POW(baseT base, expoT expo)
{
    return (expo != 0 )? base * POW(base, expo -1) : 1;
}


namespace internal{

/*!
 * \brief Helper structure used to (conditionnaly) validate indices.
 * The template parameter enable allows to control the validation, without
 * overload when the validation is disabled.
 *
 * Throw std::out_of_range exception when the validation is enabled and an
 * invalid case is detected
 */
template <bool enable>
struct IndexValidator{
    /*!
     * \brief Check if n is within [0:gsize] if enable=true
     * \throw std::out_of_range
     */
    template <class IndexT, class SizeT>
    static inline
    constexpr IndexT validate(const IndexT& n,
                              const SizeT& gsize);

    /*!
     * \brief Test used internally to validate the input index
     */
    template <class IndexT, class SizeT>
    static inline
    constexpr bool _test(const IndexT& n, const SizeT& gsize){
        return n < gsize;
    }
};

template <>
template <class IndexT, class SizeT>
constexpr IndexT
IndexValidator<true>::validate(const IndexT& n,
                               const SizeT&  gsize)
{
    return IndexValidator<true>::_test(n, gsize) ?
                n :
                throw std::out_of_range(
                    std::string("[IndexValidator] ") +
                    std::to_string(n)                +
                    std::string(" bigger than ")     +
                    std::to_string(gsize) );
}

template <>
template <class IndexT, class SizeT>
constexpr IndexT
IndexValidator<false>::validate(const IndexT& n,
                                const SizeT& /*gsize*/)
{
    return n;
}
} // namespace GlobalRegistration::Utils::internal

/*!
 * \brief Convert a normalized n-d vector to a linear index in a uniform regular grid
 * This function is recursive, and unrolled at compile time (loop over n).
 *
 * \param coord Input coordinates defined in the normalized n-hypercube.
 * \param cdim  Working dimension, must be called with n.
 * \param gsize Dimension of the grid, must be consistent in all dimensions
 *
 * \tparam validate Enable or disable the range validation
 * \tparam PointT Type of the input points (deduced)
 * \tparam IndexT Index type (deduced)
 * \tparam IndexT Size type (deduced)
 *
 * \see internal::IndexValidator for the validation procedure
 */
template<bool validate, class ndIndexT, class IndexT, class SizeT>
constexpr inline IndexT
UnrollIndexLoop(const ndIndexT& coord,
                IndexT        cdim,
                SizeT         gsize){
  return (cdim != 0)
    ? ( internal::IndexValidator<validate>::validate(IndexT(coord[cdim]), gsize)*POW(gsize, cdim) +
        UnrollIndexLoop<validate>(coord, cdim-1, gsize) )
    : internal::IndexValidator<validate>::validate(IndexT(coord[cdim]), gsize);
}


/*!
 * \brief Convert a normalized n-d vector to a linear index in a uniform regular
 * grid, moved by moved by an offset defined as a integer move in the n-d grid.
 * This function is recursive, and unrolled at compile time (loop over n). In
 * addition, it allows to offset the input coordinates.
 *
 * \param coord Input coordinates defined in the normalized n-hypercube.
 * \param cdim  Working dimension, must be called with n.
 * \param gsize Dimension of the grid, must be consistent in all dimensions
 *
 * \see UnrollIndexLoop<class PointT>
 */
template<bool validate, class ndIndexT, class IndexT, class SizeT>
constexpr inline IndexT
UnrollIndexLoop(const ndIndexT& coord,
                const ndIndexT& offset,
                IndexT          cdim,
                SizeT           gsize){
  return (cdim != 0)
    ? ( internal::IndexValidator<validate>::validate(IndexT(offset[cdim]+std::floor(coord[cdim])), gsize)*POW(gsize, cdim) +
        UnrollIndexLoop<validate>(coord, offset, cdim-1, gsize) )
    : internal::IndexValidator<validate>::validate(IndexT(offset[cdim]+std::floor(coord[cdim])), gsize);
}

} //namespace GlobalRegistration::Utils
} //namespace Super4PCS


#endif
