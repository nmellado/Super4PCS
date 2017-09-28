Super4PCS
=========

An implementation of the Super 4-points Congruent Sets (Super 4PCS) algorithm presented in:

> Super 4PCS: Fast Global Pointcloud Registration via Smart Indexing
> Nicolas Mellado, Dror Aiger, Niloy J. Mitra
> Symposium on Geometry Processing 2014.

Authors: Nicolas Mellado, Dror Aiger

[stsimg-linux]: https://api.travis-ci.org/nmellado/Super4PCS.svg?branch=master
[stsimg-windw]: https://ci.appveyor.com/api/projects/status/reg4cmhn309w1s8k/branch/master?svg=true

| Linux  - MacOS  | Windows         |
| :----:          | :-----:         |
| ![stsimg-linux] | ![stsimg-windw] |

Documentation: [http://nmellado.github.io/Super4PCS/](http://nmellado.github.io/Super4PCS/)

Paper project page: [http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS](http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS)



## News
* 28th September 2017: Super4PCS [v1.1.3](https://github.com/nmellado/Super4PCS/releases/tag/v1.1.3) is out, providing a new Meshlab plugin, PCL wrapper, doxygen-generated documentation (online [here](http://nmellado.github.io/Super4PCS/)) and several bug fixes.
* 24th August 2017: Super4PCS [v1.1.2](https://github.com/nmellado/Super4PCS/releases/tag/v1.1.2) introduces Super4PCS CMake package, better testing (including on Windows plateform), cleaner code structure (moved demos away from library code).
* 10th August 2017: Super4PCS [v1.1.1](https://github.com/nmellado/Super4PCS/releases/tag/v1.1.1) introduces minor fixes from previous release.
* 10th August 2017: Super4PCS [v1.1.0](https://github.com/nmellado/Super4PCS/releases/tag/v1.1.0) is available, improving performances, fixing use as C++ library, and generalizing CI use.
* 17th July 2017: CI integration (Windows, MacOS and Linux) enabled. Currently, only compilation is tested, but performances monitoring will be added in upcoming release.
* 5th July 2017: Super4PCS got the SGP Software Award 2017 !
* 19th June 2017: Super4PCS v1.0.0-alpha is available.
* 18th May 2016: Super4PCS [v0.2.1-alpha](https://github.com/nmellado/Super4PCS/releases/tag/v0.2.1-alpha) is available, fixing a problematic crash introduced in previous release.
* 3rd May 2016: Super4PCS [v0.2-alpha](https://github.com/nmellado/Super4PCS/releases/tag/v0.2-alpha) is out !
* 23th March 2016: Super4PCS can now be compiled with Visual Studio 2015! Checkout the [Wiki](http://github.com/nmellado/Super4PCS/wiki) for more details.

## Paper Abstract
> Data acquisition in large-scale scenes regularly involves accumulating information across multiple scans. A common approach is to locally align scan pairs using Iterative Closest Point (ICP) algorithm (or its variants), but requires static scenes and small motion between scan pairs. This prevents accumulating data across multiple scan sessions and/or different acquisition modalities (e.g., stereo, depth scans). Alternatively, one can use a global registration algorithm allowing scans to be in arbitrary initial poses. The state-of-the-art global registration algorithm, 4PCS, however has a quadratic time complexity in the number of data points. This vastly limits its applicability to acquisition of large environments. We present Super 4PCS for global pointcloud registration that is optimal, i.e., runs in linear time (in the number of data points) and is also output sensitive in the complexity of the alignment problem based on the (unknown) overlap across scan pairs. Technically, we map the algorithm as an ‘instance problem’ and solve it efficiently using a smart indexing data organization. The algorithm is simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in significant speedup over alternative approaches and allows unstructured efficient acquisition of scenes at scales previously not possible. Complete source code and datasets are available for research use at http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.
