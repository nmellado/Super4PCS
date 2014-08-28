Super4PCS
=========

Copyright 2014 Nicolas Mellado

Authors: Nicolas Mellado, Dror Aiger


##Overview
An implementation of the Super 4-points Congruent Sets (Super 4PCS) 
algorithm presented in:

>Super 4PCS: Fast Global Pointcloud Registration via Smart Indexing
>Nicolas Mellado, Dror Aiger, Niloy J. Mitra
>Symposium on Geometry Processing 2014.

Data acquisition in large-scale scenes regularly involves accumulating information across multiple scans. A common approach is to locally align scan pairs using Iterative Closest Point (ICP) algorithm (or its variants), but requires static scenes and small motion between scan pairs. This prevents accumulating data across multiple scan sessions and/or different acquisition modalities (e.g., stereo, depth scans). Alternatively, one can use a global registration algorithm allowing scans to be in arbitrary initial poses. The state-of-the-art global registration algorithm, 4PCS, however has a quadratic time complexity in the number of data points. This vastly limits its applicability to acquisition of large environments. We present Super 4PCS for global pointcloud registration that is optimal, i.e., runs in linear time (in 
the number of data points) and is also output sensitive in the complexity of the alignment problem based on the (unknown) overlap across scan pairs. Technically, we map the algorithm as an ‘instance problem’ and solve it efficiently using a smart indexing data organization. The algorithm is simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in significant speedup over alternative approaches and allows unstructured efficient acquisition of scenes at scales previously not possible. Complete source code and datasets are available for research use at http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.

##Development state
I am currently working on the source code to clean it and define a stable API. More interesting updates will come soon !


##Compilation and usage
###Dependencies:
* [Eigen](http://eigen.tuxfamily.org/)
* [LibANN](http://www.cs.umd.edu/~mount/ANN/)
* [OpenCV](http://opencv.org/)
* [Chealpix](http://healpix.jpl.nasa.gov/html/csubnode4.htm), [Ubuntu package](http://packages.ubuntu.com/km/source/utopic/all/misc/chealpix), [Debian package](https://packages.debian.org/source/jessie/chealpix)
* [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html), a dependency of chealpix. Both ubuntu and debian provide a virtual package 'cfitsio-dev'.


###Compilation
The simplest solution is to use cmake, from the source directory:
```
mkdir build
cd build
cmake -DANN_DIR=/your/path/to/ann_1.1.2/ ..
```
For now ANN doesn't come with a cmake package, so you need to set its path by hand (this will be fixed later). Chealpix sources are now part of the Super4PCS repository and are compiled automatically.

Not recommended: you can also edit and use the included Makefile to compile using GCC with C++11 support enabled.

###Compilation error
You may encounter compilation issues when compiling with old versions of Eigen and C++11:
```
/somePath/eigen3/Eigen/src/Core/util/Macros.h:252:35: error: unable to find string literal operator ‘operator""X’
 #define EIGEN_ASM_COMMENT(X)  asm("#"X)
                                    ^
/somePath/eigen3/Eigen/src/Core/products/GeneralBlockPanelKernel.h:604:1: note: in expansion of macro ‘EIGEN_ASM_COMMENT’
 EIGEN_ASM_COMMENT("mybegin2");
 ^
```
You can easily fix it by:
* update and use an up-to-date version of [Eigen](http://eigen.tuxfamily.org/). You may need to call
```
cmake -DEIGEN3_INCLUDE_DIR=/your/path/to/eigen/ ..
```
* Not recommended: modify one line in Eigen sources (Instructions [here](https://sourceforge.net/p/pagmo/mailman/message/30074799/)).


###Usage
A good starting point is to call the enclosed script `./run-example.sh`.

##Contact
Please feel free to contact us if you have any question regarding this source code (see the [project page](http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/) for contact information).

