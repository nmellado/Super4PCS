Super4PCS
=========

Copyright 2014 Nicolas Mellado

Authors: Nicolas Mellado, Dror Aiger


##Related pages
* [Paper project page](geometry.cs.ucl.ac.uk/projects/2014/super4PCS)
* [Github project page](github.com/nmellado/Super4PCS)
* [Wiki](github.com/nmellado/Super4PCS/wiki)

Please read the *wiki* to get up to date informations about the algorithm and the development state.

##Overview
An implementation of the Super 4-points Congruent Sets (Super 4PCS) 
algorithm presented in:

>Super 4PCS: Fast Global Pointcloud Registration via Smart Indexing
>Nicolas Mellado, Dror Aiger, Niloy J. Mitra
>Symposium on Geometry Processing 2014.

Data acquisition in large-scale scenes regularly involves accumulating information across multiple scans. A common approach is to locally align scan pairs using Iterative Closest Point (ICP) algorithm (or its variants), but requires static scenes and small motion between scan pairs. This prevents accumulating data across multiple scan sessions and/or different acquisition modalities (e.g., stereo, depth scans). Alternatively, one can use a global registration algorithm allowing scans to be in arbitrary initial poses. The state-of-the-art global registration algorithm, 4PCS, however has a quadratic time complexity in the number of data points. This vastly limits its applicability to acquisition of large environments. We present Super 4PCS for global pointcloud registration that is optimal, i.e., runs in linear time (in 
the number of data points) and is also output sensitive in the complexity of the alignment problem based on the (unknown) overlap across scan pairs. Technically, we map the algorithm as an ‘instance problem’ and solve it efficiently using a smart indexing data organization. The algorithm is simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in significant speedup over alternative approaches and allows unstructured efficient acquisition of scenes at scales previously not possible. Complete source code and datasets are available for research use at geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.


##Usage
A good starting point is to call the enclosed script `./run-example.sh`.

More details about parameter setting will be added soon. In the meantime please contact us to get support on parameter tuning (see [Wiki](github.com/nmellado/Super4PCS/wiki)).

This code also implements the initial 4PCS algorithm, just append -x to your command line call to run it.


##Contact
Please feel free to contact us if you have any question regarding this source code (see the [project page](geometry.cs.ucl.ac.uk/projects/2014/super4PCS/) for contact information).

