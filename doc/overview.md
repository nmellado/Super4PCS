# Overview {#mainpage}

## About Super4PCS

The Super4PCS library is a set C++ libraries, standalone applications and plugins released under the terms of the *APACHE V2 licence*, which makes it free for commercial and research use. 

It provides state of the art global registration techniques for 3d pointclouds. Still under active development, the library is currently improved by adding wrappers to existing software, and by stabilizing the API.

In a nutshell, it provides:
* two core packages: [algorithms](@ref GlobalRegistration::algorithms) and [accelerators](@ref GlobalRegistration::accelerators),
* two convenience packages: [io](@ref GlobalRegistration::io) and [utils](@ref GlobalRegistration::utils),
* several [demonstrators](@ref Super4PCS-Demos):
   * a [command line application](@ref Super4PCS-App),
   * wrappers to use the algorithms with the [Point Cloud Library](@ref Super4PCS-PCL) and [CGAL](@ref Super4PCS-CGAL) (coming soon)
   * a [meshlab plugin](@ref Super4PCS-Meshlab)
* a [test suite](@ref GlobalRegistration::Testing)

## Compilation status
[stsimg-linux]: https://api.travis-ci.org/nmellado/Super4PCS.svg?branch=master
[stsimg-windw]: https://ci.appveyor.com/api/projects/status/reg4cmhn309w1s8k/branch/master?svg=true

[Travis]: https://travis-ci.org/nmellado/Super4PCS "Travis"
[AppVeyor]: https://ci.appveyor.com/project/nmellado/super4pcs "AppVeyor"

| Linux  \& MacOS | Windows         |
| :----:          | :-----:         |
| ![stsimg-linux] | ![stsimg-windw] |
| [Travis]        | [AppVeyor]      |


## Credits
* Nicolas Mellado: conception, implementation and examples
* Dror Aiger: 4pcs implementation
* Niloy Mitra: 4pcs implementation


### Citation
```
﻿@misc {Super4PCSLibrary,
    author = {Mellado, Nicolas and Aiger, Dror and Mitra, Niloy J.},
    title = {Super 4PCS Library},
    howpublished = {https://github.com/nmellado/Super4PCS},
    year = {2017},
}
```
