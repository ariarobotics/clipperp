# CLIPPER+: A Fast Maximal Clique Algorithm for Robust Global Registration 

CLIPPER+ is an algorithm for finding maximal cliques in unweighted graphs for outlier-robust global registration. The registration problem can be 
formulated as a graph and solved by finding its maximum clique. This formulation leads to extreme robustness to outliers; however, finding the 
[maximum clique](https://en.wikipedia.org/wiki/Clique_problem) is an NP-hard problem, and therefore approximation is required in practice for large-size problems.

The performance of an approximation algorithm is evaluated by its computational complexity (the lower the runtime, the
better) and solution accuracy (how close the solution is to the maximum clique). Accordingly, the main contribution of 
CLIPPER+ is outperforming the state-of-the-art in accuracy while maintaining a relatively low runtime. CLIPPER+ builds
on prior work (PMC <sup>[[2]](#2)</sup> and CLIPPER <sup>[[3]](#3)</sup>) and prunes the graph by removing vertices that have a small core number and cannot be a part of the maximum clique. This will result in a smaller graph, on which the maximum clique can be estimated considerably faster.

Note
---------
We are actively optimizing the code and adding new features to this repo.
 

Citation
---------
If you find this code useful in your research, we ask that you cite:

- K. Fathian, T. Summers. "CLIPPER+: A Fast Maximal Clique Algorithm for Robust Global Registration," in IEEE Robotics and Automation Letters, 2024. ([Paper](https://arxiv.org/pdf/2402.15464.pdf))

```bibtex
@inproceedings{fathian2024clipper+,
  title={CLIPPER+: A Fast Maximal Clique Algorithm for Robust Global Registration},
  author={Fathian, Kaveh and Summers, Tyler},
  journal={IEEE Robotics and Automation Letters},
  year={2024},
  publisher={IEEE}
}
```

Installation
---------
### Compilation

Ensure that cmake is installed on your system.
Upon cloning the repository, run the following commands:

```bash
$ mkdir build
$ cd ./build
$ cmake ..
$ make
```

If succesful, tests can be run with the command `./test/tests`

### Build Configurations

The following cmake options are available when building CLIPPER+ :

| Option                  | Description                                      K. Fathian, T. Summers, "CLIPPER+: A Graph-Theoretic Framework for Robust Data Association," in *Proc. IEEE Int. Conf. on Robotics and Automation (ICRA)*, 2021, pp. 13828-13834. DOI: 10.1109/ICRA48506.2021.9561069. ([Paper](https://arxiv.org/pdf/2011.10202.pdf)) ([Video](https://youtu.be/QYLHueMhShY)) ([Repository](https://github.com/mit-acl/clipper))
ndings](#python-bindings) for CLIPPER+. | `OFF` |
| `BUILD_BINDINGS_MATLAB` | Creates [MATLAB bindings](#matlab-bindings) for CLIPPER+. | `OFF` |
| `DEBUG_FLAG`     | Enables debugging output | `OFF` | 
| `DEBUG_TIMING_FLAG` | Enables timing report | `OFF` |
| `DEBUG_OPTIM_FLAG ` | Enable debugging clipper optimization | `OFF` |
| `DEBUG_BUILD_PMC_HEU` | Builds PMC heuristic | `OFF` |
| `BUILD_TESTS` | Builds Google testsuites | `ON` |

You can customize these preferences by employing the -D flag during the cmake execution for each individual flag. For instance, you can use the command `cmake -D BUILD_BINDINGS_PYTHON=ON -D DEBUG_FLAG=ON` to enable Python bindings and enable basic debugging output.

### Python Bindings

If Python bindings are built using the `BUILD_BINDINGS_PYTHON` flag (See [Build Configurations](#build-configurations) above), then the generated clipper Python module will need to be installed before it can be imported in any Python 3 program. This can be done by running make:

```bash
$ cd ./build
$ make pip-install
```

Or, if the user wishes to install the module using pip directly:
```bash
$ python3 -m pip install build/bindings/python
```
A Jupyter Notebook containing a python example can be found in the [`test/python`](test/python) folder
### MATLAB Bindings

If MATLAB is installed on your computer and MATLAB bindings are requested using the cmake flag `BUILD_BINDINGS_MATLAB=ON` (See [Build Configurations](#build-configurations) above), then cmake will create a bindings/matlab folders with additional directives. By running these commands inside your build folder:
```bashs
$ cd ./bindings/matlab
$ make
```
cmake will attempt to find your MATLAB installation and subsequently generate a set of MEX files so that CLIPPER+ can be used from your MATLAB program.

MATLAB test suites can be found in the [`test/matlab`](test/matlab) folder. For future projects, ensure that all the generated MEX files are located in your current MATLAB path.

### Including as a shared library

A simple way to include `clipperplus` as a shared library in another C++ project is via `cmake`. This method will automatically clone and build `clipperplus`, making the resulting library accessible in your main project. In your project's `CMakeLists.txt` you can add

```cmake
set(CLIPPERPLUS_DIR "${CMAKE_CURRENT_BINARY_DIR}/clipperplus-download" CACHE INTERNAL "CLIPPERPLUS build dir" FORCE)
set(DEBUG_FLAG OFF CACHE BOOL "")
set(DEBUG_TIMING_FLAG OFF CACHE BOOL "")
set(DEBUG_OPTIM_FLAG OFF CACHE BOOL "")
set(BUILD_PMC_HEU OFF CACHE BOOL "")
set(BUILD_BINDINGS_MATLAB OFF CACHE BOOL "")
set(BUILD_BINDINGS_PYTHON OFF CACHE BOOL "")
set(BUILD_TESTS OFF CACHE BOOL "")
configure_file(cmake/clipperplus.cmake.in ${CLIPPERPLUS_DIR}/CMakeLists.txt IMMEDIATE @ONLY)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${CLIPPERPLUS_DIR})
execute_process(COMMAND "${CMAKE_COMMAND}" --build . WORKING_DIRECTORY ${CLIPPERPLUS_DIR})
add_subdirectory(${CLIPPERPLUS_DIR}/src ${CLIPPERPLUS_DIR}/build)
```

where `cmake/clipperplus.cmake.in` looks like

```cmake

cmake_minimum_required(VERSION 3.10)
project(clipperplus-download NONE)

include(ExternalProject)
ExternalProject_Add(clipperplus
    GIT_REPOSITORY      "https://github.com/ariarobotics/clipperp"
    GIT_TAG             main
    SOURCE_DIR          "${CMAKE_CURRENT_BINARY_DIR}/src"
    BINARY_DIR          "${CMAKE_CURRENT_BINARY_DIR}/build"
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     ""
    TEST_COMMAND        ""
)
```

Then, you can link your project with `clipperplus` using the syntax `target_link_libraries(yourproject clipperplus)`.


## Examples

See the [`examples`](examples) folder to see demos of CLIPPER+ being used for finding maximal cliques and pointcloud registration.


## Citations

<a id="1">[1]</a> K. Fathian, T. Summers. "CLIPPER+: A Fast Maximal Clique Algorithm for Robust Global Registration," in *IEEE Robotics and Automation Letters (RAL)*, 2024. ([Paper](https://arxiv.org/pdf/2402.15464.pdf))

<a id="1">[2]</a> R. A. Rossi, D. F. Gleich, A. H. Gebremedhin, M. M. Patwary, "A Fast Parallel Maximum Clique Algorithm for Large Sparse Graphs and Temporal Strong Components," *arXiv preprint arXiv:1302.6256*, 2013. ([Repository](https://github.com/ryanrossi/pmc))

<a id="1">[3]</a> P. C. Lusk, K. Fathian, J. P. How, "CLIPPER: A Graph-Theoretic Framework for Robust Data Association," in *Proc. IEEE Int. Conf. on Robotics and Automation (ICRA)*, 2021, pp. 13828-13834. DOI: 10.1109/ICRA48506.2021.9561069. ([Paper](https://arxiv.org/pdf/2011.10202.pdf)) ([Video](https://youtu.be/QYLHueMhShY)) ([Repository](https://github.com/mit-acl/clipper))

Copyright 2024, Kaveh Fathian. All rights reserved.
		
