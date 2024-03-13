/**
 * @file py_clipper.cpp
 * @brief Python bindings for CLIPPER
 * @author Parker Lusk <plusk@mit.edu>
 * @date 28 January 2021
 */

#include <cstdint>
#include <sstream>

#include <Eigen/Dense>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "clipper/clipper.h"
#include "clipper/utils.h"

#include "trampolines.h"


namespace py = pybind11;
using namespace pybind11::literals;

// ----------------------------------------------------------------------------

PYBIND11_MODULE(clipperpluspy, m)
{
  m.doc() = "CLIPPER+ is an algorithm for finding maximal cliques in unweighted graphs for outlier-robust global registration.";
  m.attr("__version__") = CLIPPERPLUS_VERSION;

  m.def("clipperplus_clique", &clipperplus_clique_wrapper,
    "adj"_a,
    "Find the densest subgraph of a weighted adjacency matrix.");
}
