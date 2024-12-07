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

// #include "clipper/clipper.h"
// #include "clipper/utils.h"

#include "wrapper.h"


namespace py = pybind11;
using namespace pybind11::literals;

// ----------------------------------------------------------------------------

PYBIND11_MODULE(clipperpluspy, m)
{
  m.doc() = "CLIPPER+ is an algorithm for finding maximal cliques in unweighted graphs for outlier-robust global registration.";
  m.attr("__version__") = CLIPPERPLUS_VERSION;

  m.def("clipperplus_clique", &Wrapper::clipperplus_clique_wrapper,
    "adj"_a,
    "Find the densest subgraph of a weighted adjacency matrix.");
  m.def("find_heuristic_clique", &Wrapper::find_heuristic_clique_wrapper,
    "adj"_a, "clique"_a,
    "Find a heuristic maximum clique in a graph.");
  m.def("clique_optimization", &Wrapper::clique_optimization_wrapper,
    "M"_a, "u0"_a,
    "Run original clipper on pruned graph");
}
