/**
 * @file trampolines.h
 * @brief Trampolines to help overriding C++ in Python
 * @author Parker Lusk <plusk@mit.edu>
 * @date 16 May 2021
 */

#pragma once

#include <pybind11/pybind11.h>
#include "clipper/invariants/abstract.h"
#include "clipperplus/clipperplus_clique.h"

template <class InvariantBase = clipper::invariants::Invariant>
class PyInvariant : public InvariantBase {
public:
    using InvariantBase::InvariantBase; // Inherit constructors
};

template <class PairwiseInvariantBase = clipper::invariants::PairwiseInvariant>
class PyPairwiseInvariant : public PyInvariant<PairwiseInvariantBase> {
public:
    using PyInvariant<PairwiseInvariantBase>::PyInvariant; // Inherit constructors
    using Datum = clipper::invariants::Datum; // for convenience
    // trampoline for virtual function
    double operator()(const Datum& ai, const Datum& aj, const Datum& bi, const Datum& bj) override {
      pybind11::gil_scoped_acquire acquire; // Acquire GIL before calling Python code
      PYBIND11_OVERRIDE_PURE_NAME(double, PairwiseInvariantBase, "__call__", operator(), ai, aj, bi, bj);
    }
};

// Wrapper function for clipperplus_clique
std::tuple<long, std::vector<int>, int> clipperplus_clique_wrapper(
  const Eigen::MatrixXd& adj)
{
  long clique_size = 0;
  std::vector<int> clique;
  int certificate = 0;

  clipperplus::clipperplus_clique(adj, clique_size, clique, certificate);

  return std::make_tuple(clique_size, clique, certificate);
};