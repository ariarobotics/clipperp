// for CLIPPER (based on Belachew) optimization-based clique algorithm
#pragma once

#include <vector>
#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include "clipperplus/utils.h"


namespace clipperplus 
{

struct Params {
  // gradient descent parameters
  double tol_u = 1e-8; // (1e-8) stop innerloop when change in u < tol
  double tol_F = 1e-9; // (1e-9) stop innerloop when change in F < tol
  int maxiniters = 200; // max num of gradient ascent steps for each d
  int maxoliters = 1000; // max num of outer loop iterations to find d
  // line search parameters
  double beta = 0.5; // backtracking step size reduction, in (0, 1)
  double sigma = 0.01; // threshold of armijo condition
  
  /* maximum number of line search iters per grad step 
    note: do NOT choose large number (>30), or otherwise 
    alpha stepsize will become zero and line search will get stuck */
  int maxlsiters = 20; 
  double minalpha = 1e-9; // minimum value of alpha for line search
  double maxalpha = 2; // maximum value of alpha for line search

  double eps = 1e-9; // (1e-9) numerical threshold around 0 (default 1e-9)
};

std::vector<long> clique_optimization(const Eigen::MatrixXd& M, const Eigen::VectorXd& u0, const Params &params);

} 
