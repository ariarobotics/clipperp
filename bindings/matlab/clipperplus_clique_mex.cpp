/**
 * @file clipperplus_clique_mex.cpp
 * @brief MATLAB/MEX binding for running CLIPPER+ 
 * @author Kaveh Fathian <kavehfathian@gmail.com>
 */
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>

#include <mex.h>
#include <Eigen/Dense>

#include "clipperplus/clipperplus_clique.h"
#include "mexutils.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
// nlhs   number of expected outputs
// plhs   array to be populated by outputs (data passed back to matlab)
// nrhs   number of inputs
// prhs   array poplulated by inputs (data passed from matlab)

bool useSparse = false;   // as determined by input adj
Eigen::MatrixXd adj;        // graph adjacency matrix
Eigen::SparseMatrix<double> adjs;  // sparse adjacency matrix

if (nrhs == 1) {
  if (mxIsSparse(prhs[0])) {
    adjs = mexMatrixToEigenSparse(prhs[0]);
    useSparse = true;
  } else if (!mxIsSparse(prhs[0])) {
    mexMatrixToEigen(prhs[0], &adj);
  } else {
    mexErrMsgIdAndTxt("clipperplus:nargin", "adj must be sparse or full matrix.");
  }
} else {
  mexErrMsgIdAndTxt("clipperplus:nargin", "function only takes the adjacency matrix as input.");
}

// number of graph nodes
const int nnodes = (useSparse) ? adjs.rows() : adj.rows();

const auto t1 = std::chrono::high_resolution_clock::now(); // timer

// run clipperplus_clique:
auto [clique, certificate] = clipperplus::find_clique(adj);
long clique_size = clique.size();

const auto t2 = std::chrono::high_resolution_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
const double elapsed = static_cast<double>(duration.count()) / 1e6;

// add 1 to clique indices (because matlab indexing starts at 1)
for (int i = 0; i < clique.size(); i++) {
    clique[i] += 1;
}

if (nlhs >= 1) {
  plhs[0] = mxCreateDoubleScalar(static_cast<double>(clique_size));
}

if (nlhs >= 2) {
    // Create a MATLAB matrix of integers
    mxArray* outputclique = mxCreateNumericMatrix(1, clique_size, mxINT32_CLASS, mxREAL);
    int* outputData = (int*)mxGetData(outputclique); // Cast to int pointer
    std::copy(clique.begin(), clique.end(), outputData);

    // Return the MATLAB matrix
    plhs[1] = outputclique;
}

if (nlhs >= 3) {
  plhs[2] = mxCreateDoubleScalar(static_cast<double>(certificate));
}

if (nlhs >= 4) {
  plhs[3] = mxCreateDoubleScalar(elapsed);
}

if (nlhs > 4) {
  mexErrMsgIdAndTxt("clipperplus:nargout", "Only up to 4 output args supported (clq_size, clq, cert, t)");
}


} //mexFunction
