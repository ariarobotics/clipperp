/**
 * @brief MATLAB/MEX binding 
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
const long nnodes = useSparse ? adjs.rows() : adj.rows();

const auto t1 = std::chrono::high_resolution_clock::now(); // timer

// run clipperplus_clique:
std::vector<int> clique; // initialize index vector of clique
clique.reserve(nnodes); // allocate memory
std::vector<long int> core_numbers(nnodes,0); // initialize vector of core numbers
long int core_bound = 0; // initialize max clique upper bound based on max kcore
int chromatic_bound = 0;; // initialize chromatic number upper bound
std::vector<int> node_colors(nnodes, 0); // initialize graph node coloring

clipperplus::Graph g(adj);
for(int i = 0; i < nnodes; i++) {
    core_numbers[i] = g.get_core_numbers()[i];
}

core_bound = g.max_core_number() + 1;
chromatic_bound = clipperplus::estimate_chormatic_number_welsh_powell(g);
clique = clipperplus::find_heuristic_clique(g, core_bound < chromatic_bound ? core_bound : chromatic_bound);

// clipperplus::clique_corenumber(adj, clique, core_numbers, 
//                               core_bound, node_colors, chromatic_bound);

long clique_size = clique.size();

const auto t2 = std::chrono::high_resolution_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
const double elapsed = static_cast<double>(duration.count()) / 1e6;

// add 1 to clique indices (because matlab indexing starts at 1)
for (int& index : clique) {
    index += 1;
}

if (nlhs >= 1) {
  plhs[0] = mxCreateDoubleScalar(static_cast<double>(clique_size));
}

if (nlhs >= 2) {
    // Create a MATLAB matrix of integers
    mxArray* outputclique = mxCreateNumericMatrix(1, clique_size, mxINT32_CLASS, mxREAL);
    auto outputData = static_cast<int*>(mxGetData(outputclique)); // Cast to int pointer 
    std::copy(clique.begin(), clique.end(), outputData);

    // Return the MATLAB matrix
    plhs[1] = outputclique;
}

if (nlhs >= 3) {
  plhs[2] = mxCreateDoubleScalar(elapsed);
}

if (nlhs > 3) {
  mexErrMsgIdAndTxt("nargout", "Only up to 3 output args supported (clq_size, clq, t)");
}


} //mexFunction
