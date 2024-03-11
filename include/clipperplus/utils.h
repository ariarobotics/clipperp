#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>

namespace clipperplus{

// find the index of an element in an std::vector of integers 
// that is equal to a given value by iterating through the 
// vector and checking each element. 
int find_index(const std::vector<int>& vec, int val);

void adjmat_to_adjlist(const Eigen::MatrixXd& adj,
                       const int& nnodes,
                       int* ei,
                       int* ej);


} 