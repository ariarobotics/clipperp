#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>


enum Certificate {NONE, PRUNING, KCORE, CHROMATIC};

std::string cert_to_string(Certificate cert);

namespace clipperplus{

// find the index of an element in an std::vector of integers 
// that is equal to a given value by iterating through the 
// vector and checking each element. 
int find_index(const std::vector<int>& vec, int val);

void adjmat_to_adjlist(const Eigen::MatrixXd& adj,
                       const int& nnodes,
                       std::vector<int>& ei,
                       std::vector<int>& ej);

} 
