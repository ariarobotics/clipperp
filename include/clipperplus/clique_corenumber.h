/**
Algorithm for computing core number of nodes in a graph, 
and finding a heuristic maximum clique in a computationally efficient manner.
Code is based on Parallel Maximum Clique Algorithm by Ryan A. Rossi (http://ryanrossi.com/pmc)

Author: kaveh fathian (kavehfathian@gmail.com)
 */
#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <Eigen/Dense>

#include "pmc/pmc.h"
#include "pmc/pmc_neigh_coloring.h"


namespace clipperplus {

// find heuristic clique using core numbers
int find_heuristic_clique(const Eigen::MatrixXd& adj, 
                    const std::vector<int>& core_numbers,
                    std::vector<int>& clique);


// get as input an adjacency matrix
unsigned long clique_corenumber(const Eigen::MatrixXd& adj,
                      std::vector<int>& clique,
                      std::vector<long>& core_numbers,
                      long& core_bound,
                      std::vector<int>& node_colors,
                      int& chromatic_bound); 


} 