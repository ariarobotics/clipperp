#pragma once

#include <Eigen/Dense>
#include <vector>

#include "pmc/pmc.h"
#include "pmc/pmc_neigh_coloring.h"

namespace pmc {

// get as input an adjacency matrix
long pmc_heuristic(const Eigen::MatrixXd& adj,
                  long& clique_size,
                  std::vector<long>& clique);    

}                       