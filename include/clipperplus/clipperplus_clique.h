#pragma once

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include "clipperplus/clique_corenumber.h"
#include "clipperplus/clique_optimization.h"
#include "clipperplus/utils.h"


namespace clipperplus {

unsigned long clipperplus_clique(const Eigen::MatrixXd& adj,
                       long& clique_size,
                       std::vector<int>& clique,
                       Certificate& certificate);

} 
