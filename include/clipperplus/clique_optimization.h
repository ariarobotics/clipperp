// for CLIPPER (based on Belachew) optimization-based clique algorithm
#pragma once

#include <vector>
#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include "clipper/clipper.h"


namespace clipperplus {

int clique_optimization(const Eigen::MatrixXd& M, 
                        const Eigen::VectorXd& u0,
                        unsigned long& clique_size,
                        std::vector<long>& clique);


} 
