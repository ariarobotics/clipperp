#pragma once

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include "clipperplus/clique_corenumber.h"
#include "clipperplus/clique_optimization.h"
#include "clipperplus/clipperplus_heuristic.h"
#include "clipperplus/utils.h"


namespace clipperplus 
{

enum class CERTIFICATE 
{
    NONE,
    HEURISTIC,
    CORE_BOUND,
    CHROMATIC_BOUND,
    CHORMATIC_WELSH
};

std::pair<std::vector<Node>, CERTIFICATE> find_clique(const Graph &graph);

unsigned long clipperplus_clique(const Eigen::MatrixXd& adj,
                       long& clique_size,
                       std::vector<int>& clique,
                       int& certificate);


} 
