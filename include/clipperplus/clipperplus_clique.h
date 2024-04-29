#pragma once

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

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
    CHROMATIC_BOUND
};

std::pair<std::vector<Node>, CERTIFICATE> find_clique(const Graph &graph);

} 
