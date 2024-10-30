#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <Eigen/Dense>

#include "clipperplus/clipperplus_graph.h"

namespace clipperplus 
{

std::vector<Node> find_heuristic_clique(
    const clipperplus::Graph &graph,
    int upper_bound = -1,
    int lower_bound = 0
);


int estimate_chormatic_number_welsh_powell(const Graph &graph);

Weight weighted_clique_size(const Graph &graph, std::vector<Node> &clique);

} 
