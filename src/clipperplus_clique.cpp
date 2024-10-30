/*
computes a maximal clique in graph, and certifies if it's maximum clique

Author: kaveh fathian (kavehfathian@gmail.com)
 */

#include "clipperplus/clipperplus_heuristic.h"
#include "clipperplus/clipperplus_clique.h"

namespace clipperplus 
{


std::pair<std::vector<Node>, CERTIFICATE> find_clique(const Graph &graph)
{
    int n = graph.size();

    // auto chromatic_welsh = estimate_chormatic_number_welsh_powell(graph);
    auto k_core_bound = graph.max_core_number();
    std::vector<Weight> core_number = graph.get_core_numbers();


    auto heuristic_clique = find_heuristic_clique(graph); 
    auto heuristic_clique_size = weighted_clique_size(graph, heuristic_clique);

    if(heuristic_clique_size == k_core_bound) {
        return {heuristic_clique, CERTIFICATE::HEURISTIC};
    }

    std::vector<Node> keep, keep_pos(n, -1);
    for(Node i = 0, j = 0; i < n; ++i) {
        if(core_number[i] > heuristic_clique_size) {
            keep.push_back(i);
            keep_pos[i] = j++;
        }
    }

    Graph graph_pruned = graph.induced(keep);
    Eigen::MatrixXd M_pruned = graph_pruned.get_adj_matrix();
    M_pruned.diagonal().setOnes();

    Eigen::VectorXd u0 = Eigen::VectorXd::Ones(keep.size());

    for(auto v : heuristic_clique) {
        assert(keep_pos[v] >= 0);
        u0(keep_pos[v]) = 0;
    }
    u0.normalize();

    M_pruned = M_pruned / M_pruned.maxCoeff();
    auto clique_optim_pruned = clipperplus::clique_optimization(M_pruned, u0, Params());
    std::vector<Node> optimal_clique;

    Weight clique_optim_pruned_size = weighted_clique_size(graph_pruned, clique_optim_pruned);

    if(clique_optim_pruned_size < heuristic_clique_size) {
        optimal_clique = heuristic_clique;
    } else {
        for(auto v : clique_optim_pruned) {
            assert(v >= 0 && v < keep.size());
            optimal_clique.push_back(keep[v]);
        }
    }


    auto certificate = CERTIFICATE::NONE;
    if(optimal_clique.size() == k_core_bound) {
        certificate = CERTIFICATE::CORE_BOUND;
    }

    return {optimal_clique, certificate};
}

} 