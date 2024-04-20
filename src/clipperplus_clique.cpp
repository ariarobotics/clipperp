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

    auto heuristic_clique = find_heuristic_clique(graph);    
    if(heuristic_clique.size() == graph.max_core_number() + 1) {
        return {heuristic_clique, CERTIFICATE::HEURISTIC};
    }

    auto max_core_number = graph.max_core_number();
    std::vector<int> core_number = graph.get_core_numbers();

    std::vector<int> keep, keep_pos(n, -1);
    for(Node i = 0, j = 0; i < n; ++i) {
        if(core_number[i] >= heuristic_clique.size()) {
            keep.push_back(i);
            keep_pos[i] = j++;
        }
    }

    Eigen::MatrixXd M_pruned = graph.get_adj_matrix()(keep, keep);
    M_pruned.diagonal().setOnes();

    Eigen::VectorXd u0 = Eigen::VectorXd::Ones(keep.size());

    for(auto v : heuristic_clique) {
        assert(keep_pos[v] >= 0);
        u0(keep_pos[v]) = 0;
    }
    u0.normalize();

    unsigned long clique_size_optim;
    std::vector<long> clique_optim_pruned;
    clipperplus::clique_optimization(M_pruned, u0, clique_size_optim, clique_optim_pruned);

    std::vector<Node> optimal_clique;
    if(clique_optim_pruned.size() < heuristic_clique.size()) {
        optimal_clique = heuristic_clique;
    } else {
        for(auto v : clique_optim_pruned) {
            assert(v >= 0 && v < keep.size());
            optimal_clique.push_back(keep[v]);
        }
    }


    auto chromatic_number = estimate_chromatic_number(graph);
    auto chromatic_welsh = estimate_chormatic_number_welsh_powell(graph);
    std::cout << chromatic_number << " vs " << chromatic_welsh << std::endl;

    auto certificate = CERTIFICATE::NONE;
    if(optimal_clique.size() == max_core_number + 1) {
        certificate = CERTIFICATE::CORE_BOUND;
    } else if(optimal_clique.size() == chromatic_number) {
        certificate = CERTIFICATE::CHROMATIC_BOUND;
    } else if(optimal_clique.size() == chromatic_welsh) {
        certificate = CERTIFICATE::CHORMATIC_WELSH;
    }

    return {optimal_clique, certificate};
}

} 