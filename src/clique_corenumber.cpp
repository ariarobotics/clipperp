/*
compute core number of graph nodes, and use them to find a miximal clique.
Code inspired by PMC Algorithm by Ryan A. Rossi (http://ryanrossi.com/pmc)
Author: kaveh fathian (fathian@ariarobotics.com)
*/

#include "clipperplus/clique_corenumber.h"

namespace clipperplus {

// find heuristic clique using core numbers
int find_heuristic_clique(const Eigen::MatrixXd& adj, 
                          const std::vector<long>& core_numbers,
                          std::vector<int>& clique) {
    const long nnodes = adj.cols(); // number of nodes

    // if graph has no edges, return trivial clique
    if (adj.sum() == 0) { 
        clique.push_back(0);
        return 1;
        }

    long max_val = 0;

    // order of nodes in decreasing core number
    std::vector<int> node_order(core_numbers.size());
    std::iota(node_order.begin(), node_order.end(), 0); // initialize node_order with sequential integers from 0 to core_numbers.size()-1
    std::stable_sort(node_order.begin(), node_order.end(), [&core_numbers](int a, int b) {
        return core_numbers[a] > core_numbers[b]; });

    for (int i : node_order) {

        if (core_numbers[i] > max_val) { 
        // subset of neighbors of node 'i' with core numbers >= max_val    
        std::vector<int> nb_idx; // index of nonzero elements on row 'i'
        nb_idx.reserve(nnodes); // reserve memory
        for (int j=0; j<nnodes; ++j) { 
            if (adj(i,j) != 0) {
                nb_idx.push_back(j);
            }
        }

        std::vector<int> select_idx; // index of core_numbers(nb_idx) >= max_val
        select_idx.reserve(nnodes); // reserve memory
        for (int k=0; k<nb_idx.size(); ++k) {
            if (core_numbers[nb_idx[k]] >= max_val) {
                select_idx.push_back(k);
            }
        }
        std::vector<int> subset_idx;
        subset_idx.reserve(nnodes); // reserve memory
        for (int idx : select_idx) {
            if (idx >= 0 && idx < nb_idx.size()) {
                subset_idx.push_back(nb_idx[idx]);
            }
        }

        if (subset_idx.empty()) {continue;}

        std::vector<int> clq;
        clq.reserve(nnodes); // reserve memory
        clq.push_back(i); // initialize clique with node 'i'

        // order selected subset of neighbors in decresing core number
        std::vector<long> core_numbers_subset_idx;
        core_numbers_subset_idx.reserve(nnodes); // reserve memory
        for (int idx : subset_idx) {
            core_numbers_subset_idx.push_back(core_numbers[idx]);        
        }
        std::vector<int> order_innerloop(core_numbers_subset_idx.size());
        std::iota(order_innerloop.begin(), order_innerloop.end(), 0); // initialize with sequential integers from 0 to size()-1
        std::stable_sort(order_innerloop.begin(), order_innerloop.end(), [&](int a, int b) {
            return core_numbers_subset_idx[a] > core_numbers_subset_idx[b];
        });

        for (int idx : order_innerloop) {
            int j = subset_idx[idx]; // node index
            // if node 'j' union with 'clq' is a clique, then add 'j' to 'clq'        
            bool allElementsTrue = true;
            for (int idx_clq : clq) {
                if (adj(j,idx_clq) == 0) {
                    allElementsTrue = false;
                    break; // No need to continue checking
                }
            }
            if (allElementsTrue) {clq.push_back(j);} //add 'j' to 'clq'
        }

        if (clq.size() > max_val) { // store better clique
            clique = clq;
            max_val = clq.size();
        }    
        } 
    } 

    return 1;
}

unsigned long clique_corenumber(const Eigen::MatrixXd& adj,
                      std::vector<int>& clique,
                      std::vector<long>& core_numbers,
                      long& core_bound,
                      std::vector<int>& node_colors,
                      int& chromatic_bound) {
    
    const long nnodes = adj.rows(); // number of graph nodes

    // create a PMC graph from the adjacency matrix
    std::vector<long> edges;
    std::vector<long> vertices;
    vertices.push_back(0);

    size_t total_num_edges = 0;
    for (size_t i=0; i<nnodes; i++) {
        for (size_t j=0; j<nnodes; j++) {
            if (adj(i,j)==1) {
                edges.push_back(j);
                total_num_edges++;
            }
        }
        vertices.push_back(total_num_edges);
    }

    // construct PMC graph
    pmc::pmc_graph G(vertices, edges);    
    pmc::input in;

    // compute k-cores 
    G.compute_cores(); // compute k-cores
    in.ub = G.get_max_core() + 1; // max clique upper bound as max k-core + 1
    core_bound = in.ub; // output
    std::vector<long> kcores = G.kcore; // k-core values for each node

    // core number ouput
    for(int i = 0; i < nnodes; i++) {
        core_numbers[i] = kcores[i] - 1;
    }
    
    #ifdef DEBUG
        std::cout << "kcores.size(): " << kcores.size() << std::endl;
        std::cout << "core_numbers.size(): " << core_numbers.size() << std::endl;
        std::cout << "max clique kcore bound: " << core_bound << std::endl;        
        std::cout << "core numbers: ";
        for (int elm : core_numbers){std::cout << elm << " ";}
        std::cout << std::endl;
    #endif
  
    // find a maximal clique using core numbers
    #ifdef DEBUG
        std::cout << "running find_heuristic_clique..."  << std::endl;
    #endif
    find_heuristic_clique(adj, core_numbers, clique);
    const unsigned long clique_size = clique.size();
    #ifdef DEBUG
        std::cout << "heuristic found clique of size: " << clique_size << std::endl;        
        std::cout << "heuristic clique: ";
        for(int elm : clique){std::cout << elm << " ";} std::cout << std::endl;
        if (in.lb == in.ub) {std::cout << "heuristic found max clique." << std::endl;}
    #endif

    // graph chromatic number (estimate)
    std::vector<int> nodeOrder(nnodes); // vector of node indices sorted by kcores in descending order
    std::iota(nodeOrder.begin(), nodeOrder.end(), 0); // set elements with sequential integers from 0 to n-1
    std::sort(nodeOrder.begin(), nodeOrder.end(), [&](int a, int b) {
        return kcores[a] > kcores[b];
    });

    for (int i = 0; i < nnodes; i++) {
        int node = nodeOrder[i];
        std::set<int> neighborColors;

        // find the colors used by neighboring nodes
        for (int j = 0; j < nnodes; j++) {
            if (adj(node,j)==1) {
                neighborColors.insert(node_colors[j]);
            }
        }

        // find the minimum available color
        int color = 1;
        while (neighborColors.count(color) > 0) {
            color++;
        }

        node_colors[node] = color;
    }

    chromatic_bound = *std::max_element(node_colors.begin(), node_colors.end());

    #ifdef DEBUG
        std::cout << "colors: ";
        for (int i : node_colors) {std::cout << i << " ";} std::cout << std::endl;
        std::cout << "max clique chromatic bound: " << chromatic_bound << std::endl;
    #endif

    // return output
    return clique_size;


} 

} 