
#include "clipperplus/clipperplus_graph.h"
#include <queue>
#include <iostream>

namespace clipperplus
{


Graph::Graph(Eigen::MatrixXd adj) : adj_matrix(std::move(adj)), adj_list(adj_matrix.rows())
{
    int nnodes = adj_matrix.rows();
    for(int i = 0; i < nnodes; ++i) {
        for(int j = 0; j < nnodes; ++j) {
            if(adj_matrix(i, j) != 0) {
                adj_list[i].push_back({j, adj_matrix(i, j)});
            }
        }
    }
}




const std::vector<Neighbor> &Graph::neighbors(Node v) const 
{
    assert(v < adj_list.size());
    return adj_list[v];
}


Weight Graph::degree(Node v) const 
{
    assert(v < adj_list.size());
    Weight sum = 0;
    for(const auto [v, w] : adj_list[v]) {
        sum += w;
    }
    return sum;
}


std::vector<Weight> Graph::degrees() const
 {
    int nnodes = adj_list.size();
    std::vector<Weight> degrees(adj_list.size());
    for(int i = 0; i < nnodes; ++i) {
        degrees[i] = degree(i);
    }
    return degrees;
}


void Graph::merge(const Graph &g)
{
    assert(adj_list.size() == g.adj_list.size());
    for(int i = 0; i < adj_list.size(); ++i) {
        adj_list[i].insert(
            adj_list[i].end(), 
            g.adj_list[i].begin(), g.adj_list[i].end()
        );

        for(auto [u, w] : g.adj_list[i]) {
            adj_matrix(i, u) = w;
        }
    }

    kcore.clear();
    kcore_ordering.clear();
}


Graph Graph::induced(const std::vector<Node> &nodes) const
{
    int n = size();

    auto g = Graph();
    g.adj_matrix = adj_matrix(nodes, nodes);


    std::vector<int> keep(size(), -1);
    for(Node i = 0; i < nodes.size(); ++i) {
        keep[nodes[i]] = i;
    }

    g.adj_list = std::vector<Neighborlist>(nodes.size());
    for(Node v = 0; v < n; ++v) {
        if(keep[v] < 0) {
            continue;
        }

        for(auto [u, w] : neighbors(v)) {
            if(keep[u] >= 0) {
                g.adj_list[keep[v]].push_back({keep[u], w});
            }
        }
    }

    return g;
}


int Graph::size() const 
{
    return adj_list.size();
}


Weight Graph::max_core_number() const
{
    if(kcore.empty()) {
        calculate_kcores();
    }
    return kcore[kcore_ordering.back()];
}


const std::vector<Weight> &Graph::get_core_numbers() const
{
    if(kcore.empty()) {
        calculate_kcores();
    }
    return kcore;
}


const std::vector<Node> &Graph::get_core_ordering() const 
{
    if(kcore.empty()) {
        calculate_kcores();
    }
    return kcore_ordering;
}


const Eigen::MatrixXd &Graph::get_adj_matrix() const
{
    return adj_matrix;
}


void Graph::calculate_kcores() const
{
    // only cause std::priority queue donsn't support decrese-key
    const int n = size();
    std::vector<bool> removed(n, false);

    auto degrees = this->degrees();
    kcore.resize(n);

    using HeapElement = std::tuple<Weight, Node>; 
    
    std::priority_queue<
        HeapElement, std::vector<HeapElement>, std::greater<HeapElement>
    > queue;

    for(auto node = 0; node < n; ++node) {
        queue.push({degrees[node], node});
    }

    Weight k = 0;
    while(!queue.empty()) {
        auto [_, node] = queue.top(); 
        queue.pop();
        if(removed[node]) { // decrece key instead of this
            continue;
        }
        removed[node] = true;

        k = std::max(k, degrees[node]);        
        kcore[node] = k;

        for(auto &[neighbor, weight] : neighbors(node)) {
            degrees[neighbor] -= weight;
            queue.push({degrees[neighbor], neighbor}); // decrece key instead
        }
    }

    kcore_ordering.clear();
    for(auto node = 0; node < n; ++node) {
        kcore_ordering.push_back(node);
    }
    std::sort(kcore_ordering.begin(), kcore_ordering.end(), 
        [this](Node a, Node b) { return kcore[a] < kcore[b]; }
    );
}

}


