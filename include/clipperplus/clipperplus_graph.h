#pragma once

#include <Eigen/Dense>
#include <memory.h>
#include <vector>



namespace clipperplus
{

using Node = unsigned int;
using Edge = std::pair<Node, Node>;
using Neighborlist = std::vector<Node>;


class Graph
{
public:
    Graph(Eigen::MatrixXd adj_matrix);
    // static Graph from_list(const std::vector<Neighborlist> &adj_list);

    int size() const;
    int degree(Node v) const;
    std::vector<int> degrees() const;

    const std::vector<Node> &neighbors(Node v) const;

    bool is_edge(Node u, Node v) const;

    void merge(const Graph &g);

    int max_core_number() const;
    const std::vector<int> &get_core_numbers() const;
    const std::vector<Node> &get_core_ordering() const;
    
    const Eigen::MatrixXd &get_adj_matrix() const;

private:
    void calculate_kcores() const;

private:
    Eigen::MatrixXd adj_matrix;
    std::vector<Neighborlist> adj_list;

    mutable std::vector<Node> kcore_ordering;
    mutable std::vector<int> kcore;
};


}
