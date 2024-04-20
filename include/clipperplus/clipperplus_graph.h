#pragma once

#include <Eigen/Dense>
#include <memory.h>
#include <vector>


namespace clipperplus
{

using Node = int;
using Edge = std::pair<Node, Node>;
using Neighborlist = std::vector<Node>;


class Graph
{
public:
    Graph() = default;
    Graph(Eigen::MatrixXd adj_matrix);
    // static Graph from_list(const std::vector<Neighborlist> &adj_list);

    int size() const;
    int degree(Node v) const;
    std::vector<int> degrees() const;

    const std::vector<Node> &neighbors(Node v) const;

    inline bool is_edge(Node u, Node v) const
    {
        return adj_matrix(u, v) != 0;
    }

    void merge(const Graph &g);
    Graph induced(const std::vector<Node> &nodes) const;

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
