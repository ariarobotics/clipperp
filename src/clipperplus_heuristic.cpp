#include "clipperplus/utils.h"
#include "clipperplus/clipperplus_heuristic.h"

namespace clipperplus
{


std::vector<Node> find_heuristic_clique(
    const clipperplus::Graph &graph,
    int upper_bound,
    int lower_bound
) {
    auto kcores = graph.get_core_numbers();
    auto ordered_nodes = graph.get_core_ordering();
    auto max_core_number = graph.max_core_number();

    std::vector<Node> max_clique;
    int max_clique_size = lower_bound;   
    
    for(auto it = ordered_nodes.rbegin(); it != ordered_nodes.rend(); ++it) {
        auto v = *it;

        if(kcores[v] < max_clique_size) {
            continue;
        }

        std::vector<Node> S;
        for(auto u : graph.neighbors(v)) {
            if(kcores[u] >= max_clique_size) {
                S.push_back(u);
            }
        }
        assert(S.size() >= max_clique_size);

        std::stable_sort(S.begin(), S.end(), [&](Node a, Node b) {
            return kcores[a] > kcores[b];
        });

        std::vector<Node> C = { v };
        for(auto u : S) {
            auto connected = std::all_of(C.begin(), C.end(), [&](Node v){
                return graph.is_edge(u, v);
            });

            if(connected) {
                C.push_back(u);
            }
        }

        if(C.size() > max_clique_size) {
            max_clique_size = C.size();
            max_clique = std::move(C);
        }

        if(max_clique_size == upper_bound) {
            break;
        }
    }

    return max_clique;
}


int estimate_chromatic_number(const Graph &graph)
{
    std::vector<int> node_color(graph.size(), 0);
    auto node_order = graph.get_core_ordering();

    auto core_numbers = graph.get_core_numbers();
    
    for(auto it = node_order.begin(); it != node_order.end(); ++it) {
        auto v = *it;
        std::set<int> neighbor_colors;

        for(auto u : graph.neighbors(v)) {
            neighbor_colors.insert(node_color[u]);
        }

        int color = 1;
        while(neighbor_colors.count(color) > 0) {
            ++color;
        }

        node_color[v] = color;
    }

    return *std::max_element(node_color.begin(), node_color.end());
}


int estimate_chormatic_number_welsh_powell(const Graph &graph)
{
    std::vector<int> node_color(graph.size(), -1);
    auto degrees = graph.degrees();

    std::vector<Node> node_order(graph.size(), 0);
    std::iota(node_order.begin(), node_order.end(), 0);
    std::sort(node_order.begin(), node_order.end(), [&degrees](Node a, Node b) { return degrees[a] > degrees[b]; });

    int color = 0;
    for(auto v : node_order) {
        std::set<int> colored; // colored using `color = i`
        if(node_color[v] != -1) {
            continue;
        }

        node_color[v] = color;
        colored.insert(v);

        for(auto u : node_order) {
            if(u == v || node_color[u] != -1) {
                continue;
            }

            auto conflicts = std::any_of(colored.begin(), colored.end(), [&](Node v){
                return graph.is_edge(u, v);
            });

            if(!conflicts) {
                node_color[u] = color;
                colored.insert(u);
            }
        }

        ++color;
    }

    return color;
}

}
