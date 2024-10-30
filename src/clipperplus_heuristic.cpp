#include <set>

#include "clipperplus/utils.h"
#include "clipperplus/clipperplus_heuristic.h"


namespace clipperplus
{

Weight weighted_clique_size(const Graph &graph, std::vector<Node> &clique)
{
    std::vector<Weight> degrees = graph.induced(clique).degrees();
    Weight min_deg = *std::min_element(degrees.begin(), degrees.end());
    return min_deg;
}


std::vector<Node> find_heuristic_clique(
    const clipperplus::Graph &graph,
    int upper_bound,
    int lower_bound
) {
    auto kcores = graph.get_core_numbers();
    auto ordered_nodes = graph.get_core_ordering();
    auto max_core_number = graph.max_core_number();

    std::vector<Node> max_clique;
    Weight max_clique_size = lower_bound;   
    
    for(auto it = ordered_nodes.rbegin(); it != ordered_nodes.rend(); ++it) {
        auto v = *it;

        if(kcores[v] < max_clique_size) {
            continue;
        }

        std::vector<Node> S;
        for(auto [u, _] : graph.neighbors(v)) {
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

        Weight c_size = weighted_clique_size(graph, C);

        if(c_size > max_clique_size) {
            max_clique_size = c_size;
            max_clique = std::move(C);
        }

        if(max_clique_size == upper_bound) {
            break;
        }
    }

    return max_clique;
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
