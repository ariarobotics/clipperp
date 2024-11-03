
#include "clipperplus/clipperplus_graph.h"
#include <iostream>
#include <numeric>
#include <queue>
#include <unordered_map>

namespace clipperplus
{


Graph::Graph(Eigen::MatrixXd adj) : adj_matrix(std::move(adj)), adj_list(adj_matrix.rows())
{
    int nnodes = adj_matrix.rows();
    for(int i = 0; i < nnodes; ++i) {
        for(int j = 0; j < nnodes; ++j) {
            if(adj_matrix(i, j) != 0) {
                adj_list[i].push_back(j);
            }
        }
    }
}




const std::vector<Node> &Graph::neighbors(Node v) const 
{
    assert(v < adj_list.size());
    return adj_list[v];
}


int Graph::degree(Node v) const 
{
    assert(v < adj_list.size());
    return adj_list[v].size();
}


std::vector<int> Graph::degrees() const
 {
    int nnodes = adj_list.size();
    std::vector<int> degrees(adj_list.size());
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

        adj_matrix(i, g.adj_list[i]).setOnes();
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

        for(auto u : neighbors(v)) {
            if(keep[u] >= 0) {
                g.adj_list[keep[v]].push_back((Node)keep[u]);
            }
        }
    }

    return g;
}


int Graph::size() const 
{
    return adj_list.size();
}


int Graph::max_core_number() const
{
    if(kcore.empty()) {
        calculate_kcores();
    }
    return kcore[kcore_ordering.back()];
}


const std::vector<int> &Graph::get_core_numbers() const
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
    int n = size();
    auto degree = degrees();
    auto max_degree = *std::max_element(degree.begin(), degree.end());

    // prepare for bucket sort by degree
    // pos[i] is the position of v_i on the sorted array
    // If you consider `pos` as an permutation `order` essentially is pos^-1
    std::vector<int> pos(n, 0);
    kcore_ordering.resize(n);

    std::vector<int> bin(max_degree + 1, 0);
    for(auto d : degree) {
        ++bin[d];
    }

    int start = 0;
    for(int d = 0; d < max_degree + 1; ++d) {
        auto num = bin[d];
        bin[d] = start;
        start += num;
    }

    // bucket sort:
    for(int v = 0; v < n; ++v) {
        pos[v] = bin[degree[v]];
        kcore_ordering[pos[v]] = v;
        ++bin[degree[v]];
    }

    for(int d = max_degree; d > 0; --d) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    // iteratively remove edges from v with lowest degree
    for(auto v : kcore_ordering) {
        for(auto u : neighbors(v)) { // remove e = (v, u)
            if(degree[v] >= degree[u]) {
                continue;
            }

            // update sorted array: pos, order bin
            // find first element in sorted array with d[w] = d[u]
            auto pos_w = bin[degree[u]];
            auto w = kcore_ordering[pos_w];

            // swap their pose and order
            if(w != u) {
                kcore_ordering[pos[u]] = w;
                kcore_ordering[pos[w]] = u;
                std::swap(pos[u], pos[w]);
            }


            ++bin[degree[u]];
            --degree[u];
        }
    }

    kcore = std::move(degree);
}


template <typename Key, typename Value>
class DecreaseKeyPriorityQueue {
private:
    struct HeapElement {
        Key key;
        Value value;

        bool operator>(const HeapElement &other) const {
            return key > other.key;
        }
    };

    std::vector<HeapElement> heap;
    std::unordered_map<Value, size_t> position_map; // Maps values to their positions in the heap
    std::function<bool(const Key&, const Key&)> comparator;

    void heapifyUp(size_t index) {
        while (index > 0) {
            size_t parent = (index - 1) / 2;
            if (comparator(heap[index].key, heap[parent].key)) {
                swapElements(index, parent);
                index = parent;
            } else {
                break;
            }
        }
    }

    void heapifyDown(size_t index) {
        size_t size = heap.size();
        while (true) {
            size_t left = 2 * index + 1;
            size_t right = 2 * index + 2;
            size_t smallest = index;

            if (left < size && comparator(heap[left].key, heap[smallest].key)) {
                smallest = left;
            }
            if (right < size && comparator(heap[right].key, heap[smallest].key)) {
                smallest = right;
            }
            if (smallest != index) {
                swapElements(index, smallest);
                index = smallest;
            } else {
                break;
            }
        }
    }

    void swapElements(size_t i, size_t j) {
        std::swap(heap[i], heap[j]);
        position_map[heap[i].value] = i;
        position_map[heap[j].value] = j;
    }

public:
    DecreaseKeyPriorityQueue(std::function<bool(const Key&, const Key&)> comp = std::less<Key>())
        : comparator(comp) {}

    bool empty() const {
        return heap.empty();
    }

    void push(const Key &key, const Value &value) {
        if (position_map.find(value) != position_map.end()) {
            throw std::invalid_argument("Value already exists in the priority queue");
        }
        heap.push_back({key, value});
        size_t index = heap.size() - 1;
        position_map[value] = index;
        heapifyUp(index);
    }

    void pop() {
        if (heap.empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        position_map.erase(heap[0].value);
        heap[0] = heap.back();
        heap.pop_back();
        if (!heap.empty()) {
            position_map[heap[0].value] = 0;
            heapifyDown(0);
        }
    }

    const HeapElement& top() const {
        if (heap.empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return heap[0];
    }

    void decreaseKey(const Value &value, const Key &newKey) {
        auto it = position_map.find(value);
        if (it == position_map.end()) {
            throw std::invalid_argument("Value not found in the priority queue");
        }
        size_t index = it->second;
        if (comparator(heap[index].key, newKey)) {
            throw std::invalid_argument("New key is not smaller than the current key");
        }
        heap[index].key = newKey;
        heapifyUp(index);
    }
};

const std::vector<int> &Graph::get_triangle_core() const
{
    if (triangle_core.empty()) {
        calculate_triangle_core();
    }
    return triangle_core;
}

const std::vector<Node> &Graph::get_triangle_core_ordering() const
{
    if (triangle_core.empty()) {
        calculate_triangle_core();
    }
    return triangle_core_ordering;
}

void Graph::calculate_triangle_core() const
{
    int n = size();
    triangle_core.resize(n);
    for (Node u = 0; u < n; ++u) {
        for (auto v : neighbors(u)) {
            if (u >= v) {
                continue;
            }
            
            for (auto w : neighbors(v)) {
                if (v >= w || !is_edge(u, w)) {
                    continue;
                }

                for (auto n : {u, v, w}) {
                    ++triangle_core[n];
                }
            }
        }
    }

    using HeapElement = std::tuple<int, Node>;
    DecreaseKeyPriorityQueue<int, Node> queue(
        [](int a, int b) { return a < b; }
    );


    for (auto node = 0; node < n; ++node) {
        queue.push(triangle_core[node], node);
    }

    int k = 0;
    std::vector<bool> removed(n, false);
    while(!queue.empty()) {
        auto [tri, node] = queue.top();
        queue.pop();
        if (removed[node]) {
            continue; 
        } 

        removed[node] = true;
            k = std::max(k, triangle_core[node]);
            triangle_core[node] = k;
            for (auto u : neighbors(node)) {
                if (removed[u]) {
                    continue;
                }
                for (auto v : neighbors(u)) {
                    if (removed[v] || v <= u) {
                        continue;
                    }
                    if (is_edge(node, v) && !removed[v]) {
                        --triangle_core[v];
                        --triangle_core[u];

                        queue.decreaseKey(v, triangle_core[v]);
                        queue.decreaseKey(u, triangle_core[u]);
                    }
                }
            }
        }

        triangle_core_ordering.clear();
        std::iota(triangle_core_ordering.begin(), triangle_core_ordering.end(), 0);
        std::sort(triangle_core_ordering.begin(), triangle_core_ordering.end(), 
            [this](Node a, Node b) { 
                return triangle_core[a] < triangle_core[b]; 
        });
    }
}





