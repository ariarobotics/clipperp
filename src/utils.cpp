
#include "clipperplus/utils.h"


std::string cert_to_string(Certificate cert){
    switch(cert){
        case NONE: return "NONE";
        case PRUNING: return "PRUNING";
        case KCORE: return "KCORE";
        case CHROMATIC: return "CHROMATIC";
        default: return "Unknown Certificate";
    }
}

namespace clipperplus {

/** find the index of an element in an std::vector of integers 
* that is equal to a given value by iterating through the 
* vector and checking each element. */
int find_index(const std::vector<int>& vec, int val) {
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] == val) {
            return i; // Return the index of the first occurrence of 'val'
        }
    }
    return -1; // Return -1 if 'val' is not found in the vector
}

// convert adjacency matrix to adjacency list
void adjmat_to_adjlist(const Eigen::MatrixXd& adj,
                       const int& nnodes, 
                       std::vector<int>& ei,
                       std::vector<int>& ej) {

    if (adj.rows() != adj.cols()) {
        throw std::invalid_argument("adjacency matrix must be square");
    }
    ei.clear();
    ej.clear();
    for (int i=0; i<nnodes-1; i++) {
        for (int j=i+1; j<nnodes; j++) {
        if (adj(i,j)==1) {
                ei.push_back(i);
                ej.push_back(j);
            }
        }
    }
} 

} 
