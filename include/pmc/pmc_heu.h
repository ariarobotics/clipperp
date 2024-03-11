/**
 ============================================================================
 Name        : Parallel Maximum Clique (PMC) Library
 Author      : Ryan A. Rossi   (rrossi@purdue.edu)
 Description : A general high-performance parallel framework for computing
               maximum cliques. The library is designed to be fast for large
               sparse graphs.

 Copyright (C) 2012-2013, Ryan A. Rossi, All rights reserved.

 Please cite the following paper if used:
   Ryan A. Rossi, David F. Gleich, Assefaw H. Gebremedhin, Md. Mostofa
   Patwary, A Fast Parallel Maximum Clique Algorithm for Large Sparse Graphs
   and Temporal Strong Components, arXiv preprint 1302.6256, 2013.

 See http://ryanrossi.com/pmc for more information.
 ============================================================================
 */

#ifndef PMC_HEU_H_
#define PMC_HEU_H_

#include "pmc_headers.h"
#include "pmc_graph.h"
#include "pmc_utils.h"
#include "pmc_input.h"
#include "pmc_vertex.h"
#include <algorithm>

namespace pmc {

class pmc_heu {
    public:
        vector<long>* E;
        vector<long>* V;
        vector<long>* K;
        vector<long>* order;
        vector<long>* degree;
        double sec;
        long ub;
        string strat;

        int num_threads;

        pmc_heu(pmc_graph& G, 
                input const& params) 
            : K(G.get_kcores()),
              order(G.get_kcore_ordering()), 
              ub(params.ub),
              strat(params.heu_strat),
              num_threads(params.threads) {
          initialize(); 
        }

        pmc_heu(pmc_graph& G, 
                int tmp_ub) 
            : K(G.get_kcores()),
              order(G.get_kcore_ordering()),
              ub(tmp_ub),
              strat("kcore") {
          initialize();
        }

        inline void initialize() {
            sec = get_time();
            srand(static_cast<unsigned int>(time(nullptr)));
        };

        int strategy(vector<int>& P);
        void set_strategy(string const& s) { strat = s; }
        long compute_heuristic(long v);

        static bool desc_heur(Vertex v,  Vertex u) {
            return (v.get_bound() > u.get_bound());
        }

        static bool incr_heur(Vertex v,  Vertex u) {
            return (v.get_bound() < u.get_bound());
        }

        int search(pmc_graph& graph, vector<long>& C_max);
        int search_cores(const pmc_graph& graph, vector<long>& C_max, int lb);
        int search_bounds(pmc_graph& graph, vector<long>& C_max);

        void branch(vector<Vertex>& P, int sz,
                int& mc, vector<long>& C, vector<short>& ind);

        void print_info(const vector<long>& C_max) const;
}; 

}; 
#endif
