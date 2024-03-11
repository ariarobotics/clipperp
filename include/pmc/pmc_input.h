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

#ifndef PMC_INPUT_H_
#define PMC_INPUT_H_

#include "pmc_headers.h"
#include "pmc_utils.h"

using namespace std;

namespace pmc {
class input {
    public:
        // instance variables
        int algorithm = 0;
        int threads = omp_get_max_threads();
        int experiment = 0;
        int lb = 0;
        long ub = 0;
        int param_ub = 0;
        int adj_limit = 20000;
        double time_limit = 60 * 60; // max time to search
        double remove_time = 4.0; // time to wait before reducing graph
        bool graph_stats = false;
        bool verbose = false;
        bool help = false;
        bool MCE = false;
        bool decreasing_order = false;
        string heu_strat = "kcore";
        string format = "mtx";
        string graph = "data/sample.mtx";
        string output = "";
        string edge_sorter = "";
        string vertex_search_order = "deg";

        input() {
            if (threads <= 0) threads = 1;
        }

        input(int argc, char *argv[]) {
            int opt;
            while ((opt=getopt(argc,argv,"i:t:f:u:l:o:e:a:r:w:h:k:dgsv")) != EOF) {
                switch (opt) {
                    case 'a':
                        algorithm = atoi(optarg);
                        if (algorithm > 9) MCE = true;
                        break;
                    case 't':
                        threads = atoi(optarg);
                        break;
                    case 'f':
                        graph = optarg;
                        break;
                    case 's':
                        graph_stats = true;
                        break;
                    case 'u':
                        param_ub = atoi(optarg); // find k-clique fast
                        lb = 2;                  // skip heuristic
                        break;
                    case 'k':
                        param_ub = atoi(optarg);
                        lb = param_ub-1;
                        break;
                    case 'l':
                        lb = atoi(optarg);
                        break;
                    case 'h':
                        heu_strat = optarg;
                        break;
                    case 'v':
                        verbose = true;
                        break;
                    case 'w':
                        time_limit = atof(optarg) * 60;  // convert minutes to seconds
                        break;
                    case 'r':
                        remove_time = atof(optarg);
                        break;
                    case 'e':
                        edge_sorter = optarg;
                        break;
                    case 'o':
                        vertex_search_order = optarg;
                        break;
                    case 'd':
                        // direction of which vertices are ordered
                        decreasing_order = true;
                        break;
                    case '?':
                        usage(argv[0]);
                        break;
                    default:
                        usage(argv[0]);
                        break;
                }
            }

            // both off, use default alg
            if (heu_strat == "0" && algorithm == -1)
                algorithm = 0;

            if (threads <= 0) threads = 1;

            if (!fexists(graph.c_str())) {
                usage(argv[0]);
                exit(-1);
            }

            FILE* fin = fopen(graph.c_str(), "r+t");
            if (fin == nullptr) {
                usage(argv[0]);
                exit(-1);
            }
            fclose(fin);

            cout << "\n\nFile Name ------------------------ " << graph.c_str() << endl;
            if (!fexists(graph.c_str()) ) {
                cout << "File not found!" << endl;
                return;
            }
            cout << "workers: " << threads <<endl;
            omp_set_num_threads(threads);
        }

};
} 
#endif
