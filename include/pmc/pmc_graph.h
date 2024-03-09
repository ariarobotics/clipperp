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

#ifndef PMC_GRAPH_H_
#define PMC_GRAPH_H_

#include <float.h>
#include <cstddef>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include <limits>
#include "math.h"
#include "pmc_headers.h"
#include "pmc_utils.h"
#include "pmc_vertex.h"


namespace pmc {
class pmc_graph {
private:
        // helper functions
        void read_mtx(const string& filename);
        void read_edges(const string& filename);
        void read_metis(const string& filename) const;

public:
        vector<long> edges;
        vector<long> vertices;
        vector<long> degree;
        long min_degree = 0;
        long max_degree = 0;
        double avg_degree = 0;
        bool is_gstats = false;
        long max_core = 0;
        string fn;
        vector<vector<bool>> adj;

        // constructor
        explicit pmc_graph(const string& filename)
                : fn(filename){
                read_graph(filename);
        };
        
        pmc_graph(bool graph_stats, const string& filename)
                : is_gstats(graph_stats),
                  fn(filename){
                read_graph(filename);
        };

        pmc_graph(const string& filename, bool make_adj)
                : fn(filename) {
                read_graph(filename);
                if (make_adj) create_adj();
        };

        pmc_graph(vector<long> vs, vector<long> es) {
                edges = std::move(es);
                vertices = std::move(vs);
                vertex_degrees();
        }

        pmc_graph(long nedges, const int *ei, const int *ej, int offset);
        explicit pmc_graph(map<int,vector<int>> const & v_map);
        
        // destructor
        ~pmc_graph() = default;

        void read_graph(const string& filename);
        void create_adj();
        void reduce_graph(int* &pruned);
        void reduce_graph(
                vector<long>& vs,
                vector<int>& es,
                int* &pruned,
                int id,
                const int& mc) const;

        size_t num_vertices() const { return vertices.size() - 1; }
        size_t num_edges() const { return edges.size()/2; }
        vector <long>* get_vertices(){ return &vertices; }
        vector<long>* get_edges(){ return &edges; }
        vector<long>* get_degree(){ return &degree; }
        vector<long> get_edges_array() const { return edges; }
        vector<long> get_vertices_array() const { return vertices; };
        vector<long> e_v;
        vector<long> e_u;
        vector<long> eid;

        long vertex_degree(long v) { return vertices[v] - vertices[v+1]; }
        long first_neigh(int v) { return vertices[v]; }
        long last_neigh(int v) { return vertices[v+1]; }

        void sum_vertex_degrees();
        void vertex_degrees();
        void update_degrees();
        void update_degrees(bool flag);
        void update_degrees(int* &pruned, const int& mc);
        double density() const; 
        long get_max_degree() const { return max_degree; }
        long get_min_degree() const { return min_degree; }
        double get_avg_degree() const { return avg_degree; }

        string get_file_extension(const string& filename) const;
        void basic_stats (double sec) const;
        void bound_stats(int alg) const;

        // vertex sorter
        void compute_ordering(vector<int>& bound, vector<long>& order) const;
        void compute_ordering(string degree, vector<long>& order);
        // edge sorters
        void degree_bucket_sort();
        void degree_bucket_sort(bool desc);

        
        vector<long> kcore;
        vector<long> kcore_order;
        vector<long>* get_kcores() { return &kcore; }
        vector<long>* get_kcore_ordering() { return &kcore_order; }
        long get_max_core() const { return max_core; }
        void update_kcores(int* &pruned);

        void compute_cores();
        void induced_cores_ordering(
                vector<long>& V,
                vector<int>& E,
                int* &pruned);

        // clique utils
        long initial_pruning(pmc_graph& G, int* &pruned, int lb);
        long initial_pruning(pmc_graph& G, int* &pruned, int lb, vector<vector<bool>> &adj);
        void order_vertices(vector<Vertex> &V, pmc_graph &G,
                const int &lb_idx, const long &lb, const string& vertex_ordering, bool decr_order);

        void print_info(const vector<int> &C_max, const double &sec) const;
        void print_break() const;
        bool time_left(const vector<int> &C_max, double sec,
                double time_limit, bool &time_expired_msg) const;
        void graph_stats(const pmc_graph& G, const int& mc, const int id, const double& sec) const;

        void reduce_graph(
                vector<long>& vs,
                vector<int>& es,
                int* &pruned,
                pmc_graph& G,
                int id,
                const int& mc) const;

        bool clique_test(pmc_graph& G, vector<int> C) const;
};

}
#endif
