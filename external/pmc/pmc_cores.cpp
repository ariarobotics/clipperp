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

#include "pmc/pmc_graph.h"

using namespace pmc;

void pmc_graph::induced_cores_ordering(
        vector<long>& V,
        vector<int>& E,
        int* &pruned) {

    long n;
    long d;
    long start;
    long num;
    long md;
    long v;
    long u;
    long w;
    long du;
    long pu;
    long pw;
    long md_end;
    n = vertices.size();

    vector <long> pos_tmp(n);
    vector <long> core_tmp(n);
    vector <long> order_tmp(n);

    md = 0;
    for(v=1; v<n; v++) {
        core_tmp[v] = V[v] - V[v-1];
        if (core_tmp[v] > md)  md = core_tmp[v];
    }

    md_end = md+1;
    vector <long> bin(md_end,0);

    for (v=1; v < n; v++)  bin[core_tmp[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    for (v=1; v<n; v++) {
        pos_tmp[v] = bin[core_tmp[v]];
        order_tmp[pos_tmp[v]] = v;
        bin[core_tmp[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    for (long i = 1; i < n; i++) {
        v=order_tmp[i];
        for (long j = V[v-1]; j < V[v]; j++) {
            u = E[j] + 1;
            if (core_tmp[u] > core_tmp[v]) {
                du = core_tmp[u];   
                pu = pos_tmp[u];
                pw = bin[du];       
                w = order_tmp[pw];
                if (u != w) {
                    pos_tmp[u] = pw;   
                    order_tmp[pu] = w;
                    pos_tmp[w] = pu;   
                    order_tmp[pw] = u;
                }
                bin[du]++;   
                core_tmp[u]--;
            }
        }
    }

    for (v=0; v<n-1; v++) {
        core_tmp[v] = core_tmp[v+1]+1;
        order_tmp[v] = order_tmp[v+1]-1;
    }

    kcore = core_tmp;
    kcore_order = order_tmp;
    bin.clear();
}




void pmc_graph::compute_cores() {
    long d;
    int start;
    int num;
    long md;
    long u;
    long w;
    long du;
    int pu;
    int pw;
    long md_end;
    long v;
    unsigned long n = vertices.size();

    vector <int> pos(n);
    if (!kcore_order.empty()) {
        vector<long> tmp(n,0);
        kcore = tmp;
        kcore_order = tmp;
    }
    else {
        kcore_order.resize(n);
        kcore.resize(n);
    }

    md = 0;
    for (v=1; v<n; v++) {
        kcore[v] = vertices[v] - vertices[v-1];
        if (kcore[v] > md)  md = kcore[v];
    }

    md_end = md+1;
    vector < int > bin(md_end,0);

    for (v=1; v < n; v++)  bin[kcore[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    // bucket sort
    for (v=1; v<n; v++) {
        pos[v] = bin[kcore[v]];
        kcore_order[pos[v]] = v;
        bin[kcore[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    // kcores
    for (int i=1; i<n; i++) {
        v=kcore_order[i];
        for (long j=vertices[v-1]; j<vertices[v]; j++) {
            u = edges[j] + 1;
            if (kcore[u] > kcore[v]) {
                du = kcore[u];   
                pu = pos[u];
                pw = bin[du];    
                w = kcore_order[pw];
                if (u != w) {
                    pos[u] = pw;   kcore_order[pu] = w;
                    pos[w] = pu;   kcore_order[pw] = u;
                }
                bin[du]++;   kcore[u]--;
            }
        }
    }

    for (v = 0; v < n-1; v++) {
        kcore[v] = kcore[v+1] + 1; // K + 1
        kcore_order[v] = kcore_order[v+1]-1;
    }
    max_core = kcore[kcore_order[num_vertices()-1]] - 1;

    bin.clear();
    pos.clear();
}
