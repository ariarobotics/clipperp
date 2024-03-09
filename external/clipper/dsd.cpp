/**
 * @file dsd.cpp
 * @brief Dense subgraph discovery using Goldberg's algorithm
 * @brief See https://github.com/MengLiuPurdue/find_densest_subgraph
 * @author Parker Lusk <plusk@mit.edu>
 * @date 13 June 2023
 */

#include "clipper/dsd.h"

/**
 * This code is due to https://github.com/MengLiuPurdue/find_densest_subgraph
 */

namespace clipper {
namespace dsd {

void init_edges( double (*edges_info)[3], const double *degree, int64_t n, int64_t m, int64_t src, int64_t dest, double g)
{
    long i;
    /*add edge from src to every other vertex*/
    for(i = m; i < m + n; i ++)
    {
        edges_info[i][0] = static_cast<double>(src);
        edges_info[i][1] = static_cast<double>(i - m + 1);
        edges_info[i][2] = static_cast<double>(m / 2);
    }
    /*add edge from every other vertex to src*/
    for(i = n + m; i < m + 2 * n; i ++)
    {
        edges_info[i][0] = static_cast<double>(i - m - n + 1);
        edges_info[i][1] = static_cast<double>(dest);
        edges_info[i][2] = static_cast<double>(m / 2 + 2 * g - degree[i - m - n]);
    }
}

// ----------------------------------------------------------------------------

void new_edge(int64_t u, int64_t v, double weight, int64_t *to, double *cap, double *flow, int64_t *next, int64_t *fin, long *nEdge)
{
    to[*nEdge] = v;
    cap[*nEdge] = weight;
    flow[*nEdge] = 0;
    next[*nEdge] = fin[u];
    fin[u] = (*nEdge) ++;
    to[*nEdge] = u;
    cap[*nEdge] = weight;
    flow[*nEdge] = weight;
    next[*nEdge] = fin[v];
    fin[v] = (*nEdge) ++;
}

// ----------------------------------------------------------------------------

bool dinic_bfs(int64_t nverts, int64_t src, int64_t dest, int64_t *dist, int64_t *Q, int64_t *fin, int64_t *next, int64_t *to, double *flow, double *cap) 
{
    int64_t st;
    int64_t en;
    int64_t i;
    int64_t u;
    int64_t v;
    while(st < en) 
    {
        u = Q[st ++];
        for(i = fin[u]; i >= 0; i = next[i]) 
        {
            v = to[i];
            if(flow[i] < cap[i] && dist[v] == -1) 
            {
                dist[v] = dist[u] +1;
                Q[en ++] = v;
            }
        }
    }
    return dist[dest] != -1;
}

// ----------------------------------------------------------------------------

double dinic_dfs(int64_t u, double fl, int64_t src, int64_t dest, int64_t *pro, int64_t *next, int64_t *to, int64_t *dist, double *cap, double *flow) 
{
    if(u == dest) return fl;
    int64_t v;
    double df;
    for(int64_t &e=pro[u]; e >= 0; e = next[e]) 
    {
        v = to[e];
        if(flow[e] < cap[e] && dist[v] == dist[u] + 1) 
        {
            if(u == src || (cap[e] - flow[e]) <= fl)
            {
                fl = cap[e] - flow[e];
            }
            df = dinic_dfs(v, fl, src, dest, pro, next, to, dist, cap, flow);
            if(df>0) 
            {
                flow[e] += df;
                flow[e^1] -= df;
                return df;
            }
        }
    }
    return 0;
}

// ----------------------------------------------------------------------------

void find_cut(int64_t u, int64_t *cut, int64_t *another_pro, int64_t *next, int64_t *to, double *flow, double *cap)
{
    cut[u] = 1;
    for(int64_t &e = another_pro[u]; e >= 0; e = next[e])
    {
        int64_t v = to[e];
        if(flow[e] < cap[e] && cut[v] == 0)
        {
            find_cut(v, cut, another_pro, next, to, flow, cap);
        }
    }
}

// ----------------------------------------------------------------------------

double max_flow(double (*edges_info)[3], int64_t nverts, int64_t nedges, int64_t src, int64_t dest, int64_t *Q, int64_t *fin, int64_t *pro, int64_t *dist, int64_t *next, int64_t *to, int64_t *cut, int64_t *another_pro, int64_t *pro3, double *flow, double *cap)
{
    long i;
    std::fill(fin, fin + nverts, -1);
    std::fill(cut, cut + nverts, 0);
    long nEdge = 0;
    for(i = 0; i < nedges; i ++)
    {
        new_edge(edges_info[i][0], edges_info[i][1], edges_info[i][2], to, cap, flow, next, fin, &nEdge);
    }

    double ret = 0;
    double df;
    while(dinic_bfs(nverts, src, dest, dist, Q, fin, next, to, flow, cap)) 
    {
        for(i = 0; i < nverts; i++) 
        {
            pro[i] = fin[i];
            another_pro[i] = fin[i];
            pro3[i] = fin[i];
        }
        while(true) 
        {
            df = dinic_dfs(src, 0, src, dest, pro, next, to, dist, cap, flow);
            if(df) ret += df;
            else break;
        }
    }
    find_cut(src, cut, another_pro, next, to, flow, cap);
    return ret;
}

// ----------------------------------------------------------------------------

/**
 * Compute the densest subgraph given a graph in edge-list format.
 *
 * The edge list format needs to be symmetric, so if (ei[k1],ej[k1]) = (r,s) then there must be (ei[k2],ej[k2]) = (s,r),
 * unless r==s.
 * 
 * @param [in] n the number of nodes of the graph
 * @param [in] m the number of edges of the graph (and also the length of arrays ei, ej, w)
 * @param [in] ei a list of sources for each edge (0 <= ei[i] <= n-1) of length m
 * @param [in] ej a list of destinations for each edge (0 <= ei[i] <= n-1) of length m
 * @param [in] w a list of non-negative weights for each edge (0 <= w[i]) of length m
 * @param [out] output a list of vertices in the densest subgraph used, 0 <= output[i] <= n-1 for 0 <= i <= outputlen-1, this array must be capable of holding n
 * @param [in/out] the valid length of the output list and the length of the set of output used. 
 * @return the density of the densest subgraph in total_edge_weight/number of vertices
 */
double densest_subgraph(int64_t n, int64_t m, int64_t *ei, int64_t *ej, double *w, int64_t *output, size_t *outputlen)
{
    int64_t nverts;
    int64_t nedges;
    int64_t src;
    int64_t dest;
    double g;
    double maxflow;
    nverts = n + 2;
    nedges = m + 2 * n;
    int64_t i;
    auto edges_info = (double (*)[3])malloc(sizeof(double) * nedges * 3);
    auto degree = (double *)malloc(sizeof(double) * n); // degree of every vertex
    std::fill(degree, degree + n, 0);
    for (i = 0; i < m; i++)
    {
        edges_info[i][0] = static_cast<double>(ei[i]); // Cast int64_t to double
        edges_info[i][1] = static_cast<double>(ej[i]); // Cast int64_t to double
        edges_info[i][2] = w[i];
        degree[static_cast<int64_t>(edges_info[i][0])] += edges_info[i][2];
        edges_info[i][0]++;
        edges_info[i][1]++;
    }

    double final_degree = 0; // the degree of the densest subgraph
    double L = 0;            // lower bound of Andrew Goldberg's algorithm
    double U = m / 2;        // upper bound of Andrew Goldberg's algorithm
    auto final_cut = (int64_t *)malloc(sizeof(int64_t) * nverts);

    int64_t iter = 0;

    /*malloc enough space to store data used to calculate max flow*/
    auto Q = (int64_t *)malloc(sizeof(int64_t) * nverts);
    auto fin = (int64_t *)malloc(sizeof(int64_t) * nverts);
    auto pro = (int64_t *)malloc(sizeof(int64_t) * nverts);
    auto another_pro = (int64_t *)malloc(sizeof(int64_t) * nverts);
    auto pro3 = (int64_t *)malloc(sizeof(int64_t) * nverts);
    auto dist = (int64_t *)malloc(sizeof(int64_t) * nverts);
    auto flow = (double *)malloc(sizeof(double) * 2 * nedges);
    auto cap = (double *)malloc(sizeof(double) * 2 * nedges);
    auto next = (int64_t *)malloc(sizeof(int64_t) * 2 * nedges);
    auto to = (int64_t *)malloc(sizeof(int64_t) * 2 * nedges);
    auto cut = (int64_t *)malloc(sizeof(int64_t) * nverts);

    /*Andrew Goldberg's algorithm*/
    while(n * (n - 1) * (U - L) >= 1)
    {
        iter ++;
        g = (U + L) / 2;
        src = 0;
        dest = nverts - 1;
        init_edges(edges_info, degree, n, m, src, dest, g);
        maxflow = max_flow(edges_info, nverts, nedges, src, dest, Q, fin, pro, dist, next, to, cut, another_pro, pro3, flow, cap);
        if (std::accumulate(cut, cut + nverts, 0) == 1) {
            U = g;
        } else {
            L = g;
            for(i = 0; i < nverts; i++) {
                final_cut[i] = cut[i];
            }
        }
    }

    final_cut[0] = 0;
    final_cut[nverts - 1] = 0;

    /*retrieve the densest subgraph from the final_cut*/
    int64_t num = 0;
    for(i = 1; i < nverts - 1; i++) {
        if (final_cut[i] != 0) {
            output[num ++] = i - 1;
            for (int64_t &e = pro3[i]; e >= 0; e = next[e]) {
                if (final_cut[to[e]] != 0) {
                    final_degree += cap[e];
                }
            }
        }
    }
    final_degree /= (2.0 * num);
    *outputlen = num;

    /*free space*/
    free(final_cut);
    free(degree);
    free(edges_info);
    free(Q);
    free(fin);
    free(pro);
    free(another_pro);
    free(pro3);
    free(dist);
    free(flow);
    free(cap);
    free(next);
    free(to);
    free(cut);
    return final_degree;
}

// ----------------------------------------------------------------------------

std::vector<long> solve(const SpAffinity& A, const std::vector<long>& _S)
{
    // allows the search for densest subgraph in A to be restricted to a
    // subgraph of A
    std::vector<long> S;
    if (!_S.empty()) {
        S = _S;
    } else {
        S.resize(A.rows());
        std::iota(S.begin(), S.end(), 0);
    }

    const int64_t n = A.rows(); // num nodes
    const int64_t m = S.size() * S.size() - S.size(); // num edges

    std::vector<int64_t> ei;
    std::vector<int64_t> ej;
    std::vector<double> w;
    ei.reserve(m);
    ej.reserve(m);
    w.reserve(m);

    for (const long i : S) {
        for (const long j : S) {
            // skip diagonals
            if (i == j) continue;

            // create a fully connected graph with weights
            ei.push_back(i);
            ej.push_back(j);

            // A is assumed symmetric and that the upper triangle is filled in.
            const double ew = (i < j) ? A.coeff(i,j) : A.coeff(j,i);
            w.push_back(ew);
        }
    }

    std::vector<int64_t> output(n, 0);
    size_t outlen = n;
    densest_subgraph(n, m, ei.data(), ej.data(), w.data(), output.data(), &outlen);

    std::vector<long> nodes;
    nodes.reserve(outlen);
    for (size_t i=0; i<outlen; ++i) {
        nodes.push_back(output[i]);
    }
    return nodes;
}

// ----------------------------------------------------------------------------

std::vector<long> solve(const Eigen::MatrixXd& A, const std::vector<long>& S)
{
    return solve(SpAffinity(A.sparseView()), S);
}

} // ns dsd
} // ns clipper
