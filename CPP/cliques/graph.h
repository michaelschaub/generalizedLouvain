#pragma once

#include <cstddef>
#include <vector>

namespace clq {

/**
 @brief  Undirected weighted graph, stored as a list of edges.

 The node ids are 0 to num_nodes() - 1. Self-loops and parallel edges are
 permitted. The Louvain move step does not use this list directly: it uses
 the CSR neighbour cache (see build_neighbour_cache in graphhelpers.h).
 */
class Graph {
public:
    explicit Graph(int num_nodes = 0) : num_nodes_(num_nodes) {}

    // Add a node and return its id.
    int add_node() { return num_nodes_++; }

    void add_edge(int u, int v, double weight)
    {
        u_.push_back(u);
        v_.push_back(v);
        weight_.push_back(weight);
    }

    void reserve_edges(std::size_t num_edges)
    {
        u_.reserve(num_edges);
        v_.reserve(num_edges);
        weight_.reserve(num_edges);
    }

    int num_nodes() const { return num_nodes_; }
    std::size_t num_edges() const { return weight_.size(); }

    int u(std::size_t edge) const { return u_[edge]; }
    int v(std::size_t edge) const { return v_[edge]; }
    double weight(std::size_t edge) const { return weight_[edge]; }

private:
    int num_nodes_;
    std::vector<int> u_;
    std::vector<int> v_;
    std::vector<double> weight_;
};

/**
 @brief  Call f(u, v, weight) for each edge, from the last edge added to the first.

 This is the order of LEMON's SmartGraph, which this library used before.
 Floating-point sums and the neighbour order depend on this order. Keep it to
 get the same results as before for a given seed.
 */
template<typename F>
void for_each_edge(const Graph &graph, F f)
{
    for (std::size_t e = graph.num_edges(); e-- > 0;) {
        f(graph.u(e), graph.v(e), graph.weight(e));
    }
}

}
