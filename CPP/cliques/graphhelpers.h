#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>
#include <math.h>
#include "graph.h"

namespace clq {

// Compressed-sparse-row adjacency cache built once per Louvain level.
// Self-loops are stored separately so the per-node move kernel does not need
// the inner branch on `u != v`.
struct NeighbourCache {
    std::vector<std::size_t> start;   // size n+1, prefix sum of non-self degrees
    std::vector<int>         node;    // size = total non-self incidences
    std::vector<double>      weight;  // same size as `node`
    std::vector<double>      selfloop_weight;  // size n, per-node self-loop weight
};

// Build a CSR neighbour cache from the edge list. The neighbours of each node
// are in the order of for_each_edge (last edge added first). A self-loop
// gives two halves of its weight, the same as a node with two edge ends.
inline NeighbourCache build_neighbour_cache(const Graph &graph) {
    const std::size_t n = static_cast<std::size_t>(graph.num_nodes());
    NeighbourCache cache;
    cache.start.assign(n + 1, 0);
    cache.selfloop_weight.assign(n, 0.0);

    // First pass: count non-self degrees and accumulate self-loop weights.
    for_each_edge(graph, [&](int u, int v, double w) {
        if (u != v) {
            ++cache.start[u + 1];
            ++cache.start[v + 1];
        } else {
            cache.selfloop_weight[u] += w / 2.0;
            cache.selfloop_weight[u] += w / 2.0;
        }
    });
    for (std::size_t i = 0; i < n; ++i) {
        cache.start[i + 1] += cache.start[i];
    }
    const std::size_t total = cache.start[n];
    cache.node.assign(total, 0);
    cache.weight.assign(total, 0.0);

    // Second pass: fill neighbour and weight arrays in CSR order.
    std::vector<std::size_t> cursor(cache.start.begin(), cache.start.begin() + n);
    for_each_edge(graph, [&](int u, int v, double w) {
        if (u == v) return;
        std::size_t k = cursor[u]++;
        cache.node[k] = v;
        cache.weight[k] = w;
        k = cursor[v]++;
        cache.node[k] = u;
        cache.weight[k] = w;
    });
    return cache;
}

// Weighted degree of each node. A self-loop counts once.
inline std::vector<double> find_weighted_degrees(const Graph &graph) {
	std::vector<double> degree(graph.num_nodes(), 0.0);
	for_each_edge(graph, [&](int u, int v, double w) {
		if (u != v) {
			degree[u] += w;
			degree[v] += w;
		} else {
			degree[u] += w / 2;
			degree[u] += w / 2;
		}
	});
	return degree;
}

inline double find_total_weight(const Graph &graph) {
	double total_weight = 0.0;
	for_each_edge(graph, [&](int u, int v, double w) {
		if (u != v) {
			total_weight = total_weight + w;
		} else {
			total_weight = total_weight + w / 2;
		}
	});
	return total_weight;
}

template<typename P, typename I>
void create_reduced_graph_from_partition(Graph &reduced_graph, const Graph &graph,
		P &partition, const std::vector<int> &new_comm_id_to_old_comm_id, I &internals) {

	// get number of communities
	int num_comm = partition.set_count();

	reduced_graph = Graph(num_comm);

	// Accumulate inter-community weights into a hash map keyed by the
	// canonical community pair (min, max). One streaming pass over the
	// original edges, then one batch insert into the reduced graph.
	std::unordered_map<std::uint64_t, double> pair_to_weight;
	pair_to_weight.reserve(static_cast<std::size_t>(num_comm) * 2u);

	for_each_edge(graph, [&](int u, int v, double w) {
		int cu = partition.find_set(u);
		int cv = partition.find_set(v);

		// internal weights already accounted for via self-loops below
		if (cu == cv) {
			return;
		}

		std::uint32_t lo = static_cast<std::uint32_t>(std::min(cu, cv));
		std::uint32_t hi = static_cast<std::uint32_t>(std::max(cu, cv));
		std::uint64_t key = (static_cast<std::uint64_t>(hi) << 32) | lo;
		pair_to_weight[key] += w;
	});

	// Insert edges in a canonical (sorted) order: unordered_map iteration order
	// is implementation-defined (libstdc++ vs libc++), and the edge order sets
	// the neighbour order at the next level, which breaks ties between equal
	// gains. Sorting keeps results for a given seed identical across platforms.
	std::vector<std::pair<std::uint64_t, double>> sorted_pairs(pair_to_weight.begin(),
			pair_to_weight.end());
	std::sort(sorted_pairs.begin(), sorted_pairs.end(),
			[](const std::pair<std::uint64_t, double> &a,
			   const std::pair<std::uint64_t, double> &b) { return a.first < b.first; });

	reduced_graph.reserve_edges(static_cast<std::size_t>(num_comm) + sorted_pairs.size());

	// add self-loops in new_graph
	for (int i = 0; i < num_comm; i++) {
		int old_comm_id = new_comm_id_to_old_comm_id[i];
		reduced_graph.add_edge(i, i, internals.comm_w_in[old_comm_id]);
	}

	for (const auto &kv : sorted_pairs) {
		int lo = static_cast<int>(kv.first & 0xffffffffull);
		int hi = static_cast<int>(kv.first >> 32);
		reduced_graph.add_edge(lo, hi, kv.second);
	}
}

}
