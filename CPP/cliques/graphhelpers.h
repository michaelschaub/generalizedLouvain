#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <lemon/maps.h>
#include <lemon/core.h>
#include <math.h>

namespace clq {

// Compressed-sparse-row adjacency cache built once per Louvain level.
// Self-loops are stored separately so the per-node move kernel does not need
// the inner branch on `graph.u(e) != graph.v(e)`.
struct NeighbourCache {
    std::vector<std::size_t> start;   // size n+1, prefix sum of non-self degrees
    std::vector<int>         node;    // size = total non-self incidences
    std::vector<double>      weight;  // same size as `node`
    std::vector<double>      selfloop_weight;  // size n, per-node self-loop weight
};

// Build a CSR neighbour cache from a LEMON graph by walking IncEdgeIt once
// per node. Self-loop weights are accumulated separately. IncEdgeIt visits
// self-loops twice (once per endpoint) so we halve the contribution to match
// the convention used elsewhere in the code.
template<typename G, typename M>
NeighbourCache build_neighbour_cache(G &graph, M &weights) {
    const std::size_t n = static_cast<std::size_t>(lemon::countNodes(graph));
    NeighbourCache cache;
    cache.start.assign(n + 1, 0);
    cache.selfloop_weight.assign(n, 0.0);

    // First pass: count non-self degrees and accumulate self-loop weights.
    for (std::size_t i = 0; i < n; ++i) {
        typename G::Node u = graph.nodeFromId(static_cast<int>(i));
        std::size_t deg = 0;
        for (typename G::IncEdgeIt e(graph, u); e != lemon::INVALID; ++e) {
            if (graph.u(e) != graph.v(e)) {
                ++deg;
            } else {
                // IncEdgeIt visits self-loops twice; accumulate half each time.
                cache.selfloop_weight[i] += weights[e] / 2.0;
            }
        }
        cache.start[i + 1] = deg;
    }
    for (std::size_t i = 0; i < n; ++i) {
        cache.start[i + 1] += cache.start[i];
    }
    const std::size_t total = cache.start[n];
    cache.node.assign(total, 0);
    cache.weight.assign(total, 0.0);

    // Second pass: fill neighbour and weight arrays in CSR order.
    std::vector<std::size_t> cursor(cache.start.begin(), cache.start.begin() + n);
    for (std::size_t i = 0; i < n; ++i) {
        typename G::Node u = graph.nodeFromId(static_cast<int>(i));
        for (typename G::IncEdgeIt e(graph, u); e != lemon::INVALID; ++e) {
            if (graph.u(e) == graph.v(e)) continue;
            typename G::Node v = graph.oppositeNode(u, e);
            std::size_t k = cursor[i]++;
            cache.node[k] = graph.id(v);
            cache.weight[k] = weights[e];
        }
    }
    return cache;
}

template<typename G>
int find_unweighted_degree(G &graph, int node_id) {
	int count = 0;
	typename G::Node n = graph.nodeFromId(node_id);

	for (typename G::IncEdgeIt e(graph, n); e != lemon::INVALID; ++e) {
		++count;
	}
	return count;
}

template<typename G, typename M, typename NO>
double find_weighted_degree(G &graph, M &weights, NO node) {
	double degree = 0.0;
	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			degree += weights[e];
		} else {
			degree += weights[e] / 2; // self-loop is iterated twice with IncEdgeIt
		}
	}
	return degree;
}

template<typename G, typename M, typename NO>
double find_weight_selfloops(G &graph, M &weights, NO node) {
	typename G::Edge edge = lemon::findEdge(graph, node, node);
	if (edge == lemon::INVALID) {
		return 0.0;
	} else {
		return double(weights[edge]);
	}
}

template<typename G, typename M>
double find_total_weight(G &graph, M &weights) {
	double total_weight = 0.0;
	for (typename G::EdgeIt e(graph); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			total_weight = total_weight + weights[e];
		} else {
			total_weight = total_weight + weights[e] / 2;
		}
	}
	return total_weight;
}

template<typename G, typename P, typename W, typename I>
void create_reduced_graph_from_partition(G &reduced_graph,
		W &reduced_weight_map, G &graph, W &weights, P &partition,
		const std::vector<int> &new_comm_id_to_old_comm_id, I &internals) {

	typedef typename G::EdgeIt EdgeIt;
	typedef typename G::Edge Edge;
	typedef typename G::Node Node;

	// get number of communities
	int num_comm = partition.set_count();

	reduced_graph.reserveNode(num_comm);

	// add self-loops in new_graph
	for (int i = 0; i < num_comm; i++) {
		Node comm_node = reduced_graph.addNode();
		Edge e = reduced_graph.addEdge(comm_node, comm_node);
		int old_comm_id = new_comm_id_to_old_comm_id[i];
		reduced_weight_map[e] = internals.comm_w_in[old_comm_id];
	}

	// Accumulate inter-community weights into a hash map keyed by the
	// canonical community pair (min, max). One streaming pass over the
	// original edges, then one batch insert into the reduced graph —
	// avoids the O(deg) lemon::findEdge per inter-community edge.
	std::unordered_map<std::uint64_t, double> pair_to_weight;
	pair_to_weight.reserve(static_cast<std::size_t>(num_comm) * 2u);

	for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
		int cu = partition.find_set(graph.id(graph.u(edge)));
		int cv = partition.find_set(graph.id(graph.v(edge)));

		// internal weights already accounted for via self-loops above
		if (cu == cv) {
			continue;
		}

		std::uint32_t lo = static_cast<std::uint32_t>(std::min(cu, cv));
		std::uint32_t hi = static_cast<std::uint32_t>(std::max(cu, cv));
		std::uint64_t key = (static_cast<std::uint64_t>(hi) << 32) | lo;
		pair_to_weight[key] += weights[edge];
	}

	// Insert edges in a canonical (sorted) order: unordered_map iteration order
	// is implementation-defined (libstdc++ vs libc++), and the edge order sets
	// the neighbour order at the next level, which breaks ties between equal
	// gains. Sorting keeps results for a given seed identical across platforms.
	std::vector<std::pair<std::uint64_t, double>> sorted_pairs(pair_to_weight.begin(),
			pair_to_weight.end());
	std::sort(sorted_pairs.begin(), sorted_pairs.end(),
			[](const std::pair<std::uint64_t, double> &a,
			   const std::pair<std::uint64_t, double> &b) { return a.first < b.first; });

	for (const auto &kv : sorted_pairs) {
		std::uint32_t lo = static_cast<std::uint32_t>(kv.first & 0xffffffffull);
		std::uint32_t hi = static_cast<std::uint32_t>(kv.first >> 32);
		Node nu = reduced_graph.nodeFromId(static_cast<int>(lo));
		Node nv = reduced_graph.nodeFromId(static_cast<int>(hi));
		Edge e = reduced_graph.addEdge(nu, nv);
		reduced_weight_map[e] = kv.second;
	}
}

// find weight to all communities excluding self-loops (to own community)
template<typename G, typename P, typename W, typename NO>
std::map<int, double> find_weight_node_to_communities(G &graph, P &partition,
		W &weights, NO node) {
	std::map<int, double> community_to_weight;

	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			double edge_weight = weights[e];
			NO opposite_node = graph.oppositeNode(node, e);
			int comm_node = partition.find_set(graph.id(opposite_node));
			community_to_weight[comm_node] += edge_weight;
		}
	}
	//clq::print_map(community_to_weight);
	return community_to_weight;
}

// find weight to specific community excluding self-loops (to own community)
template<typename G, typename P, typename W, typename NO>
double find_weight_node_to_community(G &graph, P &partition, W &weights,
		NO node, int comm_id) {

	double summed_weight = 0.0;
	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			double edge_weight = weights[e];
			NO opposite_node = graph.oppositeNode(node, e);
			int comm_node = partition.find_set(graph.id(opposite_node));
			if (comm_node == comm_id) {
				summed_weight += edge_weight;
			}
		}
	}

	return summed_weight;
}

}
