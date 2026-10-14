#pragma once

#include "graphhelpers.h"
#include "io.h"

namespace clq
{
// Define internal structure to carry statistics for generalised stability

struct LinearisedInternalsGeneric {


    // typedef for convenience
    typedef std::vector<std::vector<double>> vec_of_vec;
    unsigned int num_nodes;

    // how many outer products do we have? must be multiple of 2
    unsigned int num_null_model_vectors;
    // contains the outer product vectors of the null model
    vec_of_vec null_model_vectors;

    // loss for each community (second matrix / null model terms per community)
    vec_of_vec comm_loss_vectors;
    // gain for each community (first matrix / gain term per community)
    std::vector<double> comm_w_in;
    // mapping: node weight to each community (gain of adding node to a community, dyn. updated)
    std::vector<double> node_weight_to_communities;
    // associated list of neighbouring communities
    std::vector<unsigned int> neighbouring_communities_list;

    // CSR neighbour cache + per-node self-loop weight, built once at
    // construction. Self-loops are excluded from the CSR — looked up via
    // selfloop_weight[node_id]. The hot move kernel walks the CSR directly.
    NeighbourCache nbr;

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // simple constructor, no partition given
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    LinearisedInternalsGeneric (const Graph& graph, vec_of_vec null_model_input) :
        num_nodes (graph.num_nodes() ),
        num_null_model_vectors (null_model_input.size() ),
        null_model_vectors (null_model_input),
        comm_loss_vectors (null_model_input),
        comm_w_in ( num_nodes, 0),
        node_weight_to_communities (num_nodes, 0), neighbouring_communities_list(),
        nbr (build_neighbour_cache (graph))
    {
        if (num_null_model_vectors % 2 != 0 ) {
            clq::output ("Null model vectors must be provided as pairs!");
            exit (EXIT_FAILURE);
        }

        for (unsigned int i = 0; i < num_nodes; ++i) {
            comm_w_in[i] = nbr.selfloop_weight[i];
        }
    }

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // full constructor with reference to partition
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    template<typename P>
        LinearisedInternalsGeneric (const Graph& graph, P& partition, vec_of_vec null_model_input) :
        num_nodes (graph.num_nodes() ),
        num_null_model_vectors (null_model_input.size() ),
        null_model_vectors (null_model_input),
        comm_loss_vectors (num_null_model_vectors, std::vector<double> (num_nodes, 0) ),
        comm_w_in (num_nodes, 0),
        node_weight_to_communities (num_nodes, 0), neighbouring_communities_list(),
        nbr (build_neighbour_cache (graph))
    {
        if (num_null_model_vectors % 2 != 0 ) {
            clq::output ("Null model vectors must be provided as pairs!");
            exit (EXIT_FAILURE);
        }

        // find internal statistics based on graph, weights and partitions
        // consider all edges
        for_each_edge (graph, [&] (int node_u_id, int node_v_id, double weight) {
            // this is to distinguish within community weight with total weight
            int comm_of_node_u = partition.find_set (node_u_id);
            int comm_of_node_v = partition.find_set (node_v_id);

            // if selfloop, only half of the weight has to be considered (multiplication by two afterwards)
            if (node_u_id == node_v_id) {
                weight = weight / 2;
            }

            if (comm_of_node_u == comm_of_node_v) {
                // in case the weight stems from within the community add to internal weights
                comm_w_in[comm_of_node_u] += 2 * weight;
            }
        });

        // setup internal loss vectors (summing up all the null model terms per group)
        for (unsigned int i = 0; i < num_nodes; ++i) {
            int comm_id = partition.find_set (i);

            for (unsigned int k = 0; k < num_null_model_vectors; ++k) {
                comm_loss_vectors[k][comm_id] += null_model_vectors[k][i];
            }
        }
    }
};


/**
 @brief  isolate a node into its singleton set & update internals
 */
template<typename P>
void isolate_and_update_internals (int node_id, LinearisedInternalsGeneric& internals, P& partition)
{
    int comm_id = partition.find_set (node_id);

    // reset weights
    while (!internals.neighbouring_communities_list.empty() ) {
        unsigned int old_neighbour =
            internals.neighbouring_communities_list.back();
        internals.neighbouring_communities_list.pop_back();
        internals.node_weight_to_communities[old_neighbour] = 0;
    }

    // CSR walk over node's non-self neighbours.
    const std::size_t beg = internals.nbr.start[node_id];
    const std::size_t end = internals.nbr.start[node_id + 1];
    for (std::size_t e = beg; e < end; ++e) {
        int comm_node = partition.find_set (internals.nbr.node[e]);
        if (internals.node_weight_to_communities[comm_node] == 0) {
            internals.neighbouring_communities_list.push_back (comm_node);
        }
        internals.node_weight_to_communities[comm_node] += internals.nbr.weight[e];
    }

    for (unsigned int j = 0; j < internals.num_null_model_vectors; ++j) {
        internals.comm_loss_vectors[j][comm_id] -= internals.null_model_vectors[j][node_id];
    }

    internals.comm_w_in[comm_id] -= 2 * internals.node_weight_to_communities[comm_id]
                                    + internals.nbr.selfloop_weight[node_id];

    partition.isolate_node (node_id);
}


/**
 @brief  insert a node into the best set & update internals
 */
template<typename P>
void insert_and_update_internals (int node_id, LinearisedInternalsGeneric& internals, P& partition,
                                  int best_comm)
{
    // update loss
    for (unsigned int j = 0; j < internals.num_null_model_vectors; ++j) {
        internals.comm_loss_vectors[j][best_comm] += internals.null_model_vectors[j][node_id];
    }

    // update gain
    internals.comm_w_in[best_comm] += 2 * internals.node_weight_to_communities[best_comm]
                                      + internals.nbr.selfloop_weight[node_id];

    partition.add_node_to_set (node_id, best_comm);
}

}
