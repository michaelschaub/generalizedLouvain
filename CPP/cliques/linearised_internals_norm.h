#pragma once

#include "graphhelpers.h"


namespace clq
{

struct LinearisedInternals {
    typedef std::vector<double> range_map;
    unsigned int num_nodes;
    double two_m;
    range_map node_to_w;
    range_map comm_w_tot;
    range_map comm_w_in;
    std::vector<double> node_weight_to_communities;
    std::vector<unsigned int> neighbouring_communities_list;

    // CSR neighbour cache + per-node self-loop weight (see graphhelpers.h).
    NeighbourCache nbr;

    explicit LinearisedInternals( const Graph& graph ) :
        num_nodes( graph.num_nodes() ), node_to_w( num_nodes, 0 ),
        comm_w_tot( num_nodes, 0 ), comm_w_in( num_nodes, 0 ),
        node_weight_to_communities( num_nodes, 0 ),
        neighbouring_communities_list(),
        nbr( build_neighbour_cache( graph ) )
    {
        two_m = 2 * find_total_weight( graph );

        std::vector<double> degree = find_weighted_degrees( graph );
        for( unsigned int i = 0; i < num_nodes; ++i ) {
            comm_w_tot[i] = node_to_w[i] = degree[i];
            comm_w_in[i] = nbr.selfloop_weight[i];
        }
    }

    template<typename P>
    LinearisedInternals( const Graph& graph, P& partition ) :
        num_nodes( graph.num_nodes() ), node_to_w( num_nodes, 0 ),
        comm_w_tot( num_nodes, 0 ), comm_w_in( num_nodes, 0 ),
        node_weight_to_communities( num_nodes, 0 ),
        neighbouring_communities_list(),
        nbr( build_neighbour_cache( graph ) )
    {
        two_m = 2 * find_total_weight( graph );

        // find internal statistics based on graph, weights and partitions
        // consider all edges
        for_each_edge( graph, [&]( int node_u_id, int node_v_id, double weight ) {
            int comm_of_node_u = partition.find_set( node_u_id );
            int comm_of_node_v = partition.find_set( node_v_id );

            if( node_u_id == node_v_id ) {
                weight = weight / 2;
            }

            node_to_w[node_u_id] += weight;
            node_to_w[node_v_id] += weight;
            comm_w_tot[comm_of_node_u] += weight;
            comm_w_tot[comm_of_node_v] += weight;

            if( comm_of_node_u == comm_of_node_v ) {
                comm_w_in[comm_of_node_u] += 2 * weight;
            }
        } );
    }
};

/**
 @brief  isolate a node into its singleton set & update internals
 */
template<typename P>
void isolate_and_update_internals( int node_id, LinearisedInternals& internals, P& partition )
{
    int comm_id = partition.find_set( node_id );

    // reset weights
    while( !internals.neighbouring_communities_list.empty() ) {
        unsigned int old_neighbour =
            internals.neighbouring_communities_list.back();
        internals.neighbouring_communities_list.pop_back();
        internals.node_weight_to_communities[old_neighbour] = 0;
    }

    // CSR walk over node's non-self neighbours.
    const std::size_t beg = internals.nbr.start[node_id];
    const std::size_t end = internals.nbr.start[node_id + 1];
    for( std::size_t e = beg; e < end; ++e ) {
        int comm_node = partition.find_set( internals.nbr.node[e] );
        if( internals.node_weight_to_communities[comm_node] == 0 ) {
            internals.neighbouring_communities_list.push_back( comm_node );
        }
        internals.node_weight_to_communities[comm_node] += internals.nbr.weight[e];
    }

    internals.comm_w_tot[comm_id] -= internals.node_to_w[node_id];
    internals.comm_w_in[comm_id] -= 2
                                    * internals.node_weight_to_communities[comm_id]
                                    + internals.nbr.selfloop_weight[node_id];

    partition.isolate_node( node_id );
}

/**
 @brief  insert a node into the best set & update internals
 */
template<typename P>
void insert_and_update_internals( int node_id, LinearisedInternals& internals, P& partition,
                                  int best_comm )
{
    internals.comm_w_tot[best_comm] += internals.node_to_w[node_id];
    internals.comm_w_in[best_comm] += 2
                                      * internals.node_weight_to_communities[best_comm]
                                      + internals.nbr.selfloop_weight[node_id];

    partition.add_node_to_set( node_id, best_comm );
}

}
