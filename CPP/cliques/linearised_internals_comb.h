#pragma once

#include "graphhelpers.h"

namespace clq
{
// Define internal structure to carry statistics for combinatorial stability

struct LinearisedInternalsComb {

    typedef std::vector<double> range_map;

    unsigned int num_nodes; // number of nodes in graph
    unsigned int num_nodes_init; // number of nodes in initial graph
    double two_m; // 2 times total weight
    range_map node_to_nr_nodes_init; // mapping: node_id to number of nodes it
    // represented in the original graph
    range_map comm_tot_nodes; // total number of nodes (wrt original graph) per community
    range_map comm_w_in; // weight inside each community

    std::vector<double> node_weight_to_communities; //mapping: node weight to each community
    std::vector<unsigned int> neighbouring_communities_list; // associated list of neighbouring communities

    // CSR neighbour cache + per-node self-loop weight (see graphhelpers.h).
    NeighbourCache nbr;

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // simple constructor, no partition given; graph should be the original graph in this case
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    explicit LinearisedInternalsComb( const Graph& graph ) :
        num_nodes( graph.num_nodes() ), num_nodes_init( num_nodes ),
        node_to_nr_nodes_init( num_nodes, 1 ), comm_tot_nodes( num_nodes,
                1 ), comm_w_in( num_nodes, 0 ),
        node_weight_to_communities( num_nodes, 0 ),
        neighbouring_communities_list(),
        nbr( build_neighbour_cache( graph ) )
    {
        two_m = 2 * find_total_weight( graph );

        for( unsigned int i = 0; i < num_nodes; ++i ) {
            comm_w_in[i] = nbr.selfloop_weight[i];
        }
    }
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    // full constructor with reference to original partition
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    template<typename P>
    LinearisedInternalsComb( const Graph& graph, P& partition, P& partition_init ) :
        num_nodes( graph.num_nodes() ),
        node_to_nr_nodes_init( num_nodes, 0 ), comm_tot_nodes( num_nodes, 0 ),
        comm_w_in( num_nodes, 0 ),
        node_weight_to_communities( num_nodes, 0 ),
        neighbouring_communities_list(),
        nbr( build_neighbour_cache( graph ) )
    {
        two_m = 2 * find_total_weight( graph );
        num_nodes_init = partition_init.element_count(); // get original number of nodes

        // Get number of nodes incorporated into each supernode;
        // assumes partition is ordered/relabeled accordingly
        for( unsigned int i = 0; i < num_nodes_init; ++i ) {
            int old_comm_id = partition_init.find_set( i );
            node_to_nr_nodes_init[old_comm_id]++;
            int new_comm_id = partition.find_set( old_comm_id );
            comm_tot_nodes[new_comm_id]++;
        }

        for_each_edge( graph, [&]( int node_u_id, int node_v_id, double weight ) {
            int comm_of_node_u = partition.find_set( node_u_id );
            int comm_of_node_v = partition.find_set( node_v_id );

            if( node_u_id == node_v_id ) {
                weight = weight / 2;
            }

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
void isolate_and_update_internals( int node_id, LinearisedInternalsComb& internals, P& partition )
{
    int comm_id = partition.find_set( node_id );

    while( !internals.neighbouring_communities_list.empty() ) {
        unsigned int old_neighbour = internals.neighbouring_communities_list.back();
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

    internals.comm_tot_nodes[comm_id] -= internals.node_to_nr_nodes_init[node_id];
    internals.comm_w_in[comm_id] -= 2 * internals.node_weight_to_communities[comm_id]
                                    + internals.nbr.selfloop_weight[node_id];

    partition.isolate_node( node_id );
}

/**
 @brief  insert a node into the best set & update internals
 */
template<typename P>
void insert_and_update_internals( int node_id, LinearisedInternalsComb& internals, P& partition,
                                  int best_comm )
{
    internals.comm_tot_nodes[best_comm] += internals.node_to_nr_nodes_init[node_id];
    internals.comm_w_in[best_comm] += 2 * internals.node_weight_to_communities[best_comm]
                                      + internals.nbr.selfloop_weight[node_id];

    partition.add_node_to_set( node_id, best_comm );
}

}
