#pragma once

#include <algorithm>
#include <random>
#include <vector>
#include "graphhelpers.h"
#include "stability_gen.h"
#include "stability.h"
#include "linearised_internals_generic.h"
#include "linearised_internals_norm.h"
#include "linearised_internals_comb.h"


namespace clq
{

// For convenience let's get rid of some cumbersome notation..
typedef std::vector<std::vector<double>> vec2;

struct BestMove {
    unsigned int best_comm;
    double net_delta;  // change in global Q if best_comm != original; 0 otherwise
};

/**
 @brief  find_best_comm_move -- find the community with the largest gain in quality
 for a particular node, and return the resulting net change in global quality.

 The change in global Q from one full isolate→find_best→insert step equals
 `gain_at_best - gain_at_original` (both evaluated against the post-isolation
 internals). When the chosen community is the original, the delta is exactly
 zero — the partition is unchanged.
 */
template <typename QD, typename L>
BestMove find_best_comm_move( QD compute_quality_diff, unsigned int node_id, L &internals,
                              unsigned int comm_id )
{
    unsigned int num_neighbour_communities = internals.neighbouring_communities_list.size();

    unsigned int best_comm = comm_id;
    double best_gain = 0;
    double gain_at_orig = 0;  // gain of re-inserting back into the original community
                              // (0 when orig is empty after isolation, i.e. not in the neighbour list)

    for( unsigned int k = 0; k < num_neighbour_communities; ++k ) {
        unsigned int comm_id_neighbour = internals.neighbouring_communities_list[k];
        double gain = compute_quality_diff( internals, comm_id_neighbour, node_id );

        if( comm_id_neighbour == comm_id ) {
            gain_at_orig = gain;
        }

        if( gain > best_gain ) {
            best_comm = comm_id_neighbour;
            best_gain = gain;
        } else if( gain == best_gain && comm_id == comm_id_neighbour ) {
            best_comm = comm_id;
        }
    }

    double net_delta = ( best_comm == comm_id ) ? 0.0 : ( best_gain - gain_at_orig );
    return { best_comm, net_delta };
}


/**
 @brief fold_partition_into_orginal_graph_size -- given a partition for an aggregated graph, create 
 the corresponding partition for the graph in the original size.

 @param[in]  optimal_partitions -- vector of partitions for the different levels of the 
                                   Louvain algorithm   
 @param[in]  partition -- partition of aggregated graph, to be transformed/expanded to the original
                          graph size

 @param[out] (partition) -- expanded partition with size matching the original graph
 */
template <typename P>
P fold_partition_into_orginal_graph_size( std::vector<P>& optimal_partitions, P partition )
{
    int hierarchy = optimal_partitions.size();

    // if optimal partition is empty there is nothing to fold back..
    if( hierarchy == 0 ) {
        return partition;
    } else {
        // get size of partition one level below, i.e. number of nodes in original graph
        unsigned int original_number_of_nodes = optimal_partitions[hierarchy - 1].element_count();
        // create new empty partition of this size
        P partition_original_nodes( original_number_of_nodes );

        // loop over nodes one level below
        int old_comm, new_comm;

        for( unsigned int id = 0; id < original_number_of_nodes; id++ ) {
            // get the communities for each node one level below
            old_comm = optimal_partitions[hierarchy - 1].find_set( id );
            // use this as node_id in the current partition as old community id
            // is equivalent to new node id and read out new community id
            new_comm = partition.find_set( old_comm );
            // include pair (node, new community) id in the newly created partition
            partition_original_nodes.add_node_to_set( id, new_comm );
        }

        return partition_original_nodes;
    }
}


/**
 @brief create_reduced_null_model_vec -- created an aggregated set of 'null-model' vectors 
 to be used for a corresponding aggregated graph

 @param[in]  new_comm_id_to_old_comm_id -- mapping from new to old community ids
 @param[in]  internals -- internals of the current level
 @param[in]  num_null_model_vectors -- number of null model vectors to be aggregated.
 @param[in]  num_nodes_reduced_graph -- number of nodes in the reduced graph
 
 @param[out] (vec2) -- vector of vectors,  containing the reduced null model vectors
 */
vec2 create_reduced_null_model_vec( const std::vector<int>& new_comm_id_to_old_comm_id,
                                    const LinearisedInternalsGeneric& internals,
                                    unsigned int num_null_model_vectors, unsigned int num_nodes_reduced_graph )
{
    vec2 reduced_null_model_vec( num_null_model_vectors, std::vector<double> ( num_nodes_reduced_graph, 0 ) );

    for( unsigned int k = 0; k < num_nodes_reduced_graph; ++k ) {
        int old_comm_id = new_comm_id_to_old_comm_id[k];
        for( unsigned int j = 0; j < num_null_model_vectors; ++j ) {
            reduced_null_model_vec[j][k] = internals.comm_loss_vectors[j][old_comm_id];
        }
    }

    return reduced_null_model_vec;
}



/**
 @brief  generalise Louvain method - greedy algorithm to find community structure of a network.

 @param[in]  graph    --    weighted graph to partition
 @param[in]  null_model_vec -- null model vectors
 @param[in]  compute_quality -- functor to compute the quality function
 @param[in]  compute_quality_diff -- functor to compute change in quality function
 @param[in]  initial_partition -- partition to start from
 @param[in]  optimal_partitions -- set of optimal_partitions at each level of the algorithm (initially empty)
 @param[in]  minimum_improve -- minimum_improvement necessary in each iteration / stopping criterion

 @param[out]  (double) quality of the optimal partition

 */
template<typename P, typename QF, typename QFDIFF>
double find_optimal_partition_louvain_gen( const Graph& graph, const vec2& null_model_vec,
        QF compute_quality, QFDIFF compute_quality_diff, P initial_partition,
        std::vector<P>& optimal_partitions, double minimum_improve,
        std::mt19937& rng )
{
    // set up initial partitions
    P partition( initial_partition );

    bool do_construct_new_graph = false;

    LinearisedInternalsGeneric internals( graph, partition, null_model_vec );

    // Build the visit order (initially one slot per node, randomised below).
    // Start from the highest node id, the order of LEMON's NodeIt that this
    // code used before: the shuffle then gives the same order for a seed.
    std::vector<int> nodes_ordered_randomly;
    nodes_ordered_randomly.reserve( static_cast<std::size_t>( graph.num_nodes() ) );

    for( int node_id = graph.num_nodes() - 1; node_id >= 0; --node_id ) {
        nodes_ordered_randomly.push_back( node_id );
    }

    // Portable Fisher-Yates: bit-identical across libstdc++/libc++ given the
    // same mt19937 state. std::shuffle's internal distribution is
    // implementation-defined and would diverge between platforms.
    for( std::size_t i = nodes_ordered_randomly.size(); i > 1; --i ) {
        std::size_t j = rng() % i;
        std::swap( nodes_ordered_randomly[j], nodes_ordered_randomly[i - 1] );
    }

    // Queue-based fast-move (Leiden-style pruning): start with every node
    // dirty, in shuffled order, and re-queue a node's neighbours whenever it
    // moves. This is deterministic given the seed since the queue is FIFO and
    // the initial fill respects the shuffle. Runs strictly less work than the
    // classical "sweep all nodes" loop yet converges to the same fixed point.
    const std::size_t n = nodes_ordered_randomly.size();
    std::vector<unsigned char> in_queue( n, 1u );
    std::vector<int> queue;
    queue.reserve( n );
    for( int node_id : nodes_ordered_randomly ) {
        queue.push_back( node_id );
    }
    std::size_t head = 0;

    while( head < queue.size() ) {
        int node_id = queue[head++];
        in_queue[node_id] = 0u;

        unsigned int comm_id = partition.find_set( node_id );

        isolate_and_update_internals( node_id, internals, partition );
        BestMove move = find_best_comm_move( compute_quality_diff,
                                             static_cast<unsigned int>( node_id ),
                                             internals, comm_id );
        insert_and_update_internals( node_id, internals, partition, move.best_comm );

        if( move.best_comm != comm_id ) {
            do_construct_new_graph = true;
            // Only re-enqueue neighbours when the move's improvement clears
            // the noise floor, mirroring the original `> minimum_improve`
            // termination criterion (now applied per move rather than per
            // sweep).
            if( move.net_delta > minimum_improve ) {
                const std::size_t beg = internals.nbr.start[node_id];
                const std::size_t end = internals.nbr.start[node_id + 1];
                for( std::size_t e = beg; e < end; ++e ) {
                    int nbr_id = internals.nbr.node[e];
                    if( !in_queue[nbr_id] ) {
                        in_queue[nbr_id] = 1u;
                        queue.push_back( nbr_id );
                    }
                }
            }
        }
    }

    ////////////////////////////////////////////////////////////
    // Start Second phase - create reduced graph with self loops
    ////////////////////////////////////////////////////////////

    // 1) Normalise partition IDs and store next level in optimal partitions found by Louvain
    std::vector<int> new_comm_id_to_old_comm_id = partition.normalise_ids();

    // 2) If there has actually been some movement, then we need to assemble a new graph
    if( do_construct_new_graph == true ) {
        P partition_original_nodes = fold_partition_into_orginal_graph_size( optimal_partitions, partition );
        optimal_partitions.push_back( partition_original_nodes );

        Graph reduced_graph;
        create_reduced_graph_from_partition( reduced_graph, graph, partition,
                                             new_comm_id_to_old_comm_id, internals );

        unsigned int num_nodes_reduced_graph = reduced_graph.num_nodes();
        P reduced_partition( num_nodes_reduced_graph );
        reduced_partition.initialise_as_singletons();
        auto reduced_null_model_vec = create_reduced_null_model_vec( new_comm_id_to_old_comm_id, internals,
                                      static_cast<unsigned int>( null_model_vec.size() ),
                                      num_nodes_reduced_graph );

        return find_optimal_partition_louvain_gen( reduced_graph, reduced_null_model_vec,
                compute_quality, compute_quality_diff, reduced_partition, optimal_partitions, minimum_improve, rng );
    } else {
        if( optimal_partitions.empty() ) {
            optimal_partitions.push_back( partition );
        }
        // Final quality is computed exactly once, at the bottom of the
        // recursion — the per-sweep recompute that the original code did is
        // unnecessary because the stopping condition only needs the per-move
        // delta, which we track inline.
        return compute_quality( internals );
    }
}


/**
 @brief  standard Louvain method - greedy algorithm to find community structure of a network.
 In contrast to the generalised variant no null model vectors are provided but they are 
 created directly from the data.

 @param[in]  graph    --    weighted graph to partition
 @param[in]  compute_quality -- functor to compute the quality function
 @param[in]  compute_quality_diff -- functor to compute change in quality function
 @param[in]  initial_partition -- partition to start from
 @param[in]  optimal_partitions -- set of optimal_partitions at each level of the algorithm (initially empty)
 @param[in]  minimum_improve -- minimum_improvement necessary in each iteration / stopping criterion

 @param[out]  (double) quality of the optimal partition

 */
template<typename P, typename QF, typename QFDIFF>
double find_optimal_partition_louvain( const Graph& graph,
        QF compute_quality, QFDIFF compute_quality_diff, P initial_partition,
        std::vector<P>& optimal_partitions, double minimum_improve,
        std::mt19937& rng )
{
    P partition( initial_partition );
    P partition_init( initial_partition );

    bool do_construct_new_graph = false;

    if( !optimal_partitions.empty() ) {
        partition_init = optimal_partitions.back();
    }

    auto internals = clq::gen_internals( compute_quality, graph, partition, partition_init );

    std::vector<int> nodes_ordered_randomly;
    nodes_ordered_randomly.reserve( static_cast<std::size_t>( graph.num_nodes() ) );

    for( int node_id = graph.num_nodes() - 1; node_id >= 0; --node_id ) {
        nodes_ordered_randomly.push_back( node_id );
    }

    for( std::size_t i = nodes_ordered_randomly.size(); i > 1; --i ) {
        std::size_t j = rng() % i;
        std::swap( nodes_ordered_randomly[j], nodes_ordered_randomly[i - 1] );
    }

    const std::size_t n = nodes_ordered_randomly.size();
    std::vector<unsigned char> in_queue( n, 1u );
    std::vector<int> queue;
    queue.reserve( n );
    for( int node_id : nodes_ordered_randomly ) {
        queue.push_back( node_id );
    }
    std::size_t head = 0;

    while( head < queue.size() ) {
        int node_id = queue[head++];
        in_queue[node_id] = 0u;

        unsigned int comm_id = partition.find_set( node_id );

        isolate_and_update_internals( node_id, internals, partition );
        BestMove move = find_best_comm_move( compute_quality_diff,
                                             static_cast<unsigned int>( node_id ),
                                             internals, comm_id );
        insert_and_update_internals( node_id, internals, partition, move.best_comm );

        if( move.best_comm != comm_id ) {
            do_construct_new_graph = true;
            if( move.net_delta > minimum_improve ) {
                const std::size_t beg = internals.nbr.start[node_id];
                const std::size_t end = internals.nbr.start[node_id + 1];
                for( std::size_t e = beg; e < end; ++e ) {
                    int nbr_id = internals.nbr.node[e];
                    if( !in_queue[nbr_id] ) {
                        in_queue[nbr_id] = 1u;
                        queue.push_back( nbr_id );
                    }
                }
            }
        }
    }

    ////////////////////////////////////////////////////////////
    // Start Second phase - create reduced graph with self loops
    ////////////////////////////////////////////////////////////

    std::vector<int> new_comm_id_to_old_comm_id = partition.normalise_ids();

    if( do_construct_new_graph == true ) {
        P partition_original_nodes = fold_partition_into_orginal_graph_size( optimal_partitions, partition );
        optimal_partitions.push_back( partition_original_nodes );

        Graph reduced_graph;
        create_reduced_graph_from_partition( reduced_graph, graph, partition,
                                             new_comm_id_to_old_comm_id, internals );

        unsigned int num_nodes_reduced_graph = reduced_graph.num_nodes();
        P reduced_partition( num_nodes_reduced_graph );
        reduced_partition.initialise_as_singletons();

        return find_optimal_partition_louvain( reduced_graph, compute_quality,
                compute_quality_diff, reduced_partition, optimal_partitions, minimum_improve, rng );
    } else {
        if( optimal_partitions.empty() ) {
            optimal_partitions.push_back( partition );
        }
        return compute_quality( internals );
    }
}


// Backwards-compatible overloads: if no RNG is supplied, fall back to a
// thread-local std::mt19937 seeded from std::random_device. This preserves
// the old call signature for existing callers.
template<typename P, typename QF, typename QFDIFF>
double find_optimal_partition_louvain_gen( const Graph& graph, const vec2& null_model_vec,
        QF compute_quality, QFDIFF compute_quality_diff, P initial_partition,
        std::vector<P>& optimal_partitions, double minimum_improve )
{
    static thread_local std::mt19937 default_rng{ std::random_device{}() };
    return find_optimal_partition_louvain_gen( graph, null_model_vec,
            compute_quality, compute_quality_diff, initial_partition,
            optimal_partitions, minimum_improve, default_rng );
}

template<typename P, typename QF, typename QFDIFF>
double find_optimal_partition_louvain( const Graph& graph,
        QF compute_quality, QFDIFF compute_quality_diff, P initial_partition,
        std::vector<P>& optimal_partitions, double minimum_improve )
{
    static thread_local std::mt19937 default_rng{ std::random_device{}() };
    return find_optimal_partition_louvain( graph,
            compute_quality, compute_quality_diff, initial_partition,
            optimal_partitions, minimum_improve, default_rng );
}


}// end namespace..
