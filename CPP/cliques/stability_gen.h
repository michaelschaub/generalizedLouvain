#pragma once

#include <vector>
#include "compensated_sum.h"
#include "graphhelpers.h"
#include "linearised_internals_generic.h"

namespace clq
{

/**
 @brief  A functor for evaluating the generic stability of a weighted graph
 */

struct find_linearised_generic_stability {
    // Markov time state variable
    double markov_time;

    // Constructor with default Markov time of 1
    find_linearised_generic_stability (double markov_time = 1.0) :
        markov_time (markov_time)
    {
    }

    template<typename I>
    double operator () (I& internals)
    {
        // first part equals 1-t
        double q = 1 - markov_time;
        // Compensated sum: the error does not increase with the number of terms.
        CompensatedSum q2;

        const unsigned int num_null_model_vec = internals.num_null_model_vectors;
        const unsigned int num_nodes = internals.num_nodes;

        // Loop over all possible community indices. Loss-vector entries for
        // empty communities are zero, so they contribute zero to the sum —
        // walking past them costs only memory bandwidth. (The `if (gain != 0)
        // { ... } else { ... }` dead branch in the previous implementation
        // was bug-compatible with always evaluating the inner sum.)
        for (unsigned int i = 0; i < num_nodes; i++) {
            q2.add (markov_time * internals.comm_w_in[i]);
            for (unsigned int j = 0; j < num_null_model_vec; j += 2) {
                q2.add (- internals.comm_loss_vectors[j][i]
                        * internals.comm_loss_vectors[j + 1][i]);
            }
        }

        return q + q2.value();
    }

};


/**
 @brief  Functor for finding stability gain (normalised Laplacian) with for weighted graph
 */
struct linearised_generic_stability_gain {

    //state variable Markov time
    double markov_time;

    //constructor
    linearised_generic_stability_gain (double mtime = 1.0) :
        markov_time (mtime)
    {
    }

    template<typename I>
    double operator () (I& internals, int comm_id_neighbour, int node_id)
    {
        // compute loss factor incurred by moving node...
        double comm_loss = 0;

        //clq::output(comm_loss,internals.num_null_model_vectors,"loss\n");
        for (unsigned int i = 0; i < internals.num_null_model_vectors; i = i + 2) {
            //clq::output("DEBUG", i, node_id, comm_id_neighbour);
            comm_loss +=
                internals.null_model_vectors[i][node_id] * internals.comm_loss_vectors[i + 1][comm_id_neighbour]
                + internals.null_model_vectors[i + 1][node_id] * internals.comm_loss_vectors[i][comm_id_neighbour];
        }

        //clq::output(comm_loss,internals.num_null_model_vectors,"loss\n");

        // gain resulting from adding to community
        double w_node_to_comm = internals.node_weight_to_communities[comm_id_neighbour];
        return markov_time * w_node_to_comm * 2 - comm_loss;
    }
};

}
