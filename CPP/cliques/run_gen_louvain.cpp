#include "stability_gen.h"
#include "louvain_gen.h"
#include "vector_partition.h"
#include "io.h"

#include <ctime>
#include <random>
#include <vector>

// Call generic Louvain optimisation from command line
// Inputs: graph, null vectors (both as txt files), resolution parameter (time) and random seed (both optional)
// Example: ./run_gen_louvain.sh graph_file null_model_file time random_seed
int main(int argc, char *argv []) {

    // initialise random number generator from input
    unsigned int seed = (argc < 5) ? static_cast<unsigned int>(std::time(nullptr))
                                   : static_cast<unsigned int>(std::atoi(argv[4]));
    std::mt19937 rng(seed);

    // initialise markov time from input
    double current_markov_time = 1;
    if (argc > 3 ){
        current_markov_time = atof(argv[3]);
    }
    clq::output("Markov time: ", current_markov_time);

    // initialise graph
    clq::Graph input_graph;

    
    clq::read_edgelist_weighted_graph(argv[1], input_graph);
    
    int num_nodes = input_graph.num_nodes();
    std::vector<std::vector<double>> null_model = clq::read_null_model(argv[2],num_nodes);
    clq::print_2d_vector(null_model);

    
    // define type for vector partition initialise start partition, qualitiy functions etc.
    typedef clq::VectorPartition partition;

    partition start_partition(num_nodes);
    start_partition.initialise_as_singletons();
    clq::print_partition_list(start_partition);


    clq::find_linearised_generic_stability quality(current_markov_time);
    clq::linearised_generic_stability_gain quality_gain(current_markov_time);
    
    std::vector<partition> optimal_partitions;

    // call Louvain
    clq::output("Start Louvain");
    double stability = clq::find_optimal_partition_louvain_gen<partition>(
                           input_graph, null_model, quality, quality_gain,
                           start_partition, optimal_partitions, 1e-18, rng);
    
    // store output partitions in file
    clq::partitions_to_file("optimal_partitions.dat", optimal_partitions);

    // display final result
    partition best_partition = optimal_partitions.back();
    clq::print_partition_list(best_partition);
    clq::output("Stability: ", stability);

    return 0;
}
