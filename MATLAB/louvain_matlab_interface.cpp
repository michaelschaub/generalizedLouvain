#include <lemon/smart_graph.h>
#include <map>
#include "../CPP/cliques/louvain_gen.h"
#include "../CPP/cliques/stability.h"
#include "../CPP/cliques/stability_gen.h"
#include "../CPP/cliques/vector_partition.h"
#include <vector>
#include <string>
#include <stdlib.h>

// some redefinitions needed here fixing a bug in case MATLAB tries 
// to define its own symbols...
//namespace matlab {
//#define CHAR16_T UINT16_T // VALID!?! be careful!!
#include "mex.h"
#include "matrix.h"
//}


// GLOBAL DATA
double *data = NULL;
double precision = 1e-15; // default value here
int mode = 0;
// 0 = normalised Laplacian; 1 = combinatorial Laplacian; etc...

double m_time;
std::vector<std::vector<double>> null_vectors;
int num_iterations = 0;
unsigned int seed = -1; // random seed (to be set by call)

bool hierarchy = false;
int num_largest_dim = -1;

// parse function definition, needed by constructor
bool parse_arg(int nrhs, const mxArray *prhs[]) {

    //FIRST ARGUMENT: Graph
    if (nrhs > 0) {
        // number of columns should be 3 (n1,n2, weight)
        if (mxGetN(prhs[0]) != 3) {
            mexPrintf("Number of columns %d", mxGetN(prhs[0]));
            return false;
        }
        // data is stored in column major ordering ordering
        data = (double *) mxGetPr(prhs[0]);
        // number of rows is important for indexing
        num_largest_dim = mxGetM(prhs[0]);
    }

    //SECOND ARGUMENT: time
    if (nrhs > 1) {
        m_time = double(mxGetScalar(prhs[1]));
        //mexPrintf("input time %f",m_time);
        if (m_time < 0){
            return false;
        }
    }

    //THIRD ARGUMENT: iterations
    if (nrhs > 2) {
        num_iterations = int((double) mxGetScalar(prhs[2]));
        //std::cout << "number of iterations" << num_iterations << std::endl;
        if (num_iterations < 1) {
            return false;
        }
    }

    //FOURTH ARGUMENT: precision
    if (nrhs > 3) {
        precision = ((double) mxGetScalar(prhs[3]));
        if (precision > 1) {
            return false;
        }
    }
    //FIFTH ARGUMENT
    if (nrhs > 4) {
        // get buffer length and allocate buffer
        char *buf;
        // read out argument, take care to use correct pointers!
        mwSize buflen = mxGetN(prhs[4]) * sizeof(mxChar) + 1;
        buf = (char*) mxMalloc(buflen);
        if (!mxGetString(prhs[4], buf, buflen)) {
            //if reading successful string is created from matlab input
            const std::string input(buf);
            const std::string comparison1("normalised");
            const std::string comparison2("combinatorial");
            const std::string comparison3("generalised");

            // in case valid input change mode, else display error message
            if (!comparison1.compare(input)) {
                mode = 0;
            } else if (!comparison2.compare(input)) {
                mode = 1;
            } else if (!comparison3.compare(input)) {
                mode = 2;
            } else {
                mexErrMsgTxt("No valid stability mode specified \n");
                return false;
            }
        }
        // free read buffer
        mxFree(buf);
    }

    //SIXTH ARGUMENT: random seed passed by Matlab
    if (nrhs > 5) {
        seed = ((unsigned int) mxGetScalar(prhs[5]));
        if (seed <=0) {
            return false;
        }
    }
    
    // SEVENTH ARGUMENT -- possible null model vec input
    if (nrhs > 6){
        int num_null_model_vectors = mxGetN(prhs[6]);
        if (num_null_model_vectors % 2 != 0){
            return false;
        }
        unsigned int num_null_vec_entries = mxGetM(prhs[6]);
        null_vectors.clear();
        null_vectors = std::vector<std::vector<double>>(num_null_model_vectors,std::vector<double>(num_null_vec_entries,0));
        for (int i= 0; i < num_null_model_vectors; i++){
            for(unsigned int j = 0; j < num_null_vec_entries; j++){
                // data comes in column first order from matlab
                double * entries = (double *) mxGetPr(prhs[6]);
                null_vectors[i][j] = entries[j + num_null_vec_entries*i];

            } 
        }
    }

    //this is not implemented fully/has no effect so far.
    //SEVENTH ARGUMENT: hierarchical output
    if (nrhs > 7) {
        // get buffer length and allocate buffer
        char *buf;
        mwSize buflen = mxGetN(prhs[6]) * sizeof(mxChar) + 1;
        buf = (char*) mxMalloc(buflen);
        if (!mxGetString(prhs[6], buf, buflen)) {
            //if reading successful string is created from matlab input
            const std::string input(buf);
            const std::string comparison("h");

            // in case hierarchical output is activated print message
            if (!comparison.compare(input)) {
                hierarchy = true;
                mexPrintf("Hierarchical output from Louvain activated \n");
            }
        }
        // free read buffer
        mxFree(buf);
    }

    // SANITY CHECK for number of arguments
    if (nrhs > 8 || nrhs < 1) {
        return false;
    }
    return true;
}

// Template for reading in graph from weighted edgelist data as coming from Matlab
template<typename G, typename E>
bool read_edgelist_weighted_from_data(int num_l_dim, G &graph, E &weights) {

    // define Node class for convenience
    typedef typename G::Node Node;

    int max_node_id_seen =-1;

    // loop over complete list
    for (int i = 0; i < num_l_dim; ++i) {

        // get nodes and weights
        // column major ordering from MATLAB
        int node1_id = data[i];
        int node2_id = data[num_l_dim + i];

        if (node1_id > max_node_id_seen) {
            int difference = node1_id - max_node_id_seen;
            for (int i=0; i<difference; ++i) {
                graph.addNode();
            }
            max_node_id_seen = node1_id;
        }

        if (node2_id > max_node_id_seen) {
            int difference = node2_id - max_node_id_seen;
            for (int i=0; i<difference; ++i) {
                graph.addNode();
            }
            max_node_id_seen = node2_id;
        }

        // read in list is two-way yet undirected, but edges should only be created once
        if (node1_id > node2_id) {
            std::cout << "you provided a double sided list";
            std::cout << "this is not necessary -- some entries will be ignored";
        }

        double weight = data[2 * num_l_dim + i];
        typename G::Edge edge = graph.addEdge(graph.nodeFromId(node1_id),
                                              graph.nodeFromId(node2_id));
        weights.set(edge, weight);
    }

    return true;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // Parse arguments and return 0 if there is an error
    if (!parse_arg(nrhs, prhs)) {
        // give back zero for all outputs
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        if (nlhs > 1)
            plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
        if (nlhs > 2)
            plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
        mexErrMsgTxt("Error parsing arguments");
    }
    //create new graph and weight map
    lemon::SmartGraph mygraph;
    lemon::SmartGraph::EdgeMap<double> myweights(mygraph);

    if (!read_edgelist_weighted_from_data(num_largest_dim, mygraph, myweights)) {
        mexErrMsgTxt("Error creating graph from data");
    }

    // typedef for convenience
    typedef clq::VectorPartition partition;

    // get number of nodes 
    unsigned int num_nodes = lemon::countNodes(mygraph);

    // vector of initial partition
    partition singletons(num_nodes);
    singletons.initialise_as_singletons();

    // vector of vectors for output partitions
    std::vector<partition> optimal_partitions;
    std::vector<double> null_model(num_nodes, 0);
    //initialise stabilities
    std::vector<double> stability(num_iterations, 0);

    // randomize initial seed
    srand(seed);

    for (int i = 0; i < num_iterations; ++i) {

            // create empty vector of partitions
            std::vector<partition> hierarchical_louvain_partitions;

            // initialise quality functions
            switch (mode) {

            // normalised 
            case 0: {
                clq::find_linearised_normalised_stability quality(
                    m_time);
                clq::linearised_normalised_stability_gain quality_gain(
                    m_time);
                // now run Louvain method
                stability[i]
                = clq::find_optimal_partition_louvain<partition>(
                      mygraph, myweights, quality,
                      quality_gain, singletons,
                      hierarchical_louvain_partitions, precision);
                break;
            }

            // combinatorial
            case 1: {
                clq::find_linearised_combinatorial_stability quality(
                    m_time);
                clq::linearised_combinatorial_stability_gain quality_gain(
                    m_time);
                // now run Louvain method
                stability[i]
                = clq::find_optimal_partition_louvain<partition>(
                      mygraph, myweights, quality,
                      quality_gain, singletons,
                      hierarchical_louvain_partitions, precision);
                break;
            }

            // generalised 
            case 2: {
                clq::find_linearised_generic_stability quality(m_time);
                clq::linearised_generic_stability_gain quality_gain(m_time);
                // now run Louvain method
                stability[i]
                = clq::find_optimal_partition_louvain_gen<partition>(
                      mygraph, myweights, null_vectors, quality,
                      quality_gain, singletons,
                      hierarchical_louvain_partitions, precision);
                break;
            }
            default:
                mexErrMsgTxt("Error defining stability mode");

            }

            // last partition in vector == best partition
            partition best_partition = hierarchical_louvain_partitions.back();

            // store best partition
            optimal_partitions.push_back(best_partition);
        
    }

    //****************************************************
    //----------------------------------------------------
    // Now write data back to Matlab
    //----------------------------------------------------

    // column and row dimensions
    mwSize columns = num_iterations;
    mwSize rows = num_nodes;

    /////////////////////////////////////////
    // FIRST output: stability

    // mxReal is our data-type
    plhs[0] = mxCreateDoubleMatrix(1, columns, mxREAL);
    //Get a pointer to the data space in our newly allocated memory
    double * out1 = (double*) mxGetPr(plhs[0]);

    // write out stabilities
    for (int k = 0; k < columns; ++k) {
        out1[k] = double(stability[k]);
    }

    ////////////////////////////////////////
    // SECOND output: number of communities
    if (nlhs > 1) {

        plhs[1] = mxCreateDoubleMatrix(1, columns, mxREAL); //mxReal is our data-type

        //Get a pointer to the data space in our newly allocated memory
        double * out2 = (double*) mxGetPr(plhs[1]);

        // write out number of communities
        for (int k = 0; k < columns; ++k) {
            out2[k] = double(optimal_partitions[k].set_count());
        }

    }
    ////////////////////////////////////////
    // THIRD output: community assignments
    if (nlhs > 2) {

        // allocate storage
        plhs[2] = mxCreateDoubleMatrix(rows, columns, mxREAL);
        double * output_tab = (double*) mxGetPr(plhs[2]);

        // write out results
        for (int iteration = 0; iteration < columns; ++iteration) {
            for (int node = 0; node < rows; ++node) {
                output_tab[iteration * num_nodes + node]
                = double(optimal_partitions[iteration].find_set(node));
            }
        }
    }

}
