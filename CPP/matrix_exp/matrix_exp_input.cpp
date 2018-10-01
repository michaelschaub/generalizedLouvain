# include <vector>
# include <cmath>
# include "../cliques/math.h"
# include "../cliques/io.h"
# include <sstream>

// function to compute the matrix exponential given an input adjacency list
// output is stored as matrix in a file named EXP_TIME_INPUTFILENAME.dat
int main(int argc, char *argv []) {

    // create adjacency matrix in vectorised form
    std::vector<double> A = clq::read_edgelist_weighted(argv[1]);
    int size = std::sqrt(A.size());

    double input_time;
    if (argc < 3)
    {
        input_time = 1.0;
    } else {
        input_time = std::atof(argv[2]);
    }
    
    // assemble output filename
    std::string suffix(argv[3]);
    std::string outfilename("EXP_");
    std::ostringstream convert;
    convert << input_time;
    outfilename += convert.str() + "_" + suffix + ".dat";
    std::cout << outfilename;

    std::vector<double> result = clq::exp(A,input_time,size);
    clq::write_adj_matrix(outfilename,result);
    return 0;
}
