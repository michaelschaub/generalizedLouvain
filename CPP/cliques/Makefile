# Compiler — override on the command line, e.g. `make CXX=g++-13`.
CXX ?= g++

# Compilation flags
#  -O3              : standard high optimisation
#  -DNDEBUG         : remove asserts in hot paths
#  -march=native    : let the compiler use SIMD/cache-line widths of the host
#  -funroll-loops   : helpful for the small inner null-model loops
#  -flto            : whole-program inlining of header-only templates
#  -ffp-contract=off: do not merge a*b+c into one FMA instruction. Clang does
#                     this by default and GCC does not. Without this flag, the
#                     same seed can give different results with each compiler.
#  -std=c++17       : the code also compiles with C++11 (PyGenStability uses it)
# Note: -ffast-math is intentionally NOT enabled — modularity is a sum of
# products and aggressive re-association can perturb the partition.
CXXFLAGS ?= -O3 -DNDEBUG -march=native -funroll-loops -flto -ffp-contract=off -std=c++17 -Wall

.PHONY: all run_gen_louvain test clean

all: run_gen_louvain

run_gen_louvain:
	$(CXX) $(CXXFLAGS) run_gen_louvain.cpp -o run_gen_louvain.sh

# Build and run the tests. The test data is in this directory.
test:
	$(CXX) $(CXXFLAGS) test/test_louvain.cpp -o test/test_louvain
	./test/test_louvain .

clean:
	rm -f run_gen_louvain.sh test/test_louvain
