// Tests for the Louvain optimisation and the quality functions.
//
// Build and run from CPP/cliques with `make test`.
// Use `./test/test_louvain --print` to print the reference values of the
// current code, for example after an intended change of results.
//
// The test graphs come from std::mt19937 raw output only. The standard fixes
// this sequence; it does not fix the std:: distributions.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "../compensated_sum.h"
#include "../graph.h"
#include "../io.h"
#include "../louvain_gen.h"
#include "../stability.h"
#include "../stability_gen.h"
#include "../vector_partition.h"

namespace {

typedef clq::VectorPartition Partition;
typedef std::vector<std::vector<double>> vec2;

int num_checks = 0;
int num_failures = 0;
bool print_mode = false;

#define CHECK(cond)                                                          \
    do {                                                                     \
        ++num_checks;                                                        \
        if (!(cond)) {                                                       \
            ++num_failures;                                                  \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);      \
        }                                                                    \
    } while (0)

bool near(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

struct Edge {
    int u;
    int v;
    double w;
};

struct TestGraph {
    int num_nodes;
    std::vector<Edge> edges;
};

clq::Graph make_graph(const TestGraph &g)
{
    clq::Graph graph(g.num_nodes);
    for (const Edge &e : g.edges) {
        graph.add_edge(e.u, e.v, e.w);
    }
    return graph;
}

double run_generic(const TestGraph &g, const vec2 &null_model, double time,
                   unsigned int seed, std::vector<Partition> &levels)
{
    clq::Graph graph = make_graph(g);
    Partition start(g.num_nodes);
    start.initialise_as_singletons();
    std::mt19937 rng(seed);
    clq::find_linearised_generic_stability quality(time);
    clq::linearised_generic_stability_gain gain(time);
    return clq::find_optimal_partition_louvain_gen<Partition>(
               graph, null_model, quality, gain, start, levels, 1e-18, rng);
}

template<typename QF, typename QFDIFF>
double run_standard(const TestGraph &g, double time, unsigned int seed,
                    std::vector<Partition> &levels)
{
    clq::Graph graph = make_graph(g);
    Partition start(g.num_nodes);
    start.initialise_as_singletons();
    std::mt19937 rng(seed);
    QF quality(time);
    QFDIFF gain(time);
    return clq::find_optimal_partition_louvain<Partition>(
               graph, quality, gain, start, levels, 1e-18, rng);
}

// Generic quality of a partition, from the internals of the library.
double library_generic_quality(const TestGraph &g, const vec2 &null_model, double time,
                               Partition partition)
{
    clq::Graph graph = make_graph(g);
    clq::LinearisedInternalsGeneric internals(graph, partition, null_model);
    clq::find_linearised_generic_stability quality(time);
    return quality(internals);
}

double library_normalised_quality(const TestGraph &g, double time, Partition partition)
{
    clq::Graph graph = make_graph(g);
    clq::find_linearised_normalised_stability quality(time);
    return quality(graph, partition);
}

double library_combinatorial_quality(const TestGraph &g, double time, Partition partition)
{
    clq::Graph graph = make_graph(g);
    Partition partition_init(g.num_nodes);
    partition_init.initialise_as_singletons();
    clq::find_linearised_combinatorial_stability quality(time);
    return quality(graph, partition, partition_init);
}

// Generic quality of a partition, calculated independently of the library:
// Q = 1 - t + t * sum_{i,j in same community} A_ij - sum_c sum_k a_k(c) b_k(c).
double reference_generic_quality(const TestGraph &g, const vec2 &null_model, double time,
                                 const Partition &partition)
{
    double internal_weight = 0.0;
    for (const Edge &e : g.edges) {
        if (partition.find_set(e.u) == partition.find_set(e.v)) {
            internal_weight += (e.u == e.v) ? e.w : 2.0 * e.w;
        }
    }
    double loss = 0.0;
    for (std::size_t k = 0; k + 1 < null_model.size(); k += 2) {
        std::vector<double> a(g.num_nodes, 0.0), b(g.num_nodes, 0.0);
        for (int i = 0; i < g.num_nodes; ++i) {
            a[partition.find_set(i)] += null_model[k][i];
            b[partition.find_set(i)] += null_model[k + 1][i];
        }
        for (int c = 0; c < g.num_nodes; ++c) {
            loss += a[c] * b[c];
        }
    }
    return 1.0 - time + time * internal_weight - loss;
}

// Planted-partition graph with integer weights 1 to 3.
TestGraph planted_partition(int num_blocks, int block_size, double p_in, double p_out,
                            unsigned int seed)
{
    std::mt19937 rng(seed);
    TestGraph g;
    g.num_nodes = num_blocks * block_size;
    for (int i = 0; i < g.num_nodes; ++i) {
        for (int j = i + 1; j < g.num_nodes; ++j) {
            double p = (i / block_size == j / block_size) ? p_in : p_out;
            double r = rng() / 4294967296.0;
            unsigned int w = rng() % 3u;
            if (r < p) {
                g.edges.push_back(Edge{i, j, 1.0 + w});
            }
        }
    }
    return g;
}

std::vector<double> degrees(const TestGraph &g)
{
    std::vector<double> d(g.num_nodes, 0.0);
    for (const Edge &e : g.edges) {
        d[e.u] += e.w;
        d[e.v] += e.w;
    }
    return d;
}

double total_degree(const TestGraph &g)
{
    double two_m = 0.0;
    for (const Edge &e : g.edges) {
        two_m += 2.0 * e.w;
    }
    return two_m;
}

// Same graph, with all weights divided by 2m (the input of the generic quality).
TestGraph normalised_weights(const TestGraph &g)
{
    TestGraph h = g;
    double two_m = total_degree(g);
    for (Edge &e : h.edges) {
        e.w /= two_m;
    }
    return h;
}

// Null model of modularity: one pair (d / 2m, d / 2m).
vec2 modularity_null_model(const TestGraph &g)
{
    std::vector<double> d = degrees(g);
    double two_m = total_degree(g);
    for (double &x : d) {
        x /= two_m;
    }
    return vec2{d, d};
}

// Modularity null model plus an asymmetric pair (a, b), with a != b.
vec2 two_pair_null_model(const TestGraph &g, unsigned int seed)
{
    vec2 null_model = modularity_null_model(g);
    std::mt19937 rng(seed);
    std::vector<double> a(g.num_nodes), b(g.num_nodes);
    for (int i = 0; i < g.num_nodes; ++i) {
        a[i] = 0.2 * null_model[0][i] * (1.0 + rng() % 5u);
        b[i] = 0.1 * null_model[0][i] * (1.0 + rng() % 7u);
    }
    null_model.push_back(a);
    null_model.push_back(b);
    return null_model;
}

std::uint64_t partition_hash(const Partition &partition)
{
    std::uint64_t h = 1469598103934665603ull;  // FNV-1a
    for (int i = 0; i < partition.element_count(); ++i) {
        h ^= static_cast<std::uint64_t>(partition.find_set(i) + 1);
        h *= 1099511628211ull;
    }
    return h;
}

int community_count(const Partition &partition)
{
    Partition p(partition);
    return p.set_count();
}

// Reference results. Regenerate with `./test/test_louvain --print`.
struct Reference {
    const char *name;
    double quality;
    std::uint64_t hash;
    int num_communities;
};

const Reference references[] = {
    {"generic_modularity_seed0", 0.4462074687228928, 10928742383526330306ull, 6},
    {"generic_modularity_seed1", 0.446964699807465, 11147939978088631755ull, 6},
    {"generic_modularity_seed2", 0.46404726633281324, 6784018557408576015ull, 5},
    {"generic_modularity_seed3", 0.45512264467285568, 2932209541537649572ull, 5},
    {"generic_modularity_t0.3", 0.76241237136123696, 5848503960239146080ull, 20},
    {"generic_two_pairs_t0.7", 0.53002223677561244, 7949468812565776579ull, 10},
    {"generic_large_t1", 0.47503246840150387, 18396511715021903318ull, 10},
    {"normalised_t1", 0.45380260670110151, 9107574923049613561ull, 5},
    {"normalised_large_t0.5", 0.68641706740215203, 3909507305208399923ull, 15},
    {"combinatorial_t0.5", 0.65806925972396491, 14161383697746804514ull, 13},
    {"triangle_files_seed42", 4.9999999985999999, 215615651957586990ull, 2},
};

const Reference *find_reference(const std::string &name)
{
    for (const Reference &r : references) {
        if (name == r.name) {
            return &r;
        }
    }
    return nullptr;
}

void check_reference(const std::string &name, double quality, const Partition &partition)
{
    std::uint64_t hash = partition_hash(partition);
    int num_communities = community_count(partition);
    if (print_mode) {
        std::printf("    {\"%s\", %.17g, %lluull, %d},\n", name.c_str(), quality,
                    static_cast<unsigned long long>(hash), num_communities);
        return;
    }
    const Reference *r = find_reference(name);
    CHECK(r != nullptr);
    if (r == nullptr) {
        std::printf("  no reference for %s\n", name.c_str());
        return;
    }
    // The partition must be identical. The quality can change by a few ULP
    // when the summation order changes.
    CHECK(r->hash == hash);
    CHECK(r->num_communities == num_communities);
    CHECK(near(quality, r->quality, 1e-12));
    if (r->hash != hash || !near(quality, r->quality, 1e-12)) {
        std::printf("  %s: got Q = %.17g, %d communities\n", name.c_str(), quality,
                    num_communities);
    }
}

void check_levels(const std::vector<Partition> &levels, int num_nodes)
{
    CHECK(!levels.empty());
    for (const Partition &p : levels) {
        CHECK(p.element_count() == num_nodes);
    }
    for (std::size_t k = 1; k < levels.size(); ++k) {
        CHECK(community_count(levels[k]) <= community_count(levels[k - 1]));
    }
}

void test_generic(const std::string &name, const TestGraph &g, const vec2 &null_model,
                  double time, unsigned int seed)
{
    std::vector<Partition> levels;
    double q = run_generic(g, null_model, time, seed, levels);
    check_levels(levels, g.num_nodes);
    const Partition &best = levels.back();

    // The quality that Louvain gives agrees with the two other calculations.
    CHECK(near(q, library_generic_quality(g, null_model, time, best), 1e-12));
    CHECK(near(q, reference_generic_quality(g, null_model, time, best), 1e-12));

    // The same seed gives the same partition.
    std::vector<Partition> levels_again;
    double q_again = run_generic(g, null_model, time, seed, levels_again);
    CHECK(q_again == q);
    CHECK(levels_again.back().return_partition_vector() ==
          Partition(best).return_partition_vector());

    check_reference(name, q, best);
}

template<typename QF, typename QFDIFF>
void test_standard(const std::string &name, const TestGraph &g, double time,
                   unsigned int seed, bool normalised)
{
    std::vector<Partition> levels;
    double q = run_standard<QF, QFDIFF>(g, time, seed, levels);
    check_levels(levels, g.num_nodes);
    const Partition &best = levels.back();

    if (normalised) {
        CHECK(near(q, library_normalised_quality(g, time, best), 1e-12));
        // Normalised linearised stability is the generic quality with weights A / 2m
        // and the null model (d / 2m, d / 2m).
        CHECK(near(q, reference_generic_quality(normalised_weights(g), modularity_null_model(g),
                                                time, best),
                   1e-12));
    } else {
        CHECK(near(q, library_combinatorial_quality(g, time, best), 1e-12));
        // Combinatorial: null model (1 / n, 1 / n).
        std::vector<double> u(g.num_nodes, 1.0 / g.num_nodes);
        CHECK(near(q, reference_generic_quality(normalised_weights(g), vec2{u, u}, time, best),
                   1e-12));
    }
    check_reference(name, q, best);
}

void test_planted_partition_is_found()
{
    TestGraph g = planted_partition(4, 30, 0.5, 0.01, 7u);
    std::vector<Partition> levels;
    run_generic(normalised_weights(g), modularity_null_model(g), 1.0, 1u, levels);
    const Partition &best = levels.back();
    CHECK(community_count(best) == 4);
    for (int i = 0; i < g.num_nodes; ++i) {
        CHECK(best.find_set(i) == best.find_set((i / 30) * 30));
    }
}

void test_triangle_files(const std::string &dir)
{
    clq::Graph graph;
    clq::read_edgelist_weighted_graph(dir + "/triangletest.edj", graph);
    int num_nodes = graph.num_nodes();
    CHECK(num_nodes == 6);
    vec2 null_model = clq::read_null_model(dir + "/triangletest_null.edj", num_nodes);
    CHECK(null_model.size() == 2);

    Partition start(num_nodes);
    start.initialise_as_singletons();
    std::vector<Partition> levels;
    std::mt19937 rng(42u);
    clq::find_linearised_generic_stability quality(1.0);
    clq::linearised_generic_stability_gain gain(1.0);
    double q = clq::find_optimal_partition_louvain_gen<Partition>(
                   graph, null_model, quality, gain, start, levels, 1e-18, rng);
    check_levels(levels, num_nodes);
    check_reference("triangle_files_seed42", q, levels.back());
}

// The compensated sum gives the correctly rounded sum where a simple sum does not.
void test_compensated_sum()
{
    clq::CompensatedSum sum;
    double simple_sum = 0.0;
    for (int i = 0; i < 10000000; ++i) {
        sum.add(0.1);
        simple_sum += 0.1;
    }
    CHECK(sum.value() == 1e6);
    CHECK(simple_sum != 1e6);

    // A term larger than the sum: plain Kahan summation gives 0 here.
    clq::CompensatedSum cancel;
    cancel.add(1.0);
    cancel.add(1e100);
    cancel.add(1.0);
    cancel.add(-1e100);
    CHECK(cancel.value() == 2.0);
}

// Minimal internals for the generic quality functor.
struct FakeGenericInternals {
    unsigned int num_nodes;
    unsigned int num_null_model_vectors;
    std::vector<double> comm_w_in;
    vec2 comm_loss_vectors;
};

void test_generic_quality_cancellation()
{
    // Gain and loss terms that cancel: the exact quality is 1 - t + 2t.
    FakeGenericInternals internals;
    internals.num_nodes = 4;
    internals.num_null_model_vectors = 2;
    internals.comm_w_in = {1.0, 1e17, 1.0, 0.0};
    internals.comm_loss_vectors = vec2{{0.0, 0.0, 0.0, 1e9}, {0.0, 0.0, 0.0, 1e8}};
    clq::find_linearised_generic_stability quality(1.0);
    CHECK(quality(internals) == 2.0);
}

}  // namespace

int main(int argc, char *argv[])
{
    std::string data_dir = ".";
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--print") == 0) {
            print_mode = true;
        } else {
            data_dir = argv[i];
        }
    }

    TestGraph pp = planted_partition(5, 20, 0.3, 0.03, 1u);
    TestGraph pp_norm = normalised_weights(pp);
    TestGraph pp_large = planted_partition(10, 40, 0.15, 0.01, 2u);

    for (unsigned int seed = 0; seed < 4; ++seed) {
        test_generic("generic_modularity_seed" + std::to_string(seed), pp_norm,
                     modularity_null_model(pp), 1.0, seed);
    }
    test_generic("generic_modularity_t0.3", pp_norm, modularity_null_model(pp), 0.3, 5u);
    test_generic("generic_two_pairs_t0.7", pp_norm, two_pair_null_model(pp, 3u), 0.7, 6u);
    test_generic("generic_large_t1", normalised_weights(pp_large),
                 modularity_null_model(pp_large), 1.0, 8u);
    test_standard<clq::find_linearised_normalised_stability,
                  clq::linearised_normalised_stability_gain>("normalised_t1", pp, 1.0, 9u, true);
    test_standard<clq::find_linearised_normalised_stability,
                  clq::linearised_normalised_stability_gain>("normalised_large_t0.5", pp_large,
                          0.5, 10u, true);
    test_standard<clq::find_linearised_combinatorial_stability,
                  clq::linearised_combinatorial_stability_gain>("combinatorial_t0.5", pp, 0.5,
                          11u, false);
    test_compensated_sum();
    test_generic_quality_cancellation();
    test_planted_partition_is_found();
    test_triangle_files(data_dir);

    if (print_mode) {
        return 0;
    }
    std::printf("%d checks, %d failures\n", num_checks, num_failures);
    return num_failures == 0 ? 0 : 1;
}
