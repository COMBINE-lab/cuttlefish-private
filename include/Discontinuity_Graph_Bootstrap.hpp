
#ifndef DISCONTINUITY_GRAPH_BOOTSTRAP_HPP
#define DISCONTINUITY_GRAPH_BOOTSTRAP_HPP



#include "Kmer.hpp"

#include <cstdint>
#include <cstddef>
#include <string>


namespace cuttlefish
{

// =============================================================================
// Bootstrapper of the corresponding discontinuity-graph of a de Bruijn graph
// (of `k`-mers) from its compacted variant.
template <uint16_t k>
class Discontinuity_Graph_Bootstrap
{
private:

    const std::string cdbg_path;    // Path to the compacted de Bruijn graph FASTA file.
    const uint16_t l;   // Minimizer size.
    const uint64_t minimizer_seed;  // Seed used in computing minimizer-hashes.

    // k-mer (super-)label of the ϕ-vertex in the discontinuity graph.  // TODO: do better.
    static constexpr const char phi_label[] =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const Kmer<k> phi;  // ϕ k-mer connected to each chain-end in the discontinuity graph.


public:

    // Constructs a graph-bootstrapper from the compacted de Bruijn graph at the
    // path `cdbg_path`. `l`-minimizers are used for the graph, and the seed-
    // value `seed` is used in minimizer-hashing.
    Discontinuity_Graph_Bootstrap(const std::string& cdbg_path, uint16_t l, uint64_t seed = 0);

    // Generates the discontinuity-graph with `part_count` vertex-partitions at
    // the path-prefix `path`.
    void generate(std::size_t part_count, const std::string& path) const;
};

}



#endif
