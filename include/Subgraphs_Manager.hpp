
#ifndef SUBGRAPHS_MANAGER_HPP
#define SUBGRAPHS_MANAGER_HPP



#include "Super_Kmer_Bucket.hpp"
#include "HyperLogLog.hpp"
#include "Directed_Vertex.hpp"
#include "DNA_Utility.hpp"
#include "Discontinuity_Graph.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Character_Buffer.hpp"
#include "utility.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <limits>
#include <atomic>
#include <string>
#include <vector>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Manager for the subgraphs of the de Bruijn graph—manages their super k-mer
// based sequence representations in buckets, constructs them from these
// representations, and contracts them into their compacted form. `Colored_`
// denotes whether the vertices in the graph have associated colors.
template <uint16_t k, bool Colored_>
class Subgraphs_Manager
{
private:

    const std::string path_pref;    // Path-prefix to the super k-mer buckets.
    const std::size_t graph_count_; // Number of subgraphs; it needs to be a power of 2.
    const uint16_t l;   // `l`-minimizer size to partition the graph.

    typedef Super_Kmer_Bucket<Colored_> bucket_t;
    std::vector<Padded_Data<bucket_t>> subgraph_bucket; // Super k-mer buckets for the subgraphs.

    std::vector<HyperLogLog> HLL;   // `HLL[g]` is the cardinality-estimator for subgraph `g`.

    Discontinuity_Graph<k>& G;  // The discontinuity graph.

    std::atomic_uint64_t trivial_mtig_count_;   // Number of trivial maximal unitigs in the subgraphs (i.e. also maximal unitigs in the supergraph).
    std::atomic_uint64_t icc_count_;    // Number of trivial maximal unitigs in the subgraphs that are ICCs.

    // TODO: move out the following to some CF3-centralized location.

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded_Data<op_buf_t>> op_buf_list_t;
    op_buf_list_t& op_buf; // Worker-specific output buffers.


public:

    // Constructs a manager for the subgraphs of a de Bruijn graph which is
    // partitioned into `graph_count` subgraphs according to `l`-minimizers.
    // `logistics` is the data-logistics manager for the algorithm execution.
    // The discontinuity-graph is produced at `G` without false-phantom edges.
    // Worker-specific trivially maximal unitigs are written to the buffers in
    // `op_buf`.
    Subgraphs_Manager(const Data_Logistics& logistics, std::size_t graph_count, uint16_t l, Discontinuity_Graph<k>& G, op_buf_list_t& op_buf);

    // Returns the number of subgraphs.
    auto graph_count() const { return graph_count_; }

    // Adds a (weak) super k-mer to the subgraph `g` of the de Bruijn graph
    // with label `seq` and length `len`. The markers `l_disc` and `r_disc`
    // denote whether the left and the right ends of the (weak) super k-mer are
    // discontinuous or not.
    void add_super_kmer(std::size_t g, const char* seq, std::size_t len, bool l_disc, bool r_disc);

    // Finalizes the subgraphs for iteration—no more content should be added
    // after this.
    void finalize();

    // Returns the largest estimated size of any subgraph.
    uint64_t estimate_size_max() const;

    // Constructs and contracts each subgraph.
    void process();

    // Returns the number of trivial maximal unitigs in the subgraphs (i.e. also
    // maximal unitigs in the supergraph).
    uint64_t trivial_mtig_count() const;

    // Returns the number of trivial maximal unitigs in the subgraphs that are
    // ICCs.
    uint64_t icc_count() const;

    // Returns the subgraph ID for a minimizer with 64-bit hash value `h`.
    uint64_t graph_ID(uint64_t h) const { return h & (graph_count_ - 1); }
};


template <uint16_t k, bool Colored_>
inline void Subgraphs_Manager<k, Colored_>::add_super_kmer(const std::size_t g, const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc)
{
    assert(len >= k);

    auto& bucket = subgraph_bucket[g].data();
    bucket.add(seq, len, l_disc, r_disc);

    auto& hll = HLL[g];
    Directed_Vertex<k> v{Kmer<k>(seq)};
    constexpr auto u32_mask = std::numeric_limits<uint32_t>::max();
    std::size_t next_idx = k;
    while(true)
    {
        hll.add(v.canonical().to_u64() & u32_mask);

        if(next_idx == len)
            break;

        v.roll_forward(DNA_Utility::map_base(seq[next_idx]));
        next_idx++;
    }
}

}



#endif
