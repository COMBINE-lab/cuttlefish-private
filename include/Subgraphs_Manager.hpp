
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
#include <deque>
#include <limits>
#include <atomic>
#include <string>
#include <vector>
#include <type_traits>


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

    // TODO: centralize these params.
    static constexpr uint64_t atlas_count = 128;    // Number of subgraph atlases.
    static constexpr uint64_t graph_per_atlas = 128;    // Number of subgraphs per atlas.
    static constexpr uint64_t graph_count_ = atlas_count * graph_per_atlas; // Number of subgraphs.

    const std::string path_pref;    // Path-prefix to the super k-mer buckets.
    const uint16_t l;   // `l`-minimizer size to partition the graph.

    typedef Super_Kmer_Bucket<Colored_> bucket_t;
    std::vector<Padded<bucket_t>> atlas;    // Super k-mer buckets for the subgraph atlases.
    std::vector<Padded<bucket_t>> subgraph_bucket;  // Super k-mer buckets for the subgraphs.
    
    static constexpr std::size_t chunk_bytes = 128 * 1024;  // 128 KB chunk capacity.
    static constexpr std::size_t w_chunk_bytes = 32 * 1024; // 32 KB worker-chunk capacity. // TODO: needs to be smaller in-case graph-atlases aren't used.
    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    std::vector<Padded<chunk_t>> chunk; // Chunks to port into different buckets in an atlas.
    std::vector<Padded<std::vector<Padded<chunk_t>>>> chunk_w;  // Worker-local chunk collections to port into different buckets in an atlas.

    std::vector<Padded<HyperLogLog>> HLL;   // `HLL[g]` is the cardinality-estimator for subgraph `g`.

    Discontinuity_Graph<k, Colored_>& G_;   // The discontinuity graph.

    std::atomic_uint64_t trivial_mtig_count_;   // Number of trivial maximal unitigs in the subgraphs (i.e. also maximal unitigs in the supergraph).
    std::atomic_uint64_t icc_count_;    // Number of trivial maximal unitigs in the subgraphs that are ICCs.

    // TODO: move out the following to some CF3-centralized location.

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded<op_buf_t>> op_buf_list_t;
    op_buf_list_t& op_buf; // Worker-specific output buffers.

    const std::string color_path_pref;  // Path-prefix to the output color buckets.

    // Allocates space for the chunks and worker-local chunks of different
    // subgraphs in an atlas.
    void allocate_chunks();

    // Returns the atlas-ID of the `g`'th subgraph.
    uint64_t atlas_ID(const uint64_t g) const { return g >> log_2(atlas_count); }

    // Adds the label `seq` and length `len` to the HLL estimate of the
    // subgraph `g` of the de Bruijn graph.
    void add_to_HLL(std::size_t g, const char* seq, std::size_t len);


public:

    // Constructs a manager for the subgraphs of a de Bruijn graph which is
    // partitioned according to `l`-minimizers. `logistics` is the data-
    // logistics manager for the algorithm execution. The discontinuity-graph
    // is produced at `G` without false-phantom edges. Worker-specific
    // trivially maximal unitigs are written to the buffers in `op_buf`.
    Subgraphs_Manager(const Data_Logistics& logistics, uint16_t l, Discontinuity_Graph<k, Colored_>& G, op_buf_list_t& op_buf);

    // Returns the number of subgraphs.
    auto graph_count() const { return graph_count_; }

    // Returns the discontinuity graph.
    const auto& G() const { return G_; }

    // Adds a (weak) super k-mer to the subgraph `g` of the de Bruijn graph
    // with label `seq` and length `len`. The markers `l_disc` and `r_disc`
    // denote whether the left and the right ends of the (weak) super k-mer are
    // discontinuous or not.
    template <bool C_ = Colored_, std::enable_if_t<!C_, int> = 0>
    void add_super_kmer(std::size_t g, const char* seq, std::size_t len, bool l_disc, bool r_disc);

    // Adds a (weak) super k-mer to the subgraph `g` of the de Bruijn graph
    // with label `seq` and length `len` from source-ID `source`. The markers
    // `l_disc` and `r_disc` denote whether the left and the right ends of the
    // (weak) super k-mer are discontinuous or not.
    template <bool C_ = Colored_, std::enable_if_t<C_, int> = 0>
    void add_super_kmer(std::size_t g, const char* seq, std::size_t len, source_id_t source, bool l_disc, bool r_disc);

    // Collates the current super k-mer buffers in each subgraph per their
    // source-IDs into external-memory buckets.
    void collate_super_kmer_buffers();

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
template <bool C_, std::enable_if_t<!C_, int>>
inline void Subgraphs_Manager<k, Colored_>::add_super_kmer(const std::size_t g, const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc)
{
    assert(len >= k);

    auto& bucket = atlas[atlas_ID(g)].unwrap();
    bucket.add(seq, len, l_disc, r_disc, g);

    // add_to_HLL(g, seq, len);
}


template <uint16_t k, bool Colored_>
template <bool C_, std::enable_if_t<C_, int>>
inline void Subgraphs_Manager<k, Colored_>::add_super_kmer(const std::size_t g, const char* const seq, const std::size_t len, const source_id_t source, const bool l_disc, const bool r_disc)
{
    assert(len >= k);

    auto& bucket = atlas[atlas_ID(g)].unwrap();
    bucket.add(seq, len, source, l_disc, r_disc, g);

    // add_to_HLL(g, seq, len);
}


template <uint16_t k, bool Colored_>
inline void Subgraphs_Manager<k, Colored_>::add_to_HLL(const std::size_t g, const char* const seq, const std::size_t len)
{
    auto& hll = HLL[g].unwrap();
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
