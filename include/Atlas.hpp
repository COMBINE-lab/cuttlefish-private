
#ifndef ATLAS_HPP
#define ATLAS_HPP



#include "Super_Kmer_Chunk.hpp"
#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>


namespace cuttlefish
{

template <bool Colored_>
class Atlas
{
private:

    typedef typename Super_Kmer_Chunk<Colored_>::attribute_t attribute_t;
    typedef typename Super_Kmer_Chunk<Colored_>::label_unit_t label_unit_t;


    static constexpr uint64_t atlas_count_ = 128;   // Number of subgraph atlases.
    static constexpr uint64_t graph_per_atlas_ = 128;    // Number of subgraphs per atlas.
    static constexpr uint64_t graph_count_ = atlas_count_ * graph_per_atlas_;   // Number of subgraphs.


    const std::string path_;    // Path to the external-memory bucket.

    uint64_t size_; // Number of super k-mers in the atlas. It's not necessarily correct before closing it.

    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    chunk_t chunk;  // Super k-mer chunk for the bucket.
    std::vector<Padded<chunk_t>> chunk_w;   // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

    Spin_Lock lock; // Lock to the chunk.


public:

    // Returns the number of subgraph atlases.
    static constexpr auto atlas_count() { return atlas_count_; }

    // Returns the number of subgraphs per atlas.
    static constexpr auto graph_per_atlas() { return graph_per_atlas_; }

    // Returns the number of subgraphs.
    static constexpr auto graph_count() { return graph_count_; }

    // Returns the atlas-ID of the `g`'th subgraph.
    static auto atlas_ID(const uint64_t g) { return g >> log_2(graph_per_atlas()); }

    // Returns the graph-ID of the `g`'th subgraph within its atlas.
    static auto graph_ID(const uint64_t g) { return g & (graph_per_atlas() - 1); }

    // Constructs a super k-mer atlas for `k`-mers and `l`-minimizers, at
    // external-memory path-prefix `path`. The super chunk buffer of the atlas
    // will have a soft capacity of `chunk_cap` and each worker-local buffer
    // will have a hard capacity of `chunk_cap_per_w`.
    Atlas(uint16_t k, uint16_t l, const std::string& path, std::size_t chunk_cap, std::size_t chunk_cap_per_w);

    Atlas(Atlas&&);

    Atlas(const Atlas&) = delete;
    Atlas& operator=(const Atlas&) = delete;
    Atlas& operator=(Atlas&&) = delete;

    // Returns the number of super k-mers in the atlas. It's not necessarily
    // correct before closing it.
    auto size() const { return size_; }

    // Returns the size of the atlas in bytes. It's not necessarily correct
    // before closing the bucket.
    auto bytes() const { return size() * chunk.record_size(); }

    // Adds a super k-mer to the atlas with label `seq` and length `len`. The
    // markers `l_disc` and `r_disc` denote whether the left and the right ends
    // of the (weak) super k-mer are discontinuous or not. The associated super
    // k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, bool l_disc, bool r_disc, uint16_t g_id);

    // Adds a super k-mer to the atlas with label `seq` and length `len` from
    // source-ID `source`. The markers `l_disc` and `r_disc` denote whether the
    // left and the right ends of the (weak) super k-mer are discontinuous or
    // not. The associated super k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, source_id_t source, bool l_disc, bool r_disc, uint16_t g_id);

    // Closes the atlas—no more content should be added afterwards.
    void close();

    // Removes the atlas.
    void remove();
};

}



#endif
