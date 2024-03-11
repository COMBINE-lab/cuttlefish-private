
#ifndef SUPER_KMER_BUCKET_HPP
#define SUPER_KMER_BUCKET_HPP



#include "Super_Kmer_Chunk.hpp"
#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <cstddef>
#include <string>
#include <vector>


namespace cuttlefish
{

// =============================================================================
// A bucket of super k-mers corresponding to a subgraph of the underlying de
// Bruijn graph. `Colored_` denotes whether the super k-mers in the bucket each
// has an associated source ID.
template <bool Colored_>
class Super_Kmer_Bucket
{

private:

    const std::string path_;    // Path to the external-memory bucket.

    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    static constexpr std::size_t chunk_bytes = 4 * 1024;    // 4 KB chunk capacity.
    const std::size_t chunk_cap;    // Capacity (in number of super k-mers) of the chunk of the bucket.
    chunk_t chunk;  // Super k-mer chunk for the bucket.
    Padded_Data<chunk_t>* const chunk_w;    // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

public:

    // Constructs a super k-mer bucket for `k`-mers and `l`-minimizers, at
    // external-memory path `path`.
    Super_Kmer_Bucket(uint16_t k, uint16_t l, const std::string& path);

    ~Super_Kmer_Bucket();
};

}



#endif
