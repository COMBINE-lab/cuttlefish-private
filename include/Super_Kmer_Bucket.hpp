
#ifndef SUPER_KMER_BUCKET_HPP
#define SUPER_KMER_BUCKET_HPP



#include "Super_Kmer_Chunk.hpp"
#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <cstddef>
#include <string>
#include <vector>
#include <fstream>


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
    std::ofstream output;   // Output stream to the external-memory bucket.

    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    static constexpr std::size_t chunk_bytes = 4 * 1024;    // 4 KB chunk capacity.
    const std::size_t chunk_cap;    // Capacity (in number of super k-mers) of the chunk of the bucket.
    chunk_t chunk;  // Super k-mer chunk for the bucket.
    std::vector<Padded_Data<chunk_t>> chunk_w;  // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

    mutable Spin_Lock lock; // Lock to the chunk and the external-memory bucket.


    // Empties the local chunk of worker `w_id` to the chunk of the bucket in a
    // thread-safe manner.
    void empty_w_local_chunk(std::size_t w_id);

public:

    // Constructs a super k-mer bucket for `k`-mers and `l`-minimizers, at
    // external-memory path `path`.
    Super_Kmer_Bucket(uint16_t k, uint16_t l, const std::string& path);

    Super_Kmer_Bucket(const Super_Kmer_Bucket&) = delete;

    Super_Kmer_Bucket(Super_Kmer_Bucket&& rhs);

    // Adds a super k-mer to the bucket with label `seq` and length `len`. The
    // markers `l_disc` and `r_disc` denote whether the left and the right ends
    // of the (weak) super k-mer are discontinuous or not.
    void add(const char* seq, std::size_t len, bool l_disc, bool r_disc);

    // Closes the bucket—no more content should be added afterwards.
    void close();
};


template <bool Colored_>
inline void Super_Kmer_Bucket<Colored_>::add(const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].data();   // Worker-specific chunk.

    assert(c_w.size() < c_w.capacity());
    c_w.add(seq, len, l_disc, r_disc);

    if(c_w.size() == c_w.capacity())
        empty_w_local_chunk(w_id);
}


template <bool Colored_>
inline void Super_Kmer_Bucket<Colored_>::empty_w_local_chunk(const std::size_t w_id)
{
    auto& c_w = chunk_w[w_id].data();

    lock.lock();

    const auto break_idx = std::min(c_w.size(), chunk.free_capacity());
    chunk.append(c_w, 0, break_idx);
    if(chunk.full())
    {
        chunk.serialize(output);
        chunk.clear();

        if(break_idx < c_w.size())
            assert(chunk.capacity() >= c_w.size() - break_idx),
            chunk.append(c_w, break_idx, c_w.size());
    }

    lock.unlock();

    c_w.clear();
}

}



#endif
