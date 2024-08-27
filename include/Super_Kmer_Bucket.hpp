
#ifndef SUPER_KMER_BUCKET_HPP
#define SUPER_KMER_BUCKET_HPP



#include "Super_Kmer_Chunk.hpp"
#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// A bucket of super k-mers corresponding to a subgraph of the underlying de
// Bruijn graph. `Colored_` denotes whether the super k-mers in the bucket each
// has an associated source ID.
template <bool Colored_>
class Super_Kmer_Bucket
{
    class Iterator;
    friend class Iterator;

private:

    const std::string path_;    // Path to the external-memory bucket.
    std::ofstream output;   // Output stream to the external-memory bucket.

    uint64_t size_; // Number of super k-mers in the bucket. It's not necessarily correct before closing the bucket.

    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    static constexpr std::size_t chunk_bytes = 128 * 1024;  // 128 KB chunk capacity.
    static constexpr std::size_t w_chunk_bytes = 1 * 1024;  // 1 KB worker-chunk capacity.
    const std::size_t chunk_cap;    // Capacity (in number of super k-mers) of the chunk of the bucket.
    mutable chunk_t chunk;  // Super k-mer chunk for the bucket.    // TODO: maybe this is not required. `chunk_w[i]` can bypass this to disk.
    std::vector<Padded<chunk_t>> chunk_w;   // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

    std::vector<uint32_t> src_hist; // Frequency histogram of super k-mer sources currently committed to the chunk.

    std::vector<uint32_t> chunk_sz; // Sizes of the flushed chunks; only applicable in the colored case.

    Spin_Lock lock; // Lock to the chunk and the external-memory bucket.


    // Empties the local chunk of worker `w_id` to the chunk of the bucket in a
    // thread-safe manner.
    void empty_w_local_chunk(std::size_t w_id);

    // Flushes the super k-mer chunk to the external-memory bucket.
    void flush_chunk();

public:

    // Constructs a super k-mer bucket for `k`-mers and `l`-minimizers, at
    // external-memory path `path`.
    Super_Kmer_Bucket(uint16_t k, uint16_t l, const std::string& path);

    Super_Kmer_Bucket(const Super_Kmer_Bucket&) = delete;

    Super_Kmer_Bucket(Super_Kmer_Bucket&& rhs);

    // Returns the number of super k-mers in the bucket. It's not necessarily
    // correct before closing the bucket.
    auto size() const { return size_; }

    // Adds a super k-mer to the bucket with label `seq` and length `len`. The
    // markers `l_disc` and `r_disc` denote whether the left and the right ends
    // of the (weak) super k-mer are discontinuous or not.
    void add(const char* seq, std::size_t len, bool l_disc, bool r_disc);

    // Adds a super k-mer to the bucket with label `seq` and length `len` from
    // source-ID `source`. The markers `l_disc` and `r_disc` denote whether the
    // left and the right ends of the (weak) super k-mer are discontinuous or
    // not.
    void add(const char* seq, std::size_t len, uint32_t source, bool l_disc, bool r_disc);

    // Collates the worker-local buffers into the external-memory bucket and
    // empties them.
    void collate_buffers();

    // Closes the bucket—no more content should be added afterwards.
    void close();

    // Removes the bucket.
    void remove();

    // Returns an iterator over the super k-mers in the bucket. The bucket
    // should be closed before iteration.
    Iterator iterator() const;
};


// Iterator over super k-mer buckets.
template <bool Colored_>
class Super_Kmer_Bucket<Colored_>::Iterator
{
    friend class Super_Kmer_Bucket<Colored_>;

private:

    const Super_Kmer_Bucket<Colored_>& B;   // Bucket to iterate over.
    std::ifstream input;    // Input stream from the external-memory bucket.

    std::size_t idx;    // Current slot-index the iterator is in, i.e. next super k-mer to access.
    std::size_t chunk_start_idx;    // Index into the bucket where the current in-memory chunk starts.
    std::size_t chunk_end_idx;  // Non-inclusive index into the bucket where the current in-memory chunk ends.
    std::size_t chunk_id;   // Sequential-ID of the chunk being processed right now; only applicable in the colored case.


    // Constructs an iterator for the super k-mer bucket `B`.
    Iterator(const Super_Kmer_Bucket& B);

    // Reads in the next super k-mer chunk from the bucket and returns the
    // number of super k-mers read.
    std::size_t read_chunk();


public:

    typedef typename Super_Kmer_Chunk<Colored_>::attribute_t attribute_t;
    typedef typename Super_Kmer_Chunk<Colored_>::label_unit_t label_unit_t;

    // Return the number of 64-bit words in super k-mer encodings.
    auto super_kmer_word_count() const { return B.chunk.super_kmer_word_count(); }

    // Moves the iterator to the next super k-mer in the bucket. Iff the bucket
    // is not depleted, the associated super k-mer's attribute and label-
    // encoding are put in `att` and `label` respectively, and returns `true`.
    bool next(attribute_t& att, label_unit_t*& label);
};


template <>
inline void Super_Kmer_Bucket<false>::empty_w_local_chunk(const std::size_t w_id)
{
    auto& c_w = chunk_w[w_id].unwrap();
    if(c_w.empty())
        return;

    lock.lock();

    // TODO: do away with `chunk`, and write `c_w` directly.

    const auto break_idx = std::min(c_w.size(), chunk.free_capacity());
    chunk.append(c_w, 0, break_idx);
    if(chunk.full())
    {
        flush_chunk();

        if(break_idx < c_w.size())
            assert(chunk.capacity() >= c_w.size() - break_idx),
            chunk.append(c_w, break_idx, c_w.size());
    }

    size_ += c_w.size();

    lock.unlock();

    c_w.clear();
}


template <>
inline void Super_Kmer_Bucket<false>::add(const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap(); // Worker-specific chunk.

    assert(c_w.size() < c_w.capacity());
    c_w.add(seq, len, l_disc, r_disc);

    if(c_w.full())
        empty_w_local_chunk(w_id);
}


template <>
inline void Super_Kmer_Bucket<true>::add(const char* const seq, const std::size_t len, const uint32_t source, const bool l_disc, const bool r_disc)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap(); // Worker-specific chunk.

    c_w.add(seq, len, source, l_disc, r_disc);
    // No flush until collation is invoked explicitly from outside.
}


template <bool Colored_>
inline bool Super_Kmer_Bucket<Colored_>::Iterator::next(attribute_t& att, label_unit_t*& label)
{
    assert(idx <= B.size());

    if(CF_UNLIKELY(idx == B.size()))
    {
        B.chunk.clear();
        return false;
    }

    if(idx == chunk_end_idx)
    {
        chunk_start_idx = chunk_end_idx;
        chunk_end_idx += read_chunk();
    }

    assert(idx >= chunk_start_idx && idx < chunk_end_idx);
    B.chunk.get_super_kmer(idx - chunk_start_idx, att, label);
    idx++;

    return true;
}

}



#endif
