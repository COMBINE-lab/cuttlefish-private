
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

    typedef typename Super_Kmer_Chunk<Colored_>::attribute_t attribute_t;
    typedef typename Super_Kmer_Chunk<Colored_>::label_unit_t label_unit_t;

private:

    static constexpr uint64_t graph_per_atlas = 128;

    const uint16_t k;   // k-mer length.
    const uint16_t l;   // Minimizer size.

    const std::string path_;    // Path to the external-memory bucket.
    std::ofstream output;   // Output stream to the external-memory bucket.

    uint64_t size_; // Number of super k-mers in the bucket. It's not necessarily correct before closing the bucket.

    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    // TODO: make informed choices for the chunk-sizes based on whether atlases are used or not.
    static constexpr std::size_t chunk_bytes = 128 * 1024;  // 128 KB chunk capacity.
    static constexpr std::size_t w_chunk_bytes = 32 * 1024; // 32 KB worker-chunk capacity. // TODO: needs to smaller in-case graph-atlases aren't used.
    const std::size_t chunk_cap;    // Capacity (in number of super k-mers) of the chunk of the bucket.
    mutable chunk_t chunk;  // Super k-mer chunk for the bucket.    // TODO: maybe this is not required. `chunk_w[i]` can bypass this to disk.
    std::vector<Padded<chunk_t>> chunk_w;   // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

    std::vector<uint32_t> src_hist; // Frequency histogram of super k-mer sources currently committed to the chunk.

    std::vector<uint32_t> chunk_sz; // Sizes of the flushed chunks.

    Spin_Lock lock; // Lock to the chunk and the external-memory bucket.


    // Empties the local chunk of worker `w_id` to the chunk of the bucket in a
    // thread-safe manner.
    void empty_w_local_chunk(std::size_t w_id);

    // Flushes the super k-mer chunk to the external-memory bucket.
    void flush_chunk();

    // Shatters the super k-mer chunk `c` to the buckets in `B`.
    void shatter_chunk(const Super_Kmer_Chunk<Colored_>& c, std::vector<Padded<Super_Kmer_Bucket>>& B);

public:

    // Constructs a super k-mer bucket for `k`-mers and `l`-minimizers, at
    // external-memory path `path`.
    Super_Kmer_Bucket(uint16_t k, uint16_t l, const std::string& path);

    Super_Kmer_Bucket(const Super_Kmer_Bucket&) = delete;

    Super_Kmer_Bucket(Super_Kmer_Bucket&& rhs);

    // Allocates the worker-local chunks' memories.
    void allocate_worker_mem();

    // Deallocates the worker-local chunks' memories.
    void deallocate_worker_mem();

    // Returns the number of super k-mers in the bucket. It's not necessarily
    // correct before closing the bucket.
    auto size() const { return size_; }

    // Adds a super k-mer to the bucket with label `seq` and length `len`. The
    // markers `l_disc` and `r_disc` denote whether the left and the right ends
    // of the (weak) super k-mer are discontinuous or not. The associated super
    // k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, bool l_disc, bool r_disc, uint16_t g_id);

    // Adds a super k-mer to the bucket with label `seq` and length `len` from
    // source-ID `source`. The markers `l_disc` and `r_disc` denote whether the
    // left and the right ends of the (weak) super k-mer are discontinuous or
    // not. The associated super k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, source_id_t source, bool l_disc, bool r_disc, uint16_t g_id);

    // Adds a super k-mer to the chunk with encoding `seq` and attributes
    // `att`.
    void add(const label_unit_t* seq, const attribute_t& att);

    // Collates the worker-local buffers into the external-memory bucket and
    // empties them.
    void collate_buffers();

    // Closes the bucket—no more content should be added afterwards.
    void close();

    // Removes the bucket.
    void remove();

    // Returns an iterator over the super k-mers in the bucket. The bucket
    // should be closed before iteration.
    Iterator iterator() const { assert(chunk.empty()); return Iterator(*this); }

    // Shatters the bucket into the buckets in `B`.
    void shatter(std::vector<Padded<Super_Kmer_Bucket>>& B);
};


// Iterator over super k-mer buckets.
template <bool Colored_>
class Super_Kmer_Bucket<Colored_>::Iterator
{
    friend class Super_Kmer_Bucket<Colored_>;

private:

    const Super_Kmer_Bucket& B; // Bucket to iterate over.
    std::ifstream input;    // Input stream from the external-memory bucket.

    std::size_t idx;    // Current slot-index the iterator is in, i.e. next super k-mer to access.
    std::size_t chunk_start_idx;    // Index into the bucket where the current in-memory chunk starts.
    std::size_t chunk_end_idx;  // Non-inclusive index into the bucket where the current in-memory chunk ends.
    std::size_t chunk_id;   // Sequential-ID of the chunk being processed right now.


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
    bool next(attribute_t& att, const label_unit_t*& label);
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
inline void Super_Kmer_Bucket<false>::add(const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc, const uint16_t g_id)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap(); // Worker-specific chunk.

    assert(c_w.size() < c_w.capacity());
    c_w.add(seq, len, l_disc, r_disc, g_id);

    if(c_w.full())
        empty_w_local_chunk(w_id);
}


template <>
inline void Super_Kmer_Bucket<true>::add(const char* const seq, const std::size_t len, const source_id_t source, const bool l_disc, const bool r_disc, const uint16_t g_id)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap(); // Worker-specific chunk.

    c_w.add(seq, len, source, l_disc, r_disc, g_id);
    // No flush until collation is invoked explicitly from outside.
}


template <>
inline void Super_Kmer_Bucket<false>::add(const label_unit_t* const seq, const attribute_t& att)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap(); // Worker-specific chunk.

    assert(c_w.size() < c_w.capacity());
    c_w.add(seq, att);

    if(c_w.full())
        empty_w_local_chunk(w_id);
}


template <>
inline void Super_Kmer_Bucket<true>::add(const label_unit_t* const seq, const attribute_t& att)
{
    // TODO: implement.
    (void)seq, (void)att;
}


template <bool Colored_>
inline bool Super_Kmer_Bucket<Colored_>::Iterator::next(attribute_t& att, const label_unit_t*& label)
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
