
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
    static constexpr std::size_t chunk_bytes = 4 * 1024;    // 4 KB chunk capacity.
    const std::size_t chunk_cap;    // Capacity (in number of super k-mers) of the chunk of the bucket.
    mutable chunk_t chunk;  // Super k-mer chunk for the bucket.    // TODO: maybe this is not required. `chunk_w[i]` can bypass this to disk.
    mutable std::vector<Padded<chunk_t>> chunk_w;   // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

    std::vector<Padded<uint32_t>> latest_src_w; // Latest source-ID of committed super k-mers to the chunk, per worker.
    uint32_t latest_safe_src = 1;   // Latest source-ID of super k-mers safe to flush.
    std::vector<uint32_t> src_hist; // Frequency histogram of super k-mer sources currently committed to the chunk.

    mutable chunk_t chunk_safe; // Chunk of super k-mers safe to flush to the external-memory bucket in the colored case.
    std::vector<uint32_t> chunk_sz; // Sizes of the flushed chunks; only applicable in the colored case.

    mutable Spin_Lock lock; // Lock to the chunk and the external-memory bucket.


    // Empties the local chunk of worker `w_id` to the chunk of the bucket in a
    // thread-safe manner.
    void empty_w_local_chunk(std::size_t w_id);

    // Counting-sorts the committed super k-mers (i.e. in `chunk`) that are
    // safe to flush and appends them in order to `chunk_safe`. The unsafe
    // super k-mers are moved to the front of `chunk`. Returns the latest
    // source-ID of the sorted super k-mers.
    uint32_t counting_sort_safe_super_kmers();

    // Returns the chunk to use in reading super k-mers during iteration.
    constexpr chunk_t* reader_chunk() const { if constexpr(!Colored_) return &chunk; return &chunk_safe; }

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
        chunk.serialize(output);
        chunk.clear();

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
    auto& c_w = chunk_w[w_id].unwrap();   // Worker-specific chunk.

    assert(c_w.size() < c_w.capacity());
    c_w.add(seq, len, l_disc, r_disc);

    if(c_w.full())
        empty_w_local_chunk(w_id);
}


template <>
inline uint32_t Super_Kmer_Bucket<true>::counting_sort_safe_super_kmers()
{
    uint32_t src_th = std::numeric_limits<uint32_t>::max(); // Upper-bound of source-IDs of super k-mers in the chunk safe to flush.
    std::for_each(latest_src_w.cbegin(), latest_src_w.cend(), [&](const auto& v){ src_th = std::min(src_th, v.unwrap()); });

    src_hist.resize(src_th - latest_safe_src + 1);
    std::memset(src_hist.data(), 0, src_hist.size() * sizeof(decltype(src_hist)::value_type));
    uint64_t safe_c = 0;    // Count of super k-mers safe to flush.
    for(std::size_t i = 0; i < chunk.size(); ++i)
    {
        const auto src = chunk.att_at(i).source();
        if(src <= src_th)
            src_hist[src - latest_safe_src]++,
            safe_c++;
    }

    // assert(safe_c > 0);  // TODO: diagnose performance characteristics of this assert failing.


    // Base offset for new super k-mers into the safe chunk.
    const auto base_off = chunk_safe.size();
    chunk_safe.resize(chunk_safe.size() + safe_c);
    auto& src_off = src_hist;   // Offsets for super k-mers in the safe chunk, for counting sort.
    uint64_t pref_sum = 0;
    for(std::size_t i = 0; i < src_off.size(); ++i)
    {
        const auto temp = src_off[i];
        src_off[i] = pref_sum;
        pref_sum += temp;
    }
    assert(pref_sum == safe_c);


    std::size_t unsafe_c = 0;   // Count of super k-mers unsafe to flush.
    for(std::size_t i = 0, j; i < chunk.size(); i = j)
    {
        const auto src = chunk.att_at(i).source();
        for(j = i + 1; j < chunk.size(); ++j)
            if(chunk.att_at(j).source() != src)
                break;

        const auto stretch_sz = j - i;
        if(src <= src_th)   // Add this stretch of super k-mers to the safe chunk, counting-sorted.
        {
            const auto src_rel = src - latest_safe_src;
            chunk_safe.copy(base_off + src_off[src_rel], chunk, i, stretch_sz);
            src_off[src_rel] += stretch_sz;
            assert(src_rel + 1 == src_off.size() || src_off[src_rel] <= src_off[src_rel + 1]);
        }
        else    // Move this stretch of super k-mers toward the front.
        {
            chunk.move(unsafe_c, i, stretch_sz);
            unsafe_c += stretch_sz;
        }
    }

    chunk.resize(unsafe_c);

    return src_th;
}


template <>
inline void Super_Kmer_Bucket<true>::empty_w_local_chunk(const std::size_t w_id)
{
    auto& c_w = chunk_w[w_id].unwrap();
    if(c_w.empty())
        return;

    lock.lock();

    chunk.reserve(chunk.size() + c_w.size());
    chunk.append(c_w);
    size_ += c_w.size();

    latest_src_w[w_id].unwrap() = c_w.back_att().source();


    if(chunk.size() >= chunk_cap)
    {
        latest_safe_src = counting_sort_safe_super_kmers();

        if(chunk_safe.size() >= chunk_cap)
        {
            chunk_safe.serialize(output);
            chunk_sz.push_back(chunk_safe.size());
            chunk_safe.clear();
        }
    }

    lock.unlock();

    c_w.clear();
}


template <>
inline void Super_Kmer_Bucket<true>::add(const char* const seq, const std::size_t len, const uint32_t source, const bool l_disc, const bool r_disc)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap();   // Worker-specific chunk.

    assert(c_w.size() < c_w.capacity());
    c_w.add(seq, len, source, l_disc, r_disc);

    if(c_w.full())
        empty_w_local_chunk(w_id);
}


template <bool Colored_>
inline bool Super_Kmer_Bucket<Colored_>::Iterator::next(attribute_t& att, label_unit_t*& label)
{
    assert(idx <= B.size());

    auto* const c = B.reader_chunk();

    if(CF_UNLIKELY(idx == B.size()))
    {
        c->clear();
        return false;
    }

    if(idx == chunk_end_idx)
    {
        chunk_start_idx = chunk_end_idx;
        chunk_end_idx += read_chunk();
    }

    assert(idx >= chunk_start_idx && idx < chunk_end_idx);
    c->get_super_kmer(idx - chunk_start_idx, att, label);
    idx++;

    return true;
}

}



#endif
