
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
    // static constexpr std::size_t chunk_bytes = 4 * 1024;    // 4 KB chunk capacity.
    static constexpr std::size_t chunk_bytes = 128 * 1024;  // 128 KB chunk capacity.
    static constexpr std::size_t w_chunk_bytes = 1 * 1024;  // 1 KB worker-chunk capacity.
    const std::size_t chunk_cap;    // Capacity (in number of super k-mers) of the chunk of the bucket.
    mutable chunk_t chunk;  // Super k-mer chunk for the bucket.    // TODO: maybe this is not required. `chunk_w[i]` can bypass this to disk.
    mutable std::vector<Padded<chunk_t>> chunk_w;   // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

    std::vector<Padded<uint32_t>> latest_src_w; // Latest source-ID of committed super k-mers to the chunk, per worker.
    uint32_t safe_src_lower_bound;  // Tight lower bound of source-IDs of super k-mers in the chunk safe to flush.
    uint32_t safe_src_upper_bound;  // Tight upper bound of source-IDs of super k-mers in the chunk safe to flush.
    uint64_t safe_c;    // Count of super k-mers in the chunk safe to flush.
    std::vector<uint32_t> src_hist; // Frequency histogram of super k-mer sources currently committed to the chunk.

    mutable chunk_t chunk_safe; // Chunk of super k-mers safe to flush to the external-memory bucket in the colored case.
    std::vector<uint32_t> chunk_sz; // Sizes of the flushed chunks; only applicable in the colored case.

    mutable Spin_Lock lock; // Lock to the chunk and the external-memory bucket.


    // Empties the local chunk of worker `w_id` to the chunk of the bucket in a
    // thread-safe manner.
    void empty_w_local_chunk(std::size_t w_id);

    // Counting-sorts the committed super k-mers (i.e. in `chunk`) that are
    // safe to flush and appends them in order to `chunk_safe`. The unsafe
    // super k-mers are moved to the front of `chunk`.
    void counting_sort_safe_super_kmers();

    // Returns the ID of `src` relative to the bucket's lower bound of safe
    // source-ID.
    uint32_t rel_src(const uint32_t src) const { return src - safe_src_lower_bound; }

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
inline void Super_Kmer_Bucket<true>::counting_sort_safe_super_kmers()
{
    assert(safe_c > 0);

    chunk_safe.resize(safe_c);
    auto& src_off = src_hist;   // Offsets for super k-mers in the safe chunk, for counting sort.
    uint64_t pref_sum = 0;
    const auto safe_src_c = safe_src_upper_bound - safe_src_lower_bound + 1;
    for(std::size_t i = 0; i < safe_src_c; ++i)
    {
        const auto temp = src_off[i];
        src_off[i] = pref_sum;
        pref_sum += temp;
    }
    assert(pref_sum == safe_c);


    uint64_t unsafe_c = 0;  // Count of super k-mers in the chunk unsafe to flush.
    for(std::size_t i = 0, j; i < chunk.size(); i = j)
    {
        const auto src = chunk.att_at(i).source();
        for(j = i + 1; j < chunk.size(); ++j)
            if(chunk.att_at(j).source() != src)
                break;

        const auto stretch_sz = j - i;
        if(src <= safe_src_upper_bound) // Add this stretch of super k-mers to the safe chunk, counting-sorted.
        {
            const auto src_rel = rel_src(src);
            chunk_safe.copy(src_off[src_rel], chunk, i, stretch_sz);

            src_off[src_rel] += stretch_sz;
            assert(src_rel < safe_src_c || src_off[src_rel] <= src_off[src_rel + 1]);
        }
        else    // Move this stretch of super k-mers toward the front.
        {
            chunk.move(unsafe_c, i, stretch_sz);

            unsafe_c += stretch_sz;
        }
    }

    chunk.resize(unsafe_c);
}


template <>
inline void Super_Kmer_Bucket<true>::empty_w_local_chunk(const std::size_t w_id)
{
    auto& c_w = chunk_w[w_id].unwrap();
    if(c_w.empty())
        return;

    lock.lock();

    // Commit this worker's chunk to the bucket's chunk.
    chunk.reserve(chunk.size() + c_w.size());
    chunk.append(c_w);
    size_ += c_w.size();


    // Update bucket-chunk's source histogram for this new addition.
    const auto max_src = c_w.back_att().source();   // Maximum source-ID in the worker chunk.
    latest_src_w[w_id].unwrap() = max_src;
    if(src_hist.size() < rel_src(max_src) + 1)
        src_hist.resize(rel_src(max_src) + 1);

    for(std::size_t i = 0; i < c_w.size(); ++i)
        src_hist[rel_src(c_w.att_at(i).source())]++;


    // Update the upper bound of safe source-IDs for the current version of bucket-chunk.
    safe_src_upper_bound = std::numeric_limits<uint32_t>::max();
    std::for_each(latest_src_w.cbegin(), latest_src_w.cend(), [&](const auto& v){ safe_src_upper_bound = std::min(safe_src_upper_bound, v.unwrap()); });

    safe_c = 0;
    for(std::size_t i = safe_src_lower_bound; i <= safe_src_upper_bound; ++i)
        safe_c += src_hist[rel_src(i)];


    if(safe_c >= chunk_cap / 2) // Enough safe super k-mers exists to warrant expensive works: sort and data relocations.
    {
        counting_sort_safe_super_kmers();

        // Move frequencies of the unsafe sources to the front.
        const auto safe_src_c = safe_src_upper_bound - safe_src_lower_bound + 1;
        const auto unsafe_src_c = src_hist.size() - safe_src_c;
        if(unsafe_src_c > 0)
            std::memmove(src_hist.data() + 1, src_hist.data() + safe_src_c, unsafe_src_c * sizeof(decltype(src_hist)::value_type));
        src_hist.resize(1 + unsafe_src_c);  // Source-ID `safe_src_lower_bound` cannot be discarded yet.

        safe_src_lower_bound = safe_src_upper_bound;
        src_hist[rel_src(safe_src_lower_bound)] = 0;

        chunk_safe.serialize(output);
        chunk_sz.push_back(chunk_safe.size());
        chunk_safe.clear();
    }
    // else if(...)
    // TODO: in case chunk accumulates too much, empty all worker buffers into it, obtain globally latest source
    // per worker and set safe upper bound accordingly, and then filter and flush; i.e. keep unbounded increase in check.

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
