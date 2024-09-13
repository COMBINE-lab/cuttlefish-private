
#include "Super_Kmer_Bucket.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Super_Kmer_Bucket(const uint16_t k, const uint16_t l, const std::string& path):
      path_(path)
    , output(path_, std::ios::out | std::ios::binary)
    , size_(0)
    , chunk_cap(chunk_bytes / Super_Kmer_Chunk<Colored_>::record_size(k, l))
    , chunk(k, l, chunk_cap)
{
    assert(chunk_cap >= parlay::num_workers());

    const auto chunk_cap_per_w = w_chunk_bytes / Super_Kmer_Chunk<Colored_>::record_size(k, l);
    chunk_w.reserve(parlay::num_workers());
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
        chunk_w.emplace_back(chunk_t(k, l, chunk_cap_per_w));
}


template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Super_Kmer_Bucket(Super_Kmer_Bucket&& rhs):
      path_(std::move(rhs.path_))
    , output(std::move(rhs.output))
    , size_(std::move(rhs.size_))
    , chunk_cap(std::move(rhs.chunk_cap))
    , chunk(std::move(rhs.chunk))
    , chunk_w(std::move(rhs.chunk_w))
    , src_hist(std::move(rhs.src_hist))
    , chunk_sz(std::move(rhs.chunk_sz))
{}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::collate_buffers()
{
    assert(Colored_);

if constexpr(Colored_)
{
    std::size_t sz = 0; // Number of pending super k-mers in the worker-local buffers.
    auto src_min = std::numeric_limits<source_id_t>::max();
    auto src_max = std::numeric_limits<source_id_t>::min();
    std::for_each(chunk_w.cbegin(), chunk_w.cend(), [&](const auto& c)
    {
        const auto& c_w = c.unwrap();
        sz += c_w.size();
        if(!c_w.empty())
        {
            src_min = std::min(src_min, c_w.front_att().source());
            src_max = std::max(src_max, c_w.back_att().source());
        }
    });

    chunk.resize(sz);
    size_ += sz;


    if(CF_UNLIKELY(src_min == src_max)) // Special case possible with large individual sources and low worker-count.
    {
        std::size_t off = 0;
        std::for_each(chunk_w.begin(), chunk_w.end(), [&](auto& c)
        {
            auto& c_w = c.unwrap();
            chunk.copy(off, c_w, 0, c_w.size());
            off += c_w.size();
            c_w.clear();
        });

        flush_chunk();
        return;
    }


    src_hist.clear();
    src_hist.resize(src_max - src_min + 1);
    std::for_each(chunk_w.cbegin(), chunk_w.cend(), [&](const auto& c)
    {
        const auto& c_w = c.unwrap();
        for(std::size_t i = 0; i < c_w.size(); ++i)
        {
            const auto src = c_w.att_at(i).source();
            assert(src_min <= src && src <= src_max);

            src_hist[src - src_min]++;
        }
    });


    auto& src_off = src_hist;   // Offsets for super k-mers in the chunk, for counting sort.
    uint64_t pref_sum = 0;
    std::for_each(src_off.begin(), src_off.end(), [&](auto& off)
    {
        const auto temp = off;
        off = pref_sum;
        pref_sum += temp;
    });
    assert(pref_sum == sz);


    std::for_each(chunk_w.begin(), chunk_w.end(), [&](auto& c)
    {
        auto& c_w = c.unwrap();
        source_id_t src = 0;
        for(std::size_t i = 0, j; i < c_w.size(); i = j)
        {
            assert(c_w.att_at(i).source() >= src);  // Ensure super k-mers are source-sorted.
            src = c_w.att_at(i).source();

            for(j = i + 1; j < c_w.size(); ++j)
                if(c_w.att_at(j).source() != src)
                    break;

            const auto stretch_sz = j - i;
            const auto src_rel = src - src_min;
            chunk.copy(src_off[src_rel], c_w, i, stretch_sz);

            src_off[src_rel] += stretch_sz;
            src < src_max ? assert(src_off[src_rel] <= src_off[src_rel + 1]):
                            assert(src_off[src_rel] <= sz);
        }

        c_w.clear();
    });

    flush_chunk();
}
}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::flush_chunk()
{
    if(!chunk.empty())
    {
        chunk.serialize(output);
        if constexpr(Colored_)
            chunk_sz.push_back(chunk.size());

        chunk.clear();
    }
}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::close()
{
    if constexpr(!Colored_)
    {
        // TODO: no need to empty or serialize the chunks—the subsequent iteration over the bucket should
        // handle the buffered in-memory content. This would skip #buckets many syscalls.
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            empty_w_local_chunk(w_id);

        flush_chunk();
    }
    else
        collate_buffers();

    output.close();
}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::remove()
{
    if(output.is_open())
        output.close();

    if(!output || !remove_file(path_))
    {
        std::cerr << "Error removing file at " << path_ << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <bool Colored_>
typename Super_Kmer_Bucket<Colored_>::Iterator Super_Kmer_Bucket<Colored_>::iterator() const
{
    assert(chunk.empty());
    return Iterator(*this);
}


template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Iterator::Iterator(const Super_Kmer_Bucket& B):
      B(B)
    , input(B.path_, std::ios::in | std::ios::binary)
    , idx(0)
    , chunk_start_idx(0)
    , chunk_end_idx(0)
    , chunk_id(0)
{}


// TODO: inline.
template <bool Colored_>
std::size_t Super_Kmer_Bucket<Colored_>::Iterator::read_chunk()
{
    assert(chunk_end_idx < B.size());
    assert(!Colored_ || chunk_id < B.chunk_sz.size());
    const auto super_kmers_to_read = (!Colored_ ?
                                        std::min(B.size() - chunk_end_idx, B.chunk.capacity()):
                                        B.chunk_sz[chunk_id++]);

    B.chunk.deserialize(input, super_kmers_to_read);

    return super_kmers_to_read;
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Bucket<false>;
template class cuttlefish::Super_Kmer_Bucket<true>;
