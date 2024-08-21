
#include "Super_Kmer_Bucket.hpp"
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
    , safe_src_lower_bound(1)
    , safe_src_upper_bound(1)
    , safe_c(0)
    , chunk_safe(k, l, !Colored_ ? 0 : chunk_cap)
{
    assert(chunk_cap >= parlay::num_workers());

    const auto chunk_cap_per_w = w_chunk_bytes / Super_Kmer_Chunk<Colored_>::record_size(k, l);
    chunk_w.reserve(parlay::num_workers());
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
        chunk_w.emplace_back(chunk_t(k, l, chunk_cap_per_w));

    if constexpr(Colored_)
        latest_src_w.resize(parlay::num_workers(), safe_src_lower_bound);
}


template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Super_Kmer_Bucket(Super_Kmer_Bucket&& rhs):
      path_(std::move(rhs.path_))
    , output(std::move(rhs.output))
    , size_(std::move(rhs.size_))
    , chunk_cap(std::move(rhs.chunk_cap))
    , chunk(std::move(rhs.chunk))
    , chunk_w(std::move(rhs.chunk_w))
    , latest_src_w(std::move(rhs.latest_src_w))
    , safe_src_lower_bound(std::move(rhs.safe_src_lower_bound))
    , safe_src_upper_bound(std::move(rhs.safe_src_upper_bound))
    , safe_c(std::move(rhs.safe_c))
    , src_hist(std::move(rhs.src_hist))
    , chunk_safe(std::move(rhs.chunk_safe))
    , chunk_sz(std::move(rhs.chunk_sz))
{}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::close()
{

    if constexpr(!Colored_)
        // TODO: no need to empty or serialize the chunks—the subsequent iteration over the bucket should
        // handle the buffered in-memory content. This would skip #buckets many syscalls.
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            empty_w_local_chunk(w_id);
    else
    {
        uint32_t max_src = 0;   // Maximum source-ID of all super k-mers deposited to the bucket.
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
        {
            const auto& c_w = chunk_w[w_id].unwrap();
            max_src = std::max(max_src, !c_w.empty() ? c_w.back_att().source() : latest_src_w[w_id].unwrap());
        }

        if(src_hist.size() < rel_src(max_src) + 1)
            src_hist.resize(rel_src(max_src) + 1);


        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
        {
            auto& c_w = chunk_w[w_id].unwrap();
            if(!c_w.empty())
            {
                chunk.reserve(chunk.size() + c_w.size());
                chunk.append(c_w);
                size_ += c_w.size();

                for(std::size_t i = 0; i < c_w.size(); ++i)
                    src_hist[rel_src(c_w.att_at(i).source())]++;

                c_w.clear();
            }
        }

        assert(max_src >= safe_src_upper_bound);
        safe_src_upper_bound = max_src;
    }


    if constexpr(!Colored_)
    {
        if(!chunk.empty())
        {
            chunk.serialize(output);
            chunk.clear();
        }
    }
    else
    {
        safe_c = chunk.size();
        if(!chunk.empty())
            counting_sort_safe_super_kmers();
        assert(chunk.empty());

        if(!chunk_safe.empty())
        {
            chunk_safe.serialize(output);
            chunk_sz.push_back(chunk_safe.size());
            chunk_safe.clear();
        }
    }

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

    auto* const c = B.reader_chunk();
    c->deserialize(input, super_kmers_to_read);

    return super_kmers_to_read;
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Bucket<false>;
template class cuttlefish::Super_Kmer_Bucket<true>;
