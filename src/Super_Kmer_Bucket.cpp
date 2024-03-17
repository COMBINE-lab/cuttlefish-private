
#include "Super_Kmer_Bucket.hpp"

#include <cassert>
#include <utility>


namespace cuttlefish
{

template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Super_Kmer_Bucket(const uint16_t k, const uint16_t l, const std::string& path):
      path_(path)
    , output(path_, std::ios::out | std::ios::binary)
    , chunk_cap(chunk_bytes / Super_Kmer_Chunk<Colored_>::record_size(k, l))
    , chunk(k, l, chunk_cap)
{
    assert(chunk_cap >= parlay::num_workers());

    const auto chunk_cap_per_w = chunk_cap / parlay::num_workers();
    chunk_w.reserve(parlay::num_workers());
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
        chunk_w.emplace_back(chunk_t(k, l, chunk_cap_per_w));
}


template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Super_Kmer_Bucket(Super_Kmer_Bucket&& rhs):
      path_(std::move(rhs.path_))
    , output(std::move(rhs.output))
    , chunk_cap(std::move(rhs.chunk_cap))
    , chunk(std::move(rhs.chunk))
    , chunk_w(std::move(rhs.chunk_w))
{}


template <bool Colored_>
void Super_Kmer_Bucket<Colored_>::close()
{
    for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
        empty_w_local_chunk(w_id);

    output.close();
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Bucket<false>;
