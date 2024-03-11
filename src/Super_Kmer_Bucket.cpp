
#include "Super_Kmer_Bucket.hpp"

#include <cassert>


namespace cuttlefish
{

template <bool Colored_>
Super_Kmer_Bucket<Colored_>::Super_Kmer_Bucket(const uint16_t k, const uint16_t l, const std::string& path):
      path_(path)
    , chunk_cap(chunk_bytes / Super_Kmer_Chunk<Colored_>::record_size(k, l))
    , chunk(k, l, chunk_cap)
    , chunk_w(allocate<Padded_Data<chunk_t>>(parlay::num_workers()))
{
    assert(chunk_cap >= parlay::num_workers());

    const auto chunk_cap_per_w = chunk_cap / parlay::num_workers();
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
        new (chunk_w + i) Padded_Data(chunk_t(k, l, chunk_cap_per_w));
}


template <bool Colored_>
Super_Kmer_Bucket<Colored_>::~Super_Kmer_Bucket()
{
    deallocate(chunk_w);
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Bucket<false>;
