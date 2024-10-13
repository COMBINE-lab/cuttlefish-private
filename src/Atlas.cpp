
#include "Atlas.hpp"
#include "parlay/parallel.h"

#include <cstddef>


namespace cuttlefish
{

template <bool Colored_>
Atlas<Colored_>::Atlas(uint16_t k, uint16_t l, const std::string& path, std::size_t chunk_cap, std::size_t chunk_cap_per_w):
      path_(path)
    , size_(0)
    , chunk(k, l, chunk_cap)
{
    chunk_w.reserve(parlay::num_workers());
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
        chunk_w.emplace_back(chunk_t(k, l, chunk_cap_per_w));

    subgraph.reserve(graph_per_atlas());
    for(std::size_t i = 0; i < graph_per_atlas(); ++i)
        subgraph.emplace_back(k, l, path + "_G_" + std::to_string(i), subgraph_chunk_cap_bytes / chunk.record_size());
}


template <bool Colored_>
Atlas<Colored_>::Atlas(Atlas&& rhs):
      path_(std::move(rhs.path_))
    , size_(rhs.size_)
    , chunk(std::move(rhs.chunk))
    , chunk_w(std::move(rhs.chunk_w))
    , subgraph(std::move(rhs.subgraph))
{}

}



// Template-instantiations for the required instances.
template class cuttlefish::Atlas<false>;
template class cuttlefish::Atlas<true>;
