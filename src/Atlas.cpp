
#include "Atlas.hpp"
#include "parlay/parallel.h"

#include <cstddef>
#include <limits>


namespace cuttlefish
{

template <bool Colored_>
Atlas<Colored_>::Atlas(uint16_t k, uint16_t l, const std::string& path, std::size_t chunk_cap, std::size_t chunk_cap_per_w):
      path_(path)
    , size_(0)
    , chunk_cap(chunk_cap)
    , chunk(new chunk_t(k, l, chunk_cap))
    , flush_buf(!Colored_ ? new chunk_t(k, l, chunk_cap) : nullptr)
{
    chunk_w.reserve(parlay::num_workers());
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
        chunk_w.emplace_back(chunk_t(k, l, chunk_cap_per_w));

    subgraph.reserve(graph_per_atlas());
    for(std::size_t i = 0; i < graph_per_atlas(); ++i)
        subgraph.emplace_back(k, l, path + "_G_" + std::to_string(i), subgraph_chunk_cap_bytes / chunk->record_size());
}


template <bool Colored_>
Atlas<Colored_>::Atlas(Atlas&& rhs):
      path_(std::move(rhs.path_))
    , size_(rhs.size_)
    , chunk_cap(rhs.chunk_cap)
    , chunk(std::move(rhs.chunk))
    , flush_buf(std::move(rhs.flush_buf))
    , chunk_w(std::move(rhs.chunk_w))
    , subgraph(std::move(rhs.subgraph))
{}


template <bool Colored_>
void Atlas<Colored_>::flush_chunk(chunk_t& c)
{
    if(!c.empty())
    {
        for(std::size_t g = 0; g < graph_per_atlas(); ++g)
            subgraph[g].fetch_end();

        for(std::size_t i = 0; i < c.size(); ++i)
        {
            const auto g = graph_ID(c.att_at(i).g_id());
            subgraph[g].add(c.label_at(i), c.att_at(i));
        }

        c.clear();
    }
}


template <>
void Atlas<true>::flush_collated()
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

    if(CF_UNLIKELY(sz == 0))
        return;

    chunk->resize_uninit(sz);
    size_ += sz;


    if(CF_UNLIKELY(src_min == src_max)) // Special case possible with large individual sources and low worker-count.
    {
        std::size_t off = 0;
        std::for_each(chunk_w.begin(), chunk_w.end(), [&](auto& c)
        {
            auto& c_w = c.unwrap();
            chunk->copy(off, c_w, 0, c_w.size());
            off += c_w.size();
            c_w.clear();
        });

        flush_chunk(*chunk);
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


    auto& src_off = src_hist;   // Offsets for super k-mers in counting-sorted chunk.
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
            chunk->copy(src_off[src_rel], c_w, i, stretch_sz);

            src_off[src_rel] += stretch_sz;
            src < src_max ? assert(src_off[src_rel] <= src_off[src_rel + 1]):
                            assert(src_off[src_rel] <= sz);
        }

        c_w.clear();
    });

    flush_chunk(*chunk);
}


template <bool Colored_>
void Atlas<Colored_>::close()
{
    if constexpr(!Colored_)
    {
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            empty_w_local_chunk(w_id);

        flush_chunk(*chunk);
    }

    parlay::parallel_for(0, graph_per_atlas(),
    [&](const auto g)
    {
        subgraph[g].close();
    }, 1);
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Atlas<false>;
template class cuttlefish::Atlas<true>;
