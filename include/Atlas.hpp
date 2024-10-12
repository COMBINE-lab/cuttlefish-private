
#ifndef ATLAS_HPP
#define ATLAS_HPP



#include "utility.hpp"

#include <cstdint>


namespace cuttlefish
{

template <bool Colored_>
class Atlas
{
private:

    static constexpr uint64_t atlas_count_ = 128;   // Number of subgraph atlases.
    static constexpr uint64_t graph_per_atlas_ = 128;    // Number of subgraphs per atlas.
    static constexpr uint64_t graph_count_ = atlas_count_ * graph_per_atlas_;   // Number of subgraphs.

public:

    // Returns the number of subgraph atlases.
    static constexpr auto atlas_count() { return atlas_count_; }

    // Returns the number of subgraphs per atlas.
    static constexpr auto graph_per_atlas() { return graph_per_atlas_; }

    // Returns the number of subgraphs.
    static constexpr auto graph_count() { return graph_count_; }

    // Returns the atlas-ID of the `g`'th subgraph.
    static auto atlas_ID(const uint64_t g) { return g >> log_2(graph_per_atlas()); }

    // Returns the graph-ID of the `g`'th subgraph within its atlas.
    static auto graph_ID(const uint64_t g) { return g & (graph_per_atlas() - 1); }
};

}



#endif
