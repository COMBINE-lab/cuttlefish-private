
#ifndef ATLAS_HPP
#define ATLAS_HPP



#include "utility.hpp"

#include <cstdint>


namespace cuttlefish
{

static constexpr uint64_t atlas_count = 128;    // Number of subgraph atlases.
static constexpr uint64_t graph_per_atlas = 128;    // Number of subgraphs per atlas.
static constexpr uint64_t graph_count_ = atlas_count * graph_per_atlas; // Number of subgraphs.


// Returns the atlas-ID of the `g`'th subgraph.
inline auto atlas_ID(const uint64_t g) { return g >> log_2(atlas_count); }

}



#endif
