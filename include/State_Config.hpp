
#ifndef STATE_CONFIG_HPP
#define STATE_CONFIG_HPP



#include "DNA.hpp"

#include <cstdint>
#include <climits>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// Class for a full state-configuration of a vertex in a de Bruijn graph: this
// is a configuration attached to vertices in subgraphs.
class State_Config
{
private:

    // uint64_t color_hash;        // Hash of the vertex's color-set.
    uint8_t neighbor_freq[8];   // Frequency of the vertex's neighbors.
    // uint32_t last_color_ID;     // Last color-ID added to the color-hash; to encounter multi-set hashing problem for color-sets.
    uint8_t status;             // Some status information of the vertex, bit-packed:
                                    // whether it is a discontinuity vertex, whether it's been visited.

    static constexpr uint8_t discontinuity  = 0b0000'0001;  // Flag to denote a vertex as a discontinuity one.
    static constexpr uint8_t visited        = 0b0000'0010;  // Flag to denote a vertex as visited.

public:

    // Constructs an empty state.
    State_Config();

    // Adds the neighbor-encodings `front` and `back` to the associated sides of
    // of a corresponding vertex.
    void add_neighbor(DNA::Base front, DNA::Base back);

    // Returns whether the associated vertex is a discontinuity.
    bool is_discontinuity() const { return status & discontinuity; }

    // Returns whether the associated vertex is visited.
    bool is_visited() const { return status & visited; }
};


inline State_Config::State_Config():
    //   color_hash(0)
      neighbor_freq()
    // , last_color_ID(0)
    , status(0)
{}


inline void State_Config::add_neighbor(const DNA::Base front, const DNA::Base back)
{
    constexpr uint8_t max_f = (1lu << CHAR_BIT) - 1;    // Maximum supported frequency of a (k + 1)-mer.
    constexpr std::size_t back_off = 4;

    if(front != DNA::Base::N && neighbor_freq[front] < max_f)
        assert(front <= DNA::Base::T),
        neighbor_freq[front]++;

    if(back != DNA::Base::N && neighbor_freq[back_off + back] < max_f)
        assert(back <= DNA::Base::T),
        neighbor_freq[back_off + back]++;
}

}



#endif
