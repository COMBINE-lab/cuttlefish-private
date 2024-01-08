
#ifndef STATE_CONFIG_HPP
#define STATE_CONFIG_HPP



#include "DNA.hpp"
#include "globals.hpp"
#include "utility.hpp"

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

    static uint8_t f_th;    // Edge-frequency threshold.

    // uint64_t color_hash;        // Hash of the vertex's color-set.
    uint8_t edge_freq[8];   // Frequency of the vertex's neighbors.
    // uint32_t last_color_ID;     // Last color-ID added to the color-hash; to encounter multi-set hashing problem for color-sets.
    uint8_t status;             // Some status information of the vertex, bit-packed:
                                    // whether it is a discontinuity vertex, whether it's been visited.

    static constexpr uint8_t visited = 0b0000'0001; // Flag to denote a vertex as visited.
    static constexpr uint8_t discontinuity[2] = {0b0000'0010, 0b0000'0100}; // Flags to denote a vertex's sides as discontinuous.


public:

    // Constructs an empty state.
    State_Config();

    // Sets the edge-frequency threshold to `f_th`.
    static void set_edge_threshold(uint8_t f_th);

    // Adds the edge-encodings `front` and `back` to the associated sides of a
    // corresponding vertex.
    void update_edges(base_t front, base_t back);

    // Adds the edge-encoding `b` to the side `s` of a corresponding vertex.
    void update_edge(side_t s, base_t b);

    // Marks the associated vertex as visited.
    void mark_visited() { status |= visited; }

    // Marks the associated vertex as discontinuous at side `s`.
    void mark_discontinuous(side_t s);

    // Returns whether the associated vertex is visited.
    bool is_visited() const { return status & visited; }

    // Returns whether the associated vertex is discontinuous at side `s`.
    bool is_discontinuous(const side_t s) const;

    // Returns whether the associated vertex has any discontinuous side.
    bool is_discontinuity() const { return is_discontinuous(side_t::front) | is_discontinuous(side_t::back); }

    // Returns the `Base`-encoding of the edge(s) incident to the side `s` of a
    // vertex having this state.
    base_t edge_at(side_t s) const;

    // Returns `true` iff some vertex having this state is branching (i.e. has
    // multiple incident edges) at its side `s`.
    bool is_branching_side(side_t s) const;

    // Returns `true` iff some vertex having this state is empty (i.e. no
    // incident edges) at its side `s`.
    bool is_empty_side(side_t s) const;

    // Returns `true` iff some vertex having this state is isolated off the rest
    // of the underlying graph.
    bool is_isolated() const;
};


inline State_Config::State_Config():
    //   color_hash(0)
      edge_freq()
    // , last_color_ID(0)
    , status(0)
{}


inline void State_Config::update_edges(const base_t front, const base_t back)
{
    constexpr uint8_t max_f = (1lu << CHAR_BIT) - 1;    // Maximum supported frequency of a (k + 1)-mer.
    constexpr std::size_t back_off = 4;
    constexpr auto E = base_t::E;
    constexpr auto T = base_t::T;
    (void)T;

    assert(front == E || front <= T);
    if(front != E && edge_freq[front] < max_f)
        edge_freq[front]++;

    assert(back == E || back <= T);
    if(back != E && edge_freq[back_off + back] < max_f)
        edge_freq[back_off + back]++;
}


inline void State_Config::update_edge(const side_t s, const base_t b)
{
    constexpr uint8_t max_f = (1lu << CHAR_BIT) - 1;    // Maximum supported frequency of a (k + 1)-mer.
    constexpr auto E = base_t::E;
    constexpr auto T = base_t::T;
    (void)T;

    assert(b == E || b <= T);
    const std::size_t off = (s == side_t::back) * 4;
    if(b != E && edge_freq[off + b] < max_f)
        edge_freq[off + b]++;
}


inline void State_Config::mark_discontinuous(const side_t s)
{
    assert(as_int(s) < 2);
    status |= discontinuity[as_int(s)];
}


inline bool State_Config::is_discontinuous(const side_t s) const
{
    assert(as_int(s) < 2);
    return status & discontinuity[as_int(s)];
}


inline base_t State_Config::edge_at(const side_t s) const
{
    const std::size_t off = (s == side_t::back) * 4;    // Offset to the encoding of edge frequencies at side `s`.
    const auto edge_c = (edge_freq[off + 0] >= f_th) + (edge_freq[off + 1] >= f_th) +
                        (edge_freq[off + 2] >= f_th) + (edge_freq[off + 3] >= f_th);

    switch(edge_c)
    {
        case 0:
            return base_t::E;

        case 1:
            return static_cast<base_t>( ((edge_freq[off + 0] >= f_th) * base_t::A) +
                                        ((edge_freq[off + 1] >= f_th) * base_t::C) +
                                        ((edge_freq[off + 2] >= f_th) * base_t::G) +
                                        ((edge_freq[off + 3] >= f_th) * base_t::T));

        default:
            return base_t::N;
    }
}


inline bool State_Config::is_branching_side(const side_t s) const
{
    const std::size_t off = (s == side_t::back) * 4;    // Offset to the encoding of edge frequencies at side `s`.
    const auto edge_c = (edge_freq[off + 0] >= f_th) + (edge_freq[off + 1] >= f_th) +
                        (edge_freq[off + 2] >= f_th) + (edge_freq[off + 3] >= f_th);

    return edge_c > 1;
}


inline bool State_Config::is_empty_side(const side_t s) const
{
    const std::size_t off = (s == side_t::back) * 4;    // Offset to the encoding of edge frequencies at side `s`.
    const auto edge_c = (edge_freq[off + 0] >= f_th) + (edge_freq[off + 1] >= f_th) +
                        (edge_freq[off + 2] >= f_th) + (edge_freq[off + 3] >= f_th);

    return edge_c == 0;
}


inline bool State_Config::is_isolated() const
{
    return is_empty_side(side_t::back) && is_empty_side(side_t::front);
}

}



#endif
