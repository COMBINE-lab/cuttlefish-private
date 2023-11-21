
#ifndef SUBGRAPH_HPP
#define SUBGRAPH_HPP



#include "Directed_Vertex.hpp"
#include "DNA.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <unordered_map>


namespace cuttlefish
{


class Vertex_Info;


// =============================================================================
// A subgraph of the de Bruijn graph, induced by a KMC bin.
template <uint16_t k>
class Subgraph
{
private:

    const std::string graph_bin_dir_path;   // Path to the directory with all the graph KMC-bins.
    const std::size_t bin_id; // ID of the graph KMC-bin.

    typedef std::unordered_map<Kmer<k>, Vertex_Info, Kmer_Hasher<k>> map_t;
    map_t M;

    uint64_t edge_c;    // Number of edges in the graph.


public:

    // Constructs a subgraph object for the `bin_id`'th bin in the graph bin
    // directory `bin_dir_path`.
    Subgraph(const std::string& bin_dir_path, std::size_t bin_id);

    Subgraph(const Subgraph&) = delete;
    Subgraph(Subgraph&&) = delete;

    // Constructs the subgraph from the KMC bin into an internal navigable and
    // membership data structure.
    void construct();

    // Builds the compacted graph from the original graph.
    void compact();

    // Returns the size of the graph.
    std::size_t size() const;

    // Returns the number of (multi-)edges in the graph.
    std::size_t edge_count() const;
};


// Information associated to a vertex in a subgraph.
class Vertex_Info
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

    // Constructs an empty information object.
    Vertex_Info();

    // Adds the neighbor-encodings `front` and `back` to the associated sides of
    // of the corresponding vertex.
    void add_neighbor(DNA::Base front, DNA::Base back);

    // Returns whether the associated vertex is a discontinuity.
    bool is_discontinuity() const { return status & discontinuity; }

    // Returns whether the associated vertex is visited.
    bool is_visited() const { return status & visited; }
};


inline Vertex_Info::Vertex_Info():
    //   color_hash(0)
      neighbor_freq()
    // , last_color_ID(0)
    , status(0)
{}


inline void Vertex_Info::add_neighbor(const DNA::Base front, const DNA::Base back)
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
