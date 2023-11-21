
#ifndef SUBGRAPH_HPP
#define SUBGRAPH_HPP



#include "State_Config.hpp"
#include "Directed_Vertex.hpp"
#include "DNA.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <unordered_map>


namespace cuttlefish
{

// =============================================================================
// A subgraph of the de Bruijn graph, induced by a KMC bin.
template <uint16_t k>
class Subgraph
{
private:

    const std::string graph_bin_dir_path;   // Path to the directory with all the graph KMC-bins.
    const std::size_t bin_id; // ID of the graph KMC-bin.

    typedef std::unordered_map<Kmer<k>, State_Config, Kmer_Hasher<k>> map_t;
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

}




#endif
