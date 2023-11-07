
#ifndef SUBGRAPH_HPP
#define SUBGRAPH_HPP



#include "Directed_Vertex.hpp"

#include <cstdint>
#include <cstddef>
#include <string>


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


public:

    // Constructs a subgraph object for the `bin_id`'th bin in the graph bin
    // directory `bin_dir_path`.
    Subgraph(const std::string& bin_dir_path, std::size_t bin_id);

    Subgraph(const Subgraph&) = delete;
    Subgraph(Subgraph&&) = delete;

    // Loads the subgraph into an internal navigable and membership data
    // structure.
    void load();
};

}




#endif
