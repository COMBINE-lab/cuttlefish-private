
#ifndef SUBGRAPHS_MANAGER_HPP
#define SUBGRAPHS_MANAGER_HPP



#include "Discontinuity_Graph.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Character_Buffer.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Manager for the subgraphs of the de Bruijn graph—manages their super k-mer
// based sequence representations in buckets, constructs them from these
// representations, and contracts them into their compacted form.
template <uint16_t k>
class Subgraphs_Manager
{
private:

    const std::string path_pref;    // Path-prefix to the super k-mer buckets.
    const std::size_t graph_count;  // Numer of subgraphs.

    Discontinuity_Graph<k>& G;  // The discontinuity graph.

    uint64_t trivial_mtig_count_;   // Number of trivial maximal unitigs in the subgraphs (i.e. also maximal unitigs in the supergraph).
    uint64_t icc_count_;    // Number of trivial maximal unitigs in the subgraphs that are ICCs.

    // TODO: move out the following to some CF3-centralized location.

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded_Data<op_buf_t>> op_buf_list_t;
    op_buf_list_t& op_buf; // Worker-specific output buffers.


public:

    // Constructs a manager for the subgraphs of a de Bruijn graph which is
    // partitioned into `graph_count` subgraphs. `logistics` is the data
    // logistics manager for the algorithm execution. The discontinuity-graph
    // is produced at `G` without false-phantom edges. Worker-specific
    // trivially maximal unitigs are written to the buffers in `op_buf`.
    Subgraphs_Manager(const Data_Logistics& logistics, std::size_t graph_count, Discontinuity_Graph<k>& G, op_buf_list_t& op_buf);

    // Constructs and contracts each subgraph.
    void process();

    // Returns the number of trivial maximal unitigs in the subgraphs (i.e. also
    // maximal unitigs in the supergraph).
    uint64_t trivial_mtig_count() const;

    // Returns the number of trivial maximal unitigs in the subgraphs that are
    // ICCs.
    uint64_t icc_count() const;
};

}



#endif
