
#ifndef SUBGRAPHS_PROCESSOR_HPP
#define SUBGRAPHS_PROCESSOR_HPP



#include "Discontinuity_Graph.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Character_Buffer.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>


namespace cuttlefish
{

// =============================================================================
// Processor for the subgraphs of the de Bruijn graph—constructs them from their
// super k-mer based sequence representation in KMC bins, and contracts them
// into their compacted form.
template <uint16_t k>
class Subgraphs_Processor
{
private:

    const std::string bin_path_pref;    // Path-prefix to the KMC bins.
    const std::size_t bin_count;    // Numer of KMC bins, each one induces a subgraph.

    const std::string work_path;    // Path to the working space.

    Discontinuity_Graph<k>& G;  // The discontinuity graph.

    uint64_t trivial_mtig_count_;   // Number of trivial maximal unitigs in the subgraphs (i.e. also maximal unitigs in the supergraph).
    uint64_t icc_count_;    // Number of trivial maximal unitigs in the subgraphs that are ICCs.

    // TODO: move out the following to some CF3-centralized location.

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded_Data<op_buf_t>> op_buf_list_t;
    op_buf_list_t& op_buf; // Worker-specific output buffers.


public:

    // Constructs a processor for the subgraphs induced by the `bin_count` KMC
    // bins at path-prefix `bin_path_pref`. The discontinuity-graph is produced
    // at `G` without false-phantom edges. Worker-specific trivially maximal
    // unitigs are written to the buffers in `op_buf`.
    Subgraphs_Processor(const std::string& bin_path_pref, std::size_t bin_count, Discontinuity_Graph<k>& G, op_buf_list_t& op_buf);

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
