
#ifndef SUBGRAPHS_PROCESSOR_HPP
#define SUBGRAPHS_PROCESSOR_HPP



#include "Edge_Matrix.hpp"

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

    const std::string& bin_path_pref;   // Path-prefix to the KMC bins.
    const std::size_t bin_count;    // Numer of KMC bins, each one induces a subgraph.

    Edge_Matrix<k>& E;  // Edge-matrix of the discontinuity graph.


public:

    // Constructs a processor for the subgraphs induced by the `bin_count` KMC
    // bins at path-prefix `bin_path_pref`. Edge-matrix of the discontinuity-
    // graph is produced at `E`.
    Subgraphs_Processor(const std::string& bin_path_pref, std::size_t bin_count, Edge_Matrix<k>& E);

    // Constructs and contracts each subgraph.
    void process();
};

}



#endif
