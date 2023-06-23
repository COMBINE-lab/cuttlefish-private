
#ifndef CONTRACTED_GRAPH_EXPANDER_HPP
#define CONTRACTED_GRAPH_EXPANDER_HPP



#include "Edge_Matrix.hpp"

#include <cstdint>
#include <cstddef>
#include <string>


namespace cuttlefish
{

// =============================================================================
// Expander for contracted discontinuity-graphs.
template <uint16_t k>
class Contracted_Graph_Expander
{
private:

    Edge_Matrix<k>& E;  // Edge matrix of the (augmented) discontinuity graph.

    const std::string work_path;    // Path-prefix to temporary working files.


public:

    // Constructs an expander for the contracted discontinuity-graph with edge-
    // matrix `E`. Temporary files are stored at path-prefix `temp_path`.
    Contracted_Graph_Expander(Edge_Matrix<k>& E, const std::string& temp_path);
};

}



#endif
