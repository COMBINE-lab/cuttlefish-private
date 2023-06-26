
#ifndef CONTRACTED_GRAPH_EXPANDER_HPP
#define CONTRACTED_GRAPH_EXPANDER_HPP



#include "Edge_Matrix.hpp"
#include "Vertex_Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
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

    std::vector<Ext_Mem_Bucket<Vertex_Path_Info_Pair<k>>>& P_v; // `P_v[i]` contains path-info for vertices in partition `i`.

    const std::string work_path;    // Path-prefix to temporary working files.


public:

    // Constructs an expander for the contracted discontinuity-graph with edge-
    // matrix `E`. `P_v[i]` is to contain path-information for vertices at
    // partition `i`. Temporary files are stored at path-prefix `temp_path`.
    Contracted_Graph_Expander(Edge_Matrix<k>& E, std::vector<Ext_Mem_Bucket<Vertex_Path_Info_Pair<k>>>& P_v, const std::string& temp_path);
};

}



#endif
