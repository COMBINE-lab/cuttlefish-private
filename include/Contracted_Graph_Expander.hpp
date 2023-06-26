
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

    std::vector<Discontinuity_Edge<k>> D_i; // New edges introduced in contracted diagonal blocks.

    std::unordered_map<Kmer<k>, Vertex_Path_Info<k>, Kmer_Hasher<k>> M; // `M[v]` is the path-info for vertex `v`.


    // Expands the `[i, i]`'th (contracted) edge-block.
    void expand_diagonal_block(std::size_t i);

    // Infers a vertex v's path-info from that of vertex u's path-info `u_inf`
    // and returns it. The vertices are connected with an edge of weight `w`
    // through their sides `s_v` and `s_u`, respectively.
    Vertex_Path_Info<k> infer(Vertex_Path_Info<k> u_inf, side_t s_u, side_t s_v, weight_t w);


public:

    // Constructs an expander for the contracted discontinuity-graph with edge-
    // matrix `E`. `P_v[i]` is to contain path-information for vertices at
    // partition `i`. Temporary files are stored at path-prefix `temp_path`.
    Contracted_Graph_Expander(Edge_Matrix<k>& E, std::vector<Ext_Mem_Bucket<Vertex_Path_Info_Pair<k>>>& P_v, const std::string& temp_path);
};


template <uint16_t k>
Vertex_Path_Info<k> Contracted_Graph_Expander<k>::infer(const Vertex_Path_Info<k> u_inf, const side_t s_u, const side_t s_v, const weight_t w)
{
    // const auto p_v = u_inf.p(); // Path-ID.
    const weight_t r_v = (s_u == u_inf.o() ? u_inf.r() + w : u_inf.r() - w);    // Rank.
    const side_t o_v = (s_u == u_inf.o() ? inv_side(s_v) : s_v);    // Orientation.

    return Vertex_Path_Info(u_inf.p(), r_v, o_v);
}

}



#endif