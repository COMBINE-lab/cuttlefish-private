
#ifndef CONTRACTED_GRAPH_EXPANDER_HPP
#define CONTRACTED_GRAPH_EXPANDER_HPP



#include "Edge_Matrix.hpp"
#include "Discontinuity_Edge.hpp"
#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Kmer.hpp"
#include "Kmer_Hasher.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// Expander for contracted discontinuity-graphs.
template <uint16_t k>
class Contracted_Graph_Expander
{
    typedef uint32_t uni_idx_t; // Type of the index of a unitig in a bucket.

private:

    const Edge_Matrix<k>& E;    // Edge matrix of the (augmented) discontinuity graph.
    const std::size_t n_;   // Number of discontinuity-vertices.

    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v;   // `P_v[i]` contains path-info for vertices in partition `i`.
    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e; // `P_e[b]` contains path-info for edges in bucket `b`.

    const std::string work_path;    // Path-prefix to temporary working files.

    // TODO: remove `D_i` by adopting a more parallelization-amenable algorithm for diagonal contraction-expansion.
    std::vector<Discontinuity_Edge<k>> D_i; // New edges introduced in contracted diagonal blocks.

    Concurrent_Hash_Table<Kmer<k>, Path_Info<k>, Kmer_Hasher<k>> M; // `M[v]` is the path-info for vertex `v`.

    std::vector<Obj_Path_Info_Pair<Kmer<k>, k>> p_v_buf;    // Buffer to read-in path-information of vertices.


    // Loads the available path-info of meta-vertices from partition `i` into
    // the hash table `M`.
    void load_path_info(std::size_t i);

    // Expands the `[i, i]`'th (contracted) edge-block.
    void expand_diagonal_block(std::size_t i);

    // Infers a vertex v's path-info from that of vertex u's path-info `u_inf`
    // and returns it. The vertices are connected with an edge of weight `w`
    // through their sides `s_v` and `s_u`, respectively.
    Path_Info<k> infer(Path_Info<k> u_inf, side_t s_u, side_t s_v, weight_t w);

    // Computes the path-info of the edge `e` from its endpoints' path-info,
    // `u_inf` and `v_inf`, and adds the info to `e`'s path-info bucket.
    void add_edge_path_info(const Discontinuity_Edge<k>& e, Path_Info<k> u_inf, Path_Info<k> v_inf);

    // Computes the path-info of the edge `e` of form `(ϕ, v)` from `v`'s path-
    // info, `v_inf`, and adds the info to `e`'s path-info bucket.
    void add_edge_path_info(const Discontinuity_Edge<k>& e, Path_Info<k> v_inf);

    // Debug
    std::size_t og_edge_c = 0;


public:

    // Constructs an expander for the contracted discontinuity-graph with edge-
    // matrix `E` and `n` vertices. `P_v[i]` is to contain path-information for
    // vertices at partition `i`, and `P_e[b]` is to contain path-information
    // for edges at bucket `b`. Temporary files are stored at path-prefix
    // `temp_path`.
    Contracted_Graph_Expander(const Edge_Matrix<k>& E, std::size_t n, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const std::string& temp_path);

    // Expands the contracted discontinuity-graph.
    void expand();
};


template <uint16_t k>
inline Path_Info<k> Contracted_Graph_Expander<k>::infer(const Path_Info<k> u_inf, const side_t s_u, const side_t s_v, const weight_t w)
{
    // const auto p_v = u_inf.p(); // Path-ID.
    const weight_t r_v = (s_u == u_inf.o() ? u_inf.r() + w : u_inf.r() - w);    // Rank.
    const side_t o_v = (s_u == u_inf.o() ? inv_side(s_v) : s_v);    // Orientation.

    return Path_Info(u_inf.p(), r_v, o_v);
}


template <uint16_t k>
inline void Contracted_Graph_Expander<k>::add_edge_path_info(const Discontinuity_Edge<k>& e, const Path_Info<k> u_inf, const Path_Info<k> v_inf)
{
    // const auto p = u_inf.p();
    assert(!e.x_is_phi() && !e.y_is_phi());
    assert(u_inf.p() == v_inf.p());

    const auto r = std::min(u_inf.r(), v_inf.r());
    const auto o = (r == u_inf.r() ? e.o() : inv_side(e.o()));  // The ranking of the vertices in the unitig goes from u to v.

    assert(e.b() < P_e.size());
    P_e[e.b()].emplace(e.b_idx(), u_inf.p(), r, o); // `u_inf.p() == v_inf.p()`
}


template <uint16_t k>
inline void Contracted_Graph_Expander<k>::add_edge_path_info(const Discontinuity_Edge<k>& e, const Path_Info<k> v_inf)
{
    // const auto p = v_inf.p();
    assert(e.x_is_phi() && !e.y_is_phi());

    const auto r = (v_inf.r() == 1 ? 0 : v_inf.r());
    const auto o = (r == 0 ? e.o() : inv_side(e.o()));

    assert(e.b() < P_e.size());
    P_e[e.b()].emplace(e.b_idx(), v_inf.p(), r, o);
}

}



#endif
