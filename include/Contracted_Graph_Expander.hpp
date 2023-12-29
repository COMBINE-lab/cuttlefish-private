
#ifndef CONTRACTED_GRAPH_EXPANDER_HPP
#define CONTRACTED_GRAPH_EXPANDER_HPP



#include "Discontinuity_Graph.hpp"
#include "Discontinuity_Edge.hpp"
#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Kmer.hpp"
#include "Kmer_Hasher.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "globals.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// Expander for contracted discontinuity-graphs.
template <uint16_t k>
class Contracted_Graph_Expander
{
private:

    const Discontinuity_Graph<k>& G;    // The (augmented) discontinuity graph.

    // TODO: consider using padding.
    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v;   // `P_v[i]` contains path-info for vertices in partition `i`.
    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e; // `P_e[b]` contains path-info for edges in bucket `b`.

    typedef std::vector<std::vector<Obj_Path_Info_Pair<Kmer<k>, k>>> worker_local_P_v_t;    // `P_v` buffers specific to a worker.
    typedef std::vector<std::vector<Obj_Path_Info_Pair<uni_idx_t, k>>> worker_local_P_e_t;  // `P_e` buffers specific to a worker.
    std::vector<Padded_Data<worker_local_P_v_t>> P_v_w; // `P_v_w[w][j]` contains path-info of vertices in partition `j` computed by worker `w`.
    std::vector<Padded_Data<worker_local_P_e_t>> P_e_w; // `P_e_w[w][b]` contains path-info of edges in bucket `b` computed by worker `w`.

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

    // Collates worker local buffers in ID range `[beg, end)` from the
    // collection `source` into the global repository `dest`. Also clears the
    // local buffers.
    template <typename T_s_, typename T_d_> static uint64_t collate_w_local_bufs(T_s_& source, std::size_t beg, size_t end, T_d_& dest);

#ifndef NDEBUG
    std::vector<Padded_Data<uint64_t>> H_p_e_w; // `H_p_e_w[w]` contains 64-bit hash of the path-info of edges computed by worker `w`.
#endif

    // Debug
    double p_v_load_time = 0;   // Time to load vertices' path-info.
    double edge_read_time = 0;  // Time taken to read the edges.
    double map_fill_time = 0;   // Time taken to fill the hash map with already inferred vertices' path-info for a partition.

    static constexpr auto now = std::chrono::high_resolution_clock::now;    // Current time-point in nanoseconds.

    // Returns the equivalent time-duration in seconds from `d` nanoseconds.
    static constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };



public:

    // Constructs an expander for the contracted discontinuity-graph `G`.
    // `P_v[i]` is to contain path-information for vertices at partition `i`,
    // and `P_e[b]` is to contain path-information for edges at bucket `b`.
    // Temporary files are stored at path-prefix `temp_path`.
    Contracted_Graph_Expander(const Discontinuity_Graph<k>& G, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const std::string& temp_path);

    // Expands the contracted discontinuity-graph.
    void expand();
};


template <uint16_t k>
inline Path_Info<k> Contracted_Graph_Expander<k>::infer(const Path_Info<k> u_inf, const side_t s_u, const side_t s_v, const weight_t w)
{
    // const auto p_v = u_inf.p(); // Path-ID.
    const auto r_v = (s_u == u_inf.o() ? u_inf.r() + w :    // Rank.
                                        (u_inf.r() > w ? u_inf.r() - w : 0));   // Trying to expand crossing a deleted edge from an ICC.
                                                                                // This works as no vertex can have a rank `0` in the model.
    const auto o_v = (s_u == u_inf.o() ? inv_side(s_v) : s_v);  // Orientation.

    return Path_Info(u_inf.p(), r_v, o_v);
}


template <uint16_t k>
inline void Contracted_Graph_Expander<k>::add_edge_path_info(const Discontinuity_Edge<k>& e, const Path_Info<k> u_inf, const Path_Info<k> v_inf)
{
    // const auto p = u_inf.p();
    assert(!e.x_is_phi() && !e.y_is_phi());
    assert(u_inf.p() == v_inf.p());

    const auto r = std::min(u_inf.r(), v_inf.r());
    const auto o = (r > 0 ? (r == u_inf.r() ? e.o() : inv_side(e.o())) :    // Whether the ranking of the vertices in the unitig goes from `u` to `v`.
                            side_t::unspecified);   // This is a deleted edge from an ICC.

    assert(e.b() > 0 && e.b() < P_e.size());
    P_e_w[parlay::worker_id()].data()[e.b()].emplace_back(e.b_idx(), u_inf.p(), r, o);

#ifndef NDEBUG
    H_p_e_w[parlay::worker_id()].data() ^= P_e_w[parlay::worker_id()].data()[e.b()].back().path_info().hash();
#endif
}


template <uint16_t k>
inline void Contracted_Graph_Expander<k>::add_edge_path_info(const Discontinuity_Edge<k>& e, const Path_Info<k> v_inf)
{
    // const auto p = v_inf.p();
    assert(e.x_is_phi() && !e.y_is_phi());

    const auto r = (v_inf.r() == 1 ? 0 : v_inf.r());
    const auto o = (r == 0 ? e.o() : inv_side(e.o()));

    assert(e.b() > 0 && e.b() < P_e.size());
    P_e_w[parlay::worker_id()].data()[e.b()].emplace_back(e.b_idx(), v_inf.p(), r, o);

#ifndef NDEBUG
    H_p_e_w[parlay::worker_id()].data() ^= P_e_w[parlay::worker_id()].data()[e.b()].back().path_info().hash();
#endif
}

}



#endif
