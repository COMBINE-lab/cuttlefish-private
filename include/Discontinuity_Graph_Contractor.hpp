
#ifndef DISCONTINUITY_GRAPH_CONTRACTOR_HPP
#define DISCONTINUITY_GRAPH_CONTRACTOR_HPP



#include "Kmer.hpp"
#include "Kmer_Hasher.hpp"
#include "Discontinuity_Edge.hpp"
#include "Edge_Matrix.hpp"
#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Spin_Lock.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <unordered_map>
#include <tuple>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// Contractor of discontinuity-graphs.
template <uint16_t k>
class Discontinuity_Graph_Contractor
{
private:

    Edge_Matrix<k>& E;  // Edge matrix of the discontinuity-graph.
    const std::size_t n_;   // Number of discontinuity-vertices.

    typedef Obj_Path_Info_Pair<Kmer<k>, k> kmer_path_info_t;
    std::vector<Ext_Mem_Bucket<kmer_path_info_t>>& P_v; // `P_v[j]` contains path-info for vertices in partition `j`—specifically, the meta-vertices.
    std::vector<Padded_Data<std::vector<kmer_path_info_t>>> P_v_local;  // `P_v_local[t]` contains information of the meta-vertices formed by worker `t`.

    const std::string work_path;    // Path-prefix to temporary working files.

    std::vector<Discontinuity_Edge<k>> buf; // Buffer to read-in edges from the edge-matrix.

    class Other_End;
    Concurrent_Hash_Table<Kmer<k>, Other_End, Kmer_Hasher<k>> M;    // `M[v]` is the associated vertex to `v` at a given time.

    // TODO: remove `D_j` by adopting a more parallelization-amenable algorithm for diagonal contraction-expansion.
    std::vector<Discontinuity_Edge<k>> D_j; // Edges introduced in contracting a diagonal block.
    std::vector<Padded_Data<std::vector<Discontinuity_Edge<k>>>> D_c;   // `D_c[t]` contains the edges corresponding to compressed diagonal chains by worker `t`.


    // Contracts the `[j, j]`'th edge-block.
    void contract_diagonal_block(std::size_t j);

    // Forms a meta-vertex in the contracted graph with the vertex `v` belonging
    // to the vertex-partition `part`. In the contracted graph, `v` has a `w_1`
    // weighted edge incident to its side `s_1` and a `w_2` weighted edge
    // incident to the other side.
    void form_meta_vertex(Kmer<k> v, std::size_t part, side_t s_1, weight_t w_1, weight_t w_2);

    // Debug
    std::size_t meta_v_c = 0;

    static constexpr auto now = std::chrono::high_resolution_clock::now;    // Current time-point in nanoseconds.

    // Returns the equivalent time-duration in seconds from `d` nanoseconds.
    static constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


public:

    // Constructs a contractor for the discontinuity-graph with edge-matrix `E`
    // and `n` vertices. `P_v[j]` is to contain path-information for vertices at
    // partition `j`. Temporary files are stored at path-prefix `temp_path`.
    Discontinuity_Graph_Contractor(Edge_Matrix<k>& E, std::size_t n, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v, const std::string& temp_path);

    // Contracts the discontinuity-graph.
    void contract();
};


// Other endpoint associated to a vertex through an edge.
template <uint16_t k>
class Discontinuity_Graph_Contractor<k>::Other_End
{
private:

    Kmer<k> v_; // The endpoint vertex.
    side_t s_v_;    // Side of the endpoint to which the associated edge is incident to.
    bool is_phi_;   // Whether the endpoint is a ϕ vertex.
    weight_t w_;    // Weight of the associated edge.
    bool in_same_part_; // Whether the endpoints belong to the same partition.


public:

    // TODO: remove once we have our own hash-table.
    Other_End()
    {}

    // Constructs an endpoint with the vertex `v`, connected through its side
    // `s_v` to the corresponding edge. `is_phi` should be `true` iff `v` is the
    // ϕ vertex. The connecting edge has weight `w`, and `in_same_part` should
    // be `true` iff the endpoints of the edge belong to the same partition.
    Other_End(const Kmer<k>& v, const side_t s_v, const bool is_phi, const weight_t w, const std::size_t in_same_part):
          v_(v)
        , s_v_(s_v)
        , is_phi_(is_phi)
        , w_(w)
        , in_same_part_(in_same_part)
    {}

    // Returns the endpoint vertex.
    auto v() const { return v_; }

    // Returns the side of the endpoint to which the associated edge is incident
    // to.
    auto s_v() const { return s_v_; }

    // Returns whether the endpoint is a ϕ vertex.
    auto is_phi() const { return is_phi_; }

    // Returns the weight of the associated edge.
    auto w() const { return w_; }

    // Returns whether the endpoints belong to the same partition.
    auto in_same_part() const { return in_same_part_; }
};


template <uint16_t k>
inline void Discontinuity_Graph_Contractor<k>::form_meta_vertex(const Kmer<k> v, const std::size_t part, const side_t s_1, const weight_t w_1, const weight_t w_2)
{
    assert(part < P_v.size());

    P_v_local[parlay::worker_id()].data().emplace_back(v, v, (s_1 == side_t::back ? w_2 : w_1), side_t::back);  // The path-traversal exits `v` through its back.
}

}



#endif
