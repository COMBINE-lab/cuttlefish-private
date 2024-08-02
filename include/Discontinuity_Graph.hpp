
#ifndef DISCONTINUITY_GRAPH_HPP
#define DISCONTINUITY_GRAPH_HPP



#include "Kmer.hpp"
#include "Edge_Matrix.hpp"
#include "Unitig_File.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <cstddef>
#include <utility>
#include <atomic>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// A representation of a discontinuity graph of `k`-mers.
template <uint16_t k>
class Discontinuity_Graph
{
private:

    // k-mer (super-)label of the ϕ-vertex in the discontinuity graph.  // TODO: revisit; almost sure we don't need this.
    static constexpr const char phi_label[] =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    static const Kmer<k> phi_;  // ϕ k-mer connected to each chain-end in the discontinuity graph.

    Edge_Matrix<k> E_;  // Edge-matrix of the discontinuity graph.

    Unitig_Write_Distributor lmtigs;    // Distribution-manager for the writes of locally maximal unitigs' labels.

    std::atomic_uint64_t phantom_edge_count_;   // Number of potential phantom edges identified.


public:

    // Constructs a discontinuity graph object that operates with `part_count`
    // vertex-partitions, and the locally-maximal unitigs corresponding to its
    // edges are stored in `lmtig_bucket_count` buckets. `logistics` is the data
    // logistics manager for the algorithm execution.
    Discontinuity_Graph(std::size_t part_count, std::size_t lmtig_bucket_count, const Data_Logistics& logistics);

    // Returns the ϕ k-mer connected to each chain-end in the discontinuity
    // graph.
    static const Kmer<k>& phi() { return phi_; }

    // Returns the edge-matrix of the graph.
    const Edge_Matrix<k>& E() const { return E_; }

    // Returns the edge-matrix of the graph.
    Edge_Matrix<k>& E() { return E_; }

    // Returns the number of potential phantom edges identified.
    uint64_t phantom_edge_upper_bound() const;

    // Adds the edge `({(u, s_u), (v, s_v)}, 1)` to the graph. `u_is_phi` and
    // `v_is_phi` denote whether `u` and `v` are `ϕ` respectively. The locally-
    // maximal unitig corresponding to the edge is `mtig`. The edge should be
    // an original edge of the graph. Returns `(b, b_idx)`, where the deposited
    // lm-tig is put into bucket `b` at index `b_idx`.
    std::pair<std::size_t, std::size_t> add_edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, bool u_is_phi, bool v_is_phi, const Maximal_Unitig_Scratch<k>& mtig);

    // Adds the edge `({ϕ, (v, s_v)}, 1)` to the graph. The edge should be an
    // original edge of the graph.
    void add_edge(const Kmer<k>& v, side_t s_v);

    // Adds the edge `({(u, s_u), (v, s_v)}, w)` to the graph. `u_is_phi` and
    // `v_is_phi` denote whether `u` and `v` are `ϕ` respectively. The edge
    // should be a contracted edge and not an original one.
    void add_edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, weight_t w, bool u_is_phi, bool v_is_phi);

    // Increments the potential phantom edge count.
    void inc_potential_phantom_edge() { phantom_edge_count_++; }

    // Closes the lm-tig writer streams.
    void close_lmtig_stream();

    // Returns a tight upper bound of the maximum number of vertices in a
    // partition.
    std::size_t vertex_part_size_upper_bound() const;
};


template <uint16_t k> const Kmer<k> Discontinuity_Graph<k>::phi_(phi_label);


template <uint16_t k>
inline std::pair<std::size_t, std::size_t> Discontinuity_Graph<k>::add_edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const bool u_is_phi, const bool v_is_phi, const Maximal_Unitig_Scratch<k>& mtig)
{
    const auto w_id = parlay::worker_id();
    const auto b = lmtigs.file_idx(w_id);
    const std::size_t b_idx = lmtigs.unitig_count(b);
    lmtigs.add(w_id, mtig);

    E_.add(u, s_u, v, s_v, 1, b, b_idx, u_is_phi, v_is_phi);

    return {b, b_idx};
}


template <uint16_t k>
inline void Discontinuity_Graph<k>::add_edge(const Kmer<k>& v, const side_t s_v)
{
    const auto w_id = parlay::worker_id();
    const auto b = lmtigs.file_idx(w_id);
    std::size_t b_idx = lmtigs.unitig_count(b);
    lmtigs.add(w_id, s_v == side_t::front ? v : v.reverse_complement());

    E_.add(phi_, side_t::back, v, s_v, 1, b, b_idx, true, false);
}


template <uint16_t k>
inline void Discontinuity_Graph<k>::add_edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const weight_t w, const bool u_is_phi, const bool v_is_phi)
{
    // Edge-partition 0 associates to edges that do not have any corresponding lm-tig (i.e. has weight > 1).
    E_.add(u, s_u, v, s_v, w, 0, 0, u_is_phi, v_is_phi);
}

}



#endif
