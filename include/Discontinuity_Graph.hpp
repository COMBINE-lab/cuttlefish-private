
#ifndef DISCONTINUITY_GRAPH_HPP
#define DISCONTINUITY_GRAPH_HPP



#include "Kmer.hpp"
#include "Edge_Matrix.hpp"
#include "Unitig_File.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <cstddef>


namespace cuttlefish
{

// =============================================================================
// A representation of a discontinuity graph of `k`-mers.
template <uint16_t k>
class Discontinuity_Graph
{
private:

    Edge_Matrix<k> E_;  // Edge-matrix of the discontinuity graph.

    Unitig_Write_Distributor lmtigs;    // Distribution-manager for the writes of locally maximal unitigs' labels.


public:

    // Constructs a discontinuity graph object that operates with `part_count`
    // vertex-partitions, and the locally-maximal unitigs corresponding to its
    // edges are stored in `lmtig_bucket_count` buckets. Temporary working
    // files are stored at path-prefix `work_path`.
    Discontinuity_Graph(std::size_t part_count, std::size_t lmtig_bucket_count, const std::string& work_path);

    // Adds the edge `({(u, s_u), (v, s_v)}, 1)` to the graph. `u_is_phi` and
    // `v_is_phi` denote whether `u` and `v` are `ϕ` respectively. The locally-
    // maximal unitig corresponding to the edge is `mtig`. The edge should be an
    // original edge of the graph.
    void add_edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, bool u_is_phi, bool v_is_phi, const Maximal_Unitig_Scratch<k>& mtig);

    // Adds the edge `({(u, s_u), (v, s_v)}, w)` to the graph. `u_is_phi` and
    // `v_is_phi` denote whether `u` and `v` are `ϕ` respectively. The edge
    // should be a contracted edge and not an original one.
    void add_edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, weight_t w, bool u_is_phi, bool v_is_phi);

    // Closes the lm-tig writer streams.
    void close_lmtig_stream();

    // Returns the edge-matrix of the graph.
    const Edge_Matrix<k>& E() const { return E_; }
};


template <uint16_t k>
inline void Discontinuity_Graph<k>::add_edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const bool u_is_phi, const bool v_is_phi, const Maximal_Unitig_Scratch<k>& mtig)
{
    const auto w_id = parlay::worker_id();
    const auto b = lmtigs.file_idx(w_id);
    std::size_t b_idx = lmtigs.unitig_count(b);
    lmtigs.add(w_id, mtig);

    E_.add(u, s_u, v, s_v, 1, b, b_idx, u_is_phi, v_is_phi);
}


template <uint16_t k>
inline void Discontinuity_Graph<k>::add_edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const weight_t w, const bool u_is_phi, const bool v_is_phi)
{
    // Edge-partition 0 associates to edges that do not have any corresponding lm-tig (i.e. has weight > 1).
    E_.add(u, s_u, v, s_v, w, 0, 0, u_is_phi, v_is_phi);
}

}



#endif
