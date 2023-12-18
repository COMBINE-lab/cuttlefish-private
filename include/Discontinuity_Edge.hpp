
#ifndef DISCONTINUITY_EDGE_HPP
#define DISCONTINUITY_EDGE_HPP



#include "Kmer.hpp"
#include "globals.hpp"

#include <cstdint>
#include <algorithm>


namespace cuttlefish
{

// =============================================================================
// An edge `e_{u, v} = ({(u, s_u), (v, s_v)}, w, e_b)` in a discontinuity-graph
// of `k`-mers.
template <uint16_t k>
class Discontinuity_Edge
{
private:

    // TODO: try packing booleans and sides.

    Kmer<k> u_; // An endpoint of the edge.
    Kmer<k> v_; // An endpoint of the edge.
    side_t s_u_;    // Side of the vertex `u_` to which the edge is incident to.
    side_t s_v_;    // Side of the vertex `v_` to which the edge is incident to.
    bool u_is_phi_; // Whether `u_` is the ϕ vertex.
    bool v_is_phi_; // Whether `v_` is the ϕ vertex.
    weight_t weight;    // Weight of the edge.
    uint16_t bucket_id; // ID of the bucket of the unitig corresponding to the edge.
    uni_idx_t b_idx_;   // Index of the corresponding unitig within its bucket.
    side_t o_;  // Exit-orientation of the corresponding literal unitig wrt the `(u, v)` orientation of the edge.

    // k-mer (super-)label of the ϕ-vertex in the discontinuity graph.  // TODO: revisit; almost sure we don't need this.
    static constexpr const char phi_label[] =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    static const Kmer<k> phi_;  // ϕ k-mer connected to each chain-end in the discontinuity graph.


public:

    // Constructs a placeholder edge.
    Discontinuity_Edge(){}

    // Constructs an edge `{(u, s_u), (v, s_v)}` between the vertices `u` and
    // `v` that connects there sides `s_u` and `s_v` resp. It has weight `w`.
    // The locally-maximal unitig corresponding to this edge is stored in the
    // `b`'th bucket, at index `b_idx`. `u_is_phi` and `v_is_phi` denote whether
    // `u` and `v` are ϕ, respectively. `o` is the exit-orientation of the
    // corresponding literal unitig (if any) wrt to the `(u, v)` orientation.
    Discontinuity_Edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, weight_t w, uint16_t b, std::size_t b_idx, bool u_is_phi, bool v_is_phi, side_t o);

    // Returns the `u` endpoint of the edge.
    const Kmer<k>& u() const { return u_; }

    // Returns the `v` endpoint of the edge.
    const Kmer<k>& v() const { return v_; }

    // Returns the side of the `u` endpoint to which the edge is incident to.
    side_t s_u() const { return s_u_; }

    // Returns the side of the `v` endpoint to which the edge is incident to.
    side_t s_v() const { return s_v_; }

    // Returns the `u` endpoint of the edge.
    const Kmer<k>& x() const { return u(); }

    // Returns the `v` endpoint of the edge.
    const Kmer<k>& y() const { return v(); }

    // Returns the side of the `u` endpoint to which the edge is incident to.
    side_t s_x() const { return s_u(); }

    // Returns the side of the `v` endpoint to which the edge is incident to.
    side_t s_y() const { return s_v(); }

    // Returns the weight of the edge.
    weight_t w() const { return weight; }

    // Returns the ID of the bucket of this edge.
    uint16_t b() const { return bucket_id; }

    // Returns the index of the corresponding unitig within its bucket.
    std::size_t b_idx() const { return b_idx_; }

    // Returns whether `u` is the ϕ vertex.
    bool u_is_phi() const { return u_is_phi_; }

    // Returns whether `v` is the ϕ vertex.
    bool v_is_phi() const { return v_is_phi_; }

    // Returns whether `u` is the ϕ vertex.
    bool x_is_phi() const { return u_is_phi(); }

    // Returns whether `v` is the ϕ vertex.
    bool y_is_phi() const { return v_is_phi(); }

    // Returns exit-orientation of the corresponding literal unitig wrt the
    // `(u, v)` orientation of the edge.
    side_t o() const { return o_; }

    // Returns the ϕ k-mer connected to each chain-end in the discontinuity
    // graph.
    static const Kmer<k>& phi() { return phi_; }
};


template <uint16_t k>
inline Discontinuity_Edge<k>::Discontinuity_Edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const weight_t w, const uint16_t b, const std::size_t b_idx, const bool u_is_phi, const bool v_is_phi, const side_t o):
      u_(u)
    , v_(v)
    , s_u_(s_u)
    , s_v_(s_v)
    , u_is_phi_(u_is_phi)
    , v_is_phi_(v_is_phi)
    , weight(w)
    , bucket_id(b)
    , b_idx_(b_idx)
    , o_(o)
{}

}



#endif
