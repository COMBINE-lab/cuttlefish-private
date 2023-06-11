
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

    Kmer<k> u_; // An endpoint of the edge.
    Kmer<k> v_; // An endpoint of the edge.
    side_t s_u; // Side of the vertex `u_` to which the edge is incident to.
    side_t s_v; // Side of the vertex `v_` to which the edge is incident to.
    bool u_is_phi_; // Whether `u_` is the ϕ vertex.
    bool v_is_phi_; // Whether `v_` is the ϕ vertex.
    uint16_t bucket_id; // ID of the bucket of this edge.
    uint16_t weight;    // Weight of the edge.


public:

    // Construct an edge `({(u, s_u), (v, s_v)}, w, b)`. `u_is_phi` and
    // `v_is_phi` denote whether `u` and `v` are ϕ, respectively.
    Discontinuity_Edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, uint32_t w, uint16_t b, bool u_is_phi = false, bool v_is_phi = false);

    // Returns the `u` endpoint of the edge.
    const Kmer<k>& u() const { return u_; }

    // Returns the `v` endpoint of the edge.
    const Kmer<k>& v() const { return v_; }

    // Returns whether `u` is the ϕ vertex.
    bool u_is_phi() const { return u_is_phi_; }

    // Returns whether `v` is the ϕ vertex.
    bool v_is_phi() const { return v_is_phi_; }

    // Inverts the `u` and the `v` endpoints of the edge.
    void invert();
};


template <uint16_t k>
inline Discontinuity_Edge<k>::Discontinuity_Edge(const Kmer<k>& u, const cuttlefish::side_t s_u, const Kmer<k>& v, const cuttlefish::side_t s_v, const uint32_t w, const uint16_t b, bool u_is_phi, bool v_is_phi):
      u_(u)
    , v_(v)
    , s_u(s_u)
    , s_v(s_v)
    , u_is_phi_(u_is_phi)
    , v_is_phi_(v_is_phi)
    , bucket_id(b)
    , weight(w)
{}


template <uint16_t k>
inline void Discontinuity_Edge<k>::invert()
{
    std::swap(u_, v_);
    std::swap(s_u, s_v);
    std::swap(u_is_phi_, v_is_phi_);
}

}



#endif
