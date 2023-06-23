
#ifndef VERTEX_PATH_INFO_HPP
#define VERTEX_PATH_INFO_HPP



#include "Kmer.hpp"
#include "globals.hpp"

#include <cstdint>


namespace cuttlefish
{

// =============================================================================
// Path-information, i.e. path-ID, rank within the path, and orientation of a
// vertex in a discontinuity-graph.
template <uint16_t k>
class Vertex_Path_Info
{
public:

    typedef Kmer<k> path_id_t;


private:

    const path_id_t p_v_;   // The path-ID of the vertex.
    const weight_t r_v_;    // The rank of the vertex in the path.
    const side_t o_v_;  // The orientation of the vertex in its specified rank—the path traversal exits `v` through the side `o_v`.


public:

    // Constructs a path-info object for a k-mer such that its path-ID is `p_v`
    // and rank in the path is `r_v` when the path is traversed in the
    // orientation such that the traversal exits `v` through its side `o_v`.
    Vertex_Path_Info(const path_id_t p_v, const weight_t r_v, const side_t o_v):
          p_v_(p_v)
        , r_v_(r_v)
        , o_v_(o_v)
    {}

    // Returns the path-ID of the vertex.
    const path_id_t p_v() const { return p_v_; }

    // Returns the rank of the vertex in the path.
    weight_t r_v() const { return r_v_; }

    // Returns the orientation `o_v` of the vertex in its specified rank—the
    // path traversal exits the vertex through the side `o_v`.
    side_t o_v() const { return o_v_; }
};

}



#endif
