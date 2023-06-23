
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

    const path_id_t p_; // The path-ID of the vertex.
    const weight_t r_;  // The rank of the vertex in the path.
    const side_t o_;    // The orientation of the vertex in its specified rank—the path traversal exits `v` through the side `o_v`.


public:

    // Constructs a path-info object for a vertex such that its path-ID is `p`
    // and rank in the path is `r` when the path is traversed in the
    // orientation such that the traversal exits the vertex through its side
    // `o`.
    Vertex_Path_Info(const path_id_t p, const weight_t r, const side_t o):
          p_(p)
        , r_(r)
        , o_(o)
    {}

    // Returns the path-ID of the vertex.
    const path_id_t p() const { return p_; }

    // Returns the rank of the vertex in the path.
    weight_t r() const { return r_; }

    // Returns the orientation `o_v` of the vertex in its specified rank—the
    // path traversal exits the vertex through the side `o_v`.
    side_t o() const { return o_; }
};

}



#endif
