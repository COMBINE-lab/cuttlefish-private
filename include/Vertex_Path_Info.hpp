
#ifndef VERTEX_PATH_INFO_HPP
#define VERTEX_PATH_INFO_HPP



#include "Kmer.hpp"
#include "globals.hpp"

#include <cstdint>


namespace cuttlefish
{

// =============================================================================
// Path-information of a vertex in a discontinuity graph: its path-ID, rank in a
// fixed traversal of the path, and orientation in that traversal.
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


// A vertex and its path-information.
template <uint16_t k>
class Vertex_Path_Info_Pair
{
private:

    Kmer<k> v_; // The vertex.
    Vertex_Path_Info<k> path_info_; // Path-information of the vertex.


public:

    typedef typename Vertex_Path_Info<k>::path_id_t path_id_t;

    // For a vertex `v`, constructs a pairing of it with its path-info specified
    // with its path-ID `p` and rank in the path `r` when the path is traversed
    // in the orientation such that the traversal exits `v` through side `o`.
    Vertex_Path_Info_Pair(const Kmer<k> v, const path_id_t p, const weight_t r, const side_t o):
          v_(v)
        , path_info_(p, r, o)
    {}


    // Returns the vertex.
    const auto v() const { return v_; }

    // Returns the path-info of the vertex.
    const auto path_info() const { return path_info_; }
};

}



#endif
