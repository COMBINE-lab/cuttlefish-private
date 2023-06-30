
#ifndef PATH_INFO_HPP
#define PATH_INFO_HPP



#include "Kmer.hpp"
#include "globals.hpp"

#include <cstdint>


namespace cuttlefish
{

// =============================================================================
// Path-information of an object in a discontinuity graph: its path-ID, rank in
// a fixed traversal of the path, and orientation in that traversal.
template <uint16_t k>
class Path_Info
{
public:

    typedef Kmer<k> path_id_t;


private:

    path_id_t p_;   // The path-ID.
    weight_t r_;    // The rank.
    side_t o_;  // The orientation of the vertex in its specified rank—the path traversal exits `v` through the side `o`.


public:

    Path_Info()
    {}


    // Constructs a path-info object for an object such that its path-ID is `p`
    // and rank in the path is `r` when the path is traversed in the
    // orientation such that the traversal exits the object through its side
    // `o`.
    Path_Info(const path_id_t p, const weight_t r, const side_t o):
          p_(p)
        , r_(r)
        , o_(o)
    {}


    // Returns the path-ID.
    const path_id_t p() const { return p_; }

    // Returns the rank.
    weight_t r() const { return r_; }

    // Returns the orientation `o` of the object in its specified rank—the
    // path traversal exits the object through the side `o`.
    side_t o() const { return o_; }

    // Returns `true` iff this information is the same as in `rhs`.
    bool operator==(const Path_Info<k>& rhs) const { return p_ == rhs.p_ && r_ == rhs.r_ && o_ == rhs.o_; }
};


// A vertex and its path-information.
template <uint16_t k>
class Vertex_Path_Info_Pair
{
private:

    Kmer<k> v_; // The vertex.
    Path_Info<k> path_info_;    // Path-information of the vertex.


public:

    typedef typename Path_Info<k>::path_id_t path_id_t;

    Vertex_Path_Info_Pair()
    {}


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

    // Returns `true` iff this information is the same as in `rhs`.
    bool operator==(const Vertex_Path_Info_Pair<k>& rhs) const { return v_ == rhs.v_ && path_info_ == rhs.path_info_; }
};

}



#endif
