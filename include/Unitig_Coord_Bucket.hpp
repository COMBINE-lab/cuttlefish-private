
#ifndef UNITIG_COORD_BUCKET_HPP
#define UNITIG_COORD_BUCKET_HPP



#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "globals.hpp"

#include <cstdint>
#include <string>


namespace cuttlefish
{

// =============================================================================
// A bucket storing full coordinates for unitigs: for a specific unitig, it's
// containing maximal unitig's unique ID, its rank in the maximal unitig in a
// fixed traversal of the path, its orientation in that traversal, and its
// literal label.
template <uint16_t k>
class Unitig_Coord_Bucket
{
    typedef max_unitig_id_t<k> path_id_t;   // Type of a maximal unitig path's ID.

    // Coordinate information of a unitig.
    struct Unitig_Coord
    {
        Path_Info<k> path_info; // Coordinate of the unitig in the de Bruijn graph.
        uint32_t label_idx; // Index of the label of the unitig into the dump-string of the bucket.
        uni_len_t label_len;    // Length of the label of the unitig.

        Unitig_Coord(const Path_Info<k>& path_info, const uint32_t label_idx, const uni_len_t label_len):
              path_info(path_info)
            , label_idx(label_idx)
            , label_len(label_len)
        {}
    };

private:

    const std::string path_pref;    // Path-prefix to the file(s) storing the bucket.

    Ext_Mem_Bucket<Unitig_Coord> coord_bucket;  // External-memory bucket of the unitig-coordinates.
    Ext_Mem_Bucket<char> label_bucket;  // External-memory bucket of the unitig-labels.


public:

    // Constructs a unitig-coordinate bucket at path-prefix `file_path`.
    Unitig_Coord_Bucket(const std::string& path_pref);

    // Adds a unitig to the bucket with its path-information in the de Bruijn
    // graph `path_info`, label `label`, and length `len`.
    void add(const Path_Info<k>& path_info, const char* label, uni_len_t len);
};


template <uint16_t k>
inline void Unitig_Coord_Bucket<k>::add(const Path_Info<k>& path_info, const char* const label, const uni_len_t len)
{
    coord_bucket.emplace(path_info, label_bucket.size(), len);
    label_bucket.add(label, len);
}

}



#endif
