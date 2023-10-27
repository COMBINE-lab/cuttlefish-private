
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
// Coordinate information of a unitig, both in the de Bruijn graph and in the
// dump-string in a bucket.
template <uint16_t k>
class Unitig_Coord
{
    typedef uint32_t label_idx_t;

private:

    Path_Info<k> path_info; // Coordinate of the unitig in the de Bruijn graph.
    label_idx_t label_idx_; // Index of the label of the unitig into the dump-string of the bucket.
    uni_len_t label_len_;   // Length of the label of the unitig.


public:

    Unitig_Coord(const Path_Info<k>& path_info, const uint32_t label_idx, const uni_len_t label_len):
          path_info(path_info)
        , label_idx_(label_idx)
        , label_len_(label_len)
    {}

    // Returns the path-ID of the unitig.
    auto p() const { return path_info.p(); }

    // Returns the rank of the unitig in a specific traversal orientation.
    auto r() const { return path_info.r(); }

    // Returns the orientation `o` of the unitig in its specified rank—the path-
    // traversal exits the object through the side `o`.
    auto o() const { return path_info.o(); }

    // Returns the index of the label of the unitig into the dump-string of the
    // bucket.
    auto label_idx() const { return label_idx_; }

    // Returns the length of the label of the unitig.
    auto label_len() const { return label_len_; }

    // Returns `true` iff this coordinate's path-info is lexicographically
    // smaller than `rhs`'s path-info.
    bool operator<(const Unitig_Coord& rhs) const { return path_info < rhs.path_info; }
};


// =============================================================================
// A bucket storing full coordinates for unitigs: for a specific unitig, it's
// containing maximal unitig's unique ID, its rank in the maximal unitig in a
// fixed traversal of the path, its orientation in that traversal, and
// additionally its literal label.
template <uint16_t k>
class Unitig_Coord_Bucket
{
    typedef max_unitig_id_t<k> path_id_t;   // Type of a maximal unitig path's ID.

private:

    const std::string path_pref;    // Path-prefix to the file(s) storing the bucket.

    Ext_Mem_Bucket<Unitig_Coord<k>> coord_bucket;   // External-memory bucket of the unitig-coordinates.
    Ext_Mem_Bucket<char> label_bucket;  // External-memory bucket of the unitig-labels.

    std::size_t size_;  // Number of unitigs stored in the bucket.

    std::size_t label_len_; // Total length of the labels of the stored unitigs.


public:

    // Constructs a unitig-coordinate bucket at path-prefix `file_path`.
    Unitig_Coord_Bucket(const std::string& path_pref);

    // Returns the number of unitigs stored in the bucket.
    std::size_t size() const { return size_; }

    // Returns the total length of the labels of the stored unitigs.
    std::size_t label_len() const { return label_len_; }

    // Adds a unitig to the bucket with its path-information in the de Bruijn
    // graph `path_info`, label `label`, and length `len`.
    void add(const Path_Info<k>& path_info, const char* label, uni_len_t len);

    // Loads all the unitig-coordinates in the bucket to `buf`, and returns this
    // size.
    std::size_t load_coords(Unitig_Coord<k>* buf) const;

    // Loads the concatenated label string of the entire bucket into `buf`, and
    // returns its length.
    std::size_t load_labels(char* buf) const;
};


template <uint16_t k>
inline void Unitig_Coord_Bucket<k>::add(const Path_Info<k>& path_info, const char* const label, const uni_len_t len)
{
    coord_bucket.emplace(path_info, label_bucket.size(), len);
    label_bucket.add(label, len);

    size_++;
    label_len_ += len;
}

}



#endif
