
#ifndef COLOR_TABLE_HPP
#define COLOR_TABLE_HPP



#include "boost/unordered/concurrent_flat_map.hpp"

#include <cstdint>
#include <cstddef>


namespace cuttlefish
{

// Hashtable for color-sets. Keys are color-set hashes and values are color-set
// coordinates.
class Color_Table
{

    typedef uint64_t hash_t;
    typedef uint32_t coord_t;

private:

    boost::unordered::concurrent_flat_map<hash_t, coord_t> M;


public:

    // Constructs an empty color-table.
    Color_Table();

    // Reserves enough space for at least `n` elements in the table.
    void reserve(std::size_t n) { M.reserve(n); }

    // Returns `true` iff the table contains the key `h`.
    bool contains(const hash_t h) const { return M.contains(h); }

    // Tries to insert the key `h` with value `c` to the table and returns
    // `true` iff the `h` was absent in the table prior to the insertion.
    bool add(hash_t h, coord_t c) { return M.emplace(h, c); }

    // Returns the size of the table.
    auto size() const { return M.size(); }
};

}



#endif