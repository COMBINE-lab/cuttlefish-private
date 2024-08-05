
#ifndef COLOR_TABLE_HPP
#define COLOR_TABLE_HPP



#include "boost/unordered/concurrent_flat_map.hpp"

#include <cstdint>
#include <cstddef>
#include <utility>


namespace cuttlefish
{

enum class Color_Status;    // Extraction-status of a color-set.


// Coordinate of a color in the actual color-collection.
class Color_Coordinate
{
    typedef uint64_t pack_t;

private:

    // Flag to denote whether the corresponding color is in the process of
    // extraction or not.
    static constexpr pack_t in_process = (pack_t(1) << (sizeof(pack_t) * 8 - 1));

    static constexpr uint32_t idx_pos = 9;  // Position of the index (in worker-local bucket) of a color-set.

    pack_t bit_pack;    // Packed representation of the color-coordinate.

public:

    // Constructs an empty coordinate.
    constexpr Color_Coordinate(): bit_pack(0)
    {}

    // Constructs the coordinate `(w_id, idx)`.
    Color_Coordinate(uint16_t w_id, std::size_t idx): bit_pack(w_id | (pack_t(idx) << idx_pos))
    {}

    // Returns whether the corresponding color is in the process of extraction
    // or not.
    constexpr bool is_in_process() const { return bit_pack & in_process; }

    // Marks that the corresponding color is in the process of extraction.
    constexpr void mark_in_process() { bit_pack |= in_process; }

    // Returns the designated coordinate of a color that is being extracted.
    static constexpr Color_Coordinate in_process_coordinate() { Color_Coordinate c; c.mark_in_process(); return c; }
};


// Hashtable for color-sets. Keys are color-set hashes and values are color-set
// coordinates.
class Color_Table
{

    typedef uint64_t hash_t;
    typedef Color_Coordinate coord_t;

private:

    boost::unordered::concurrent_flat_map<hash_t, coord_t> M;


public:

    // Constructs an empty color-table.
    Color_Table();

    // Reserves enough space for at least `n` elements in the table.
    void reserve(std::size_t n) { M.reserve(n); }

    // Returns the size of the table.
    auto size() const { return M.size(); }

    // Returns `true` iff the table contains the key `h`.
    bool contains(const hash_t h) const { return M.contains(h); }

    // Tries to insert the key `h` with value `c` to the table and returns
    // `true` iff the `h` was absent in the table prior to the insertion.
    bool add(hash_t h, coord_t c) { return M.emplace(h, c); }

    // Marks that the color with hash `h` is in the process of extraction, if a
    // corresponding coordinate for `h` does not already exist in the table. If
    // it does, it is put in `c`. Returns the extraction-status of the color
    // prior to this invocation.
    Color_Status mark_in_process(hash_t h, Color_Coordinate& c);
};


// Extraction-status of a color-set.
enum class Color_Status
{
    undiscovered,   // has not been seen yet
    in_process,     // is in the process of extraction
    discovered,     // completely extracted
};


inline Color_Status Color_Table::mark_in_process(const hash_t h, Color_Coordinate& c)
{
    const auto r = M.insert_or_visit({h, Color_Coordinate::in_process_coordinate()},
                    [&](const auto& p)
                    {
                        c = p.second;
                    });

    typedef Color_Status status_t;
    return r ? status_t::undiscovered : (c.is_in_process() ? status_t::in_process : status_t::discovered);
}

}



#endif