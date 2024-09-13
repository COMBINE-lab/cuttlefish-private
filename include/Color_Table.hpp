
#ifndef COLOR_TABLE_HPP
#define COLOR_TABLE_HPP



#include "Color_Encoding.hpp"
#include "boost/unordered/concurrent_flat_map.hpp"

#include <cstdint>
#include <cstddef>
#include <cassert>


namespace cuttlefish
{

enum class Color_Status;    // Extraction-status of a color-set.


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

    // Marks that the color with hash `h` is in the process of extraction by
    // the `w`'th worker, if a corresponding entry for `h` does not already
    // exist in the table. If it does, it is put in `c`. Returns the
    // extraction-status of the color prior to this invocation.
    Color_Status mark_in_process(hash_t h, uint64_t w, Color_Coordinate& c);

    // Updates the key `h` with value `c` if `h` is marked as in process of
    // extraction. Returns `true` iff the update is successful.
    bool update_if_in_process(hash_t h, Color_Coordinate c);

    // Assigns `c` to the value of the key `h`.
    void assign(hash_t h, coord_t c);

    // Returns the value associated to the key `h`.
    Color_Coordinate get(hash_t h) const;
};


// Extraction-status of a color-set.
enum class Color_Status
{
    undiscovered,   // has not been seen yet
    in_process,     // is in the process of extraction
    discovered,     // completely extracted
};


inline Color_Status Color_Table::mark_in_process(const hash_t h, const uint64_t w, Color_Coordinate& c)
{
    const auto r =  M.emplace_or_cvisit(h, Color_Coordinate(w), [&](const auto& p)
                    {
                        c = p.second;
                    });

    typedef Color_Status status_t;
    return r ? status_t::undiscovered : (c.is_in_process() ? status_t::in_process : status_t::discovered);
}


inline bool Color_Table::update_if_in_process(const hash_t h, const Color_Coordinate c)
{
    assert(M.contains(h));

    bool was_in_process = false;
    const auto r =  M.visit(h, [&](auto& p)
                    {
                        auto& color = p.second;
                        was_in_process = color.is_in_process();
                        color = (was_in_process ? c : color);
                    });
    assert(r); (void)r;

    return was_in_process;
}


inline Color_Coordinate Color_Table::get(const hash_t h) const
{
    assert(M.contains(h));

    Color_Coordinate c;
    const auto r =  M.cvisit(h, [&](const auto& p)
                    {
                        c = p.second;
                    });
    assert(r); (void)r;

    return c;
}

}



#endif