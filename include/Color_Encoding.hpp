
#ifndef COLOR_ENCODING_HPP
#define COLOR_ENCODING_HPP



#include <cstdint>
#include <cassert>


namespace cuttlefish
{

// Coordinate of a color in the actual color-collection.
class Color_Coordinate
{
    typedef uint64_t pack_t;

private:

    // Flag to denote whether the corresponding color is in the process of
    // extraction or not.
    static constexpr pack_t in_process = (pack_t(1) << (sizeof(pack_t) * 8 - 1));

    static constexpr uint32_t idx_pos = 8;  // Position of the index (in worker-local bucket) of a color-set.

    static constexpr pack_t w_limit = (pack_t(1) << idx_pos);   // Maximum worker count: 2^8.
    static constexpr pack_t idx_limit  = (pack_t(1) << 32); // Maximum worker-local bucket size: 2^32.

    pack_t bit_pack;    // Packed representation of the color-coordinate.

public:

    // Constructs an empty coordinate.
    Color_Coordinate(): bit_pack(0)
    {}

    // Constructs an empty coordinate marked as being processed by worker-ID
    // `w`.
    explicit Color_Coordinate(const pack_t w_id): bit_pack(w_id | in_process) { assert(w_id < w_limit); }

    // Constructs the coordinate `(w_id, idx)`.
    Color_Coordinate(const pack_t w_id, const pack_t idx): bit_pack(w_id | (pack_t(idx) << idx_pos)) { assert(w_id < w_limit); assert(idx < idx_limit); }

    // Returns whether the corresponding color is in the process of extraction
    // or not.
    bool is_in_process() const { return bit_pack & in_process; }

    // Returns the worker-ID that marked this coordinate as processing.
    pack_t processing_worker() const { assert(is_in_process()); return bit_pack & (~in_process); }

    // Returns 40-bit packing of the color-coordinate.
    auto as_u40() const { assert(bit_pack < (pack_t(1) << 40)); return bit_pack; }
};

}



#endif
