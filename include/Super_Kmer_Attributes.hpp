
#ifndef SUPER_KMER_ATTRIBUTES_HPP
#define SUPER_KMER_ATTRIBUTES_HPP



#include <cstdint>
#include <cstddef>
#include <cassert>


namespace cuttlefish
{


// =============================================================================
// Collection of attributes of a super k-mer. `Colored_` denotes whether the
// super k-mer has an associated source ID.
template <bool Colored_>
class Super_Kmer_Attributes
{};


template <>
class Super_Kmer_Attributes<false>
{

private:

    uint16_t bit_pack;  // Packed attribute collection of a super k-mer.

    static constexpr uint32_t len_pos = 0;  // Bit-index of length in the pack.
    static constexpr uint32_t l_disc_pos = 8;   // Bit-index of the left discontinuity marker.
    static constexpr uint32_t r_disc_pos  = 9;  // Bit-index of the right discontinuity marker.
    static constexpr uint16_t len_mask = 0b1111'1111 << len_pos;
    static constexpr uint16_t l_disc_mask = (0b1 << l_disc_pos);
    static constexpr uint16_t r_disc_mask = (0b1 << r_disc_pos);


public:

    Super_Kmer_Attributes() {}

    // Constructs an attribute object with length (in bases) `len` and left /
    // right discontinuity markers `l_disc` and `r_disc`.
    Super_Kmer_Attributes(std::size_t len, bool l_disc, bool r_disc);

    // Returns the length of the super k-mer (in bases).
    uint8_t len() const { return (bit_pack & len_mask) >> len_pos; }

    // Returns whether the super k-mer is discontinuous on the left.
    bool left_discontinuous() const { return bit_pack & l_disc_mask; }

    // Returns whether the super k-mer is discontinuous on the right.
    bool right_discontinuous() const { return bit_pack & r_disc_mask; }
};


inline Super_Kmer_Attributes<false>::Super_Kmer_Attributes(const std::size_t len, const bool l_disc, const bool r_disc):
      bit_pack((len << len_pos) | (static_cast<uint16_t>(l_disc) << l_disc_pos) | (static_cast<uint16_t>(r_disc) << r_disc_pos))
{
    assert(len <= (len_mask >> len_pos));
}

}



#endif
