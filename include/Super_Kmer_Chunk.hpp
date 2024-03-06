
#ifndef SUPER_KMER_CHUNK_HPP
#define SUPER_KMER_CHUNK_HPP



#include "Super_Kmer_Attributes.hpp"
#include "Kmer_Utility.hpp"

#include <cstdint>
#include <cstddef>
#include <fstream>


namespace cuttlefish
{

// =============================================================================
// A chunk of super k-mers: their attributes and labels. `Colored_` denotes
// whether each super k-mer has an associated source ID.
template <bool Colored_>
class Super_Kmer_Chunk
{

private:

    typedef Super_Kmer_Attributes<Colored_> attribute_t;
    typedef attribute_t* attribute_buf_t;
    typedef uint64_t label_unit_t;
    typedef label_unit_t* label_buf_t;

    const std::size_t max_sup_kmer_len; // Maximum length of the (weak) super k-mers.
    const std::size_t sup_kmer_word_c;  // Number of 64-bit words in super k-mer encodings.

    const std::size_t cap_; // Maximum capacity of the chunk in number of super k-mers.
    std::size_t size_;  // Size of the chunk in number of super k-mers.

    attribute_buf_t const att_buf;  // Buffer of attributes of the super k-mers.
    label_buf_t const label_buf;    // Buffer of concatenated labels of the super k-mers.


public:

    // Constructs a super k-mer chunk for `k`-mers and `l`-minimizers, with
    // maximum capacity `cap`.
    Super_Kmer_Chunk(uint16_t k, uint16_t l, std::size_t cap);

    ~Super_Kmer_Chunk();

    // Returns the number of super k-mers in the chunk.
    std::size_t size() const { return size_; }

    // Clears the chunk.
    void clear() { size_ = 0; }

    // Serializes the chunk to the stream `os`.
    void serialize(std::ofstream& os) const;

    // Deserializes a chunk from the stream `is` with `sz` super k-mers.
    void deserialize(std::ifstream& is, std::size_t sz);
};


}



#endif
