
#ifndef SUPER_KMER_CHUNK_HPP
#define SUPER_KMER_CHUNK_HPP



#include "Super_Kmer_Attributes.hpp"
#include "Kmer_Utility.hpp"

#include <cstdint>
#include <vector>
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
    typedef std::vector<attribute_t> attribute_buf_t;
    typedef uint64_t label_unit_t;
    typedef std::vector<label_unit_t> label_buf_t;

    const std::size_t max_sup_kmer_len; // Maximum length of the (weak) super k-mers.
    const std::size_t sup_kmer_word_c;  // Number of 64-bit words in super k-mer encodings.

    attribute_buf_t att_buf;    // Buffer of attributes of the super k-mers.
    label_buf_t label_buf;  // Buffer of concatenated labels of the super k-mers.


public:

    // Constructs a super k-mer chunk for `k`-mers and `l`-minimizers.
    Super_Kmer_Chunk(uint16_t k, uint16_t l);

    // Returns the current size of the chunk in bytes.
    std::size_t bytes() const;

    // Returns the number of super k-mers in the chunk.
    std::size_t size() const { return att_buf.size(); }

    // Clears the chunk.
    void clear();

    // Serializes the chunk to the stream `os`.
    void serialize(std::ofstream& os) const;

    // Deserializes a chunk from the stream `is` with `sz` super k-mers.
    void deserialize(std::ifstream& is, std::size_t sz);
};


template <bool Colored_>
inline std::size_t Super_Kmer_Chunk<Colored_>::bytes() const
{
    return att_buf.size() * sizeof(attribute_t) + label_buf.size() * sizeof(label_unit_t);
}

}



#endif
