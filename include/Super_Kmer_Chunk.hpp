
#ifndef SUPER_KMER_CHUNK_HPP
#define SUPER_KMER_CHUNK_HPP



#include "Super_Kmer_Attributes.hpp"
#include "Kmer_Utility.hpp"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <vector>
#include <fstream>


namespace cuttlefish
{

// =============================================================================
// A chunk of weak super k-mers: their attributes and labels. `Colored_` denotes
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

    const std::size_t cap_; // Maximum capacity of the chunk in number of super k-mers.
    std::size_t size_;  // Size of the chunk in number of super k-mers.

    attribute_buf_t att_buf;    // Buffer of attributes of the super k-mers.
    label_buf_t label_buf;  // Buffer of concatenated labels of the super k-mers.


public:

    // Constructs a super k-mer chunk for `k`-mers and `l`-minimizers, with
    // maximum capacity `cap` in number of super k-mers.
    Super_Kmer_Chunk(uint16_t k, uint16_t l, std::size_t cap);

    Super_Kmer_Chunk(const Super_Kmer_Chunk&) = delete;

    Super_Kmer_Chunk(Super_Kmer_Chunk&&) = default;

    // Returns the number of super k-mers in the chunk.
    auto size() const { return size_; }

    // Returns the maximum capacity of the chunk in number of super k-mers.
    auto capacity() const {return cap_;  }

    // Returns the free capacity of the chunk in number of super k-mers.
    auto free_capacity() const { return capacity() - size(); }

    // Returns whether the chunk is full or not.
    auto full() const { return size() == capacity(); }

    // Returns the number of units, i.e. 64-bit words, in the label buffer.
    auto label_units() const { return size() * sup_kmer_word_c; }

    // Clears the chunk.
    void clear() { size_ = 0; }

    // Returns the size of a super k-mer record in bytes, that is over `k`-mers
    // and `l`-minimizers.
    static std::size_t record_size(uint16_t k, uint16_t l);

    // Serializes the chunk to the stream `os`.
    void serialize(std::ofstream& os) const;

    // Deserializes a chunk from the stream `is` with `sz` super k-mers.
    void deserialize(std::ifstream& is, std::size_t sz);

    // Adds a super k-mer to the chunk with label `seq` and length `len`. The
    // markers `l_disc` and `r_disc` denote whether the left and the right ends
    // of the (weak) super k-mer are discontinuous or not.
    void add(const char* seq, std::size_t len, bool l_disc, bool r_disc);

    // Appends the chunk `c`'s contents in the indices `[l, r)` to this chunk.
    void append(const Super_Kmer_Chunk& c, std::size_t l, std::size_t r);

    // Appends the chunk `c` to the end of this chunk.
    void append(const Super_Kmer_Chunk& c);
};


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::add(const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc)
{
    assert(size() < cap_);

    att_buf[size()] = Super_Kmer_Attributes<Colored_>(len, l_disc, r_disc);

    const auto label_off = label_units();   // Offset into the packed-encoding concatenation where to put the label.
    int64_t word_idx = sup_kmer_word_c - 1; // Index of the current word being encoded from the label; the encoding is MSB-boundary aligned for now.
    for(std::size_t b_idx = 0; b_idx < len; b_idx += 32)
    {
        label_buf[label_off + word_idx] = Kmer_Utility::encode<32>(seq + b_idx);
        word_idx--;
    }

    size_++;
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::append(const Super_Kmer_Chunk& c, const std::size_t l, const std::size_t r)
{
    assert(r >= l);

    const auto n = r - l;
    assert(size() + n <= cap_);

    std::memcpy(att_buf.data() + size(), c.att_buf.data() + l, n * sizeof(attribute_t));
    std::memcpy(label_buf.data() + label_units(), c.label_buf.data() + l * sup_kmer_word_c, n * sup_kmer_word_c * sizeof(label_unit_t));

    size_ += n;
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::append(const Super_Kmer_Chunk& c)
{
    append(c, 0, c.size());
}

}



#endif
