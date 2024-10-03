
#include "Super_Kmer_Chunk.hpp"
#include "utility.hpp"

#include <cstdlib>
#include <cassert>


namespace cuttlefish
{

template <bool Colored_>
Super_Kmer_Chunk<Colored_>::Super_Kmer_Chunk(const uint16_t k, const uint16_t l, const std::size_t cap):
      max_sup_kmer_len(2 * (k - 1) - l + 2)
    , sup_kmer_word_c((max_sup_kmer_len + 31) / 32)
    , cap_(cap)
    , size_(0)
    , att_buf(cap_)
    , label_buf(cap * sup_kmer_word_c)
{
    assert(k > l);
    assert(sup_kmer_word_c > 0);
    assert(cap_ > 0 || !Colored_);  // A super k-mer bucket's "safe"-chunk is meaningless in the uncolored case.
}


template <bool Colored_>
Super_Kmer_Chunk<Colored_>::Super_Kmer_Chunk(Super_Kmer_Chunk&& rhs):
      max_sup_kmer_len(rhs.max_sup_kmer_len)
    , sup_kmer_word_c(rhs.sup_kmer_word_c)
    , cap_(rhs.cap_)
    , size_(rhs.size_)
    , att_buf(std::move(rhs.att_buf))
    , label_buf(std::move(rhs.label_buf))
{}


template <bool Colored_>
std::size_t Super_Kmer_Chunk<Colored_>::record_size(const uint16_t k, const uint16_t l)
{
    return sizeof(attribute_t) + (((2 * (k - 1) - l + 2) + 31) / 32) * sizeof(label_unit_t);
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Chunk<false>;
template class cuttlefish::Super_Kmer_Chunk<true>;
