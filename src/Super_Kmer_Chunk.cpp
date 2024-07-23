
#include "Super_Kmer_Chunk.hpp"
#include "utility.hpp"

#include <cstdlib>
#include <iostream>
#include <cassert>


namespace cuttlefish
{

template <bool Colored_>
Super_Kmer_Chunk<Colored_>::Super_Kmer_Chunk(const uint16_t k, const uint16_t l, const std::size_t cap):
      max_sup_kmer_len(2 * (k - 1) - l + 2)
    , sup_kmer_word_c((max_sup_kmer_len + 31) / 32)
    , cap_(cap)
    , size_(0)
    , att_buf(allocate<attribute_t>(cap_))
    , label_buf(allocate<label_unit_t>(cap * sup_kmer_word_c))
{
    assert(k > l);
    assert(sup_kmer_word_c > 0);
    assert(cap_ > 0);
}


template <bool Colored_>
Super_Kmer_Chunk<Colored_>::Super_Kmer_Chunk(Super_Kmer_Chunk&& rhs):
      max_sup_kmer_len(rhs.max_sup_kmer_len)
    , sup_kmer_word_c(rhs.sup_kmer_word_c)
    , cap_(rhs.cap_)
    , size_(rhs.size_)
    , att_buf(rhs.att_buf)
    , label_buf(rhs.label_buf)
{
    // Moved objects are not really moved in C++.
    rhs.att_buf = nullptr;
    rhs.label_buf = nullptr;
}


template <bool Colored_>
Super_Kmer_Chunk<Colored_>::~Super_Kmer_Chunk()
{
    deallocate(att_buf);
    deallocate(label_buf);
}


template <bool Colored_>
std::size_t Super_Kmer_Chunk<Colored_>::record_size(const uint16_t k, const uint16_t l)
{
    return sizeof(attribute_t) + (((2 * (k - 1) - l + 2) + 31) / 32) * sizeof(label_unit_t);
}


template <bool Colored_>
void Super_Kmer_Chunk<Colored_>::serialize(std::ofstream& os) const
{
    os.write(reinterpret_cast<const char*>(att_buf), size() * sizeof(attribute_t));
    os.write(reinterpret_cast<const char*>(label_buf), label_units() * sizeof(label_unit_t));

    if(!os)
    {
        std::cerr << "Serialization of super k-mer chunk of size " << size() << " failed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <bool Colored_>
void Super_Kmer_Chunk<Colored_>::deserialize(std::ifstream& is, const std::size_t sz)
{
    assert(sz <= cap_);

    size_ = sz;

    is.read(reinterpret_cast<char*>(att_buf), size() * sizeof(attribute_t));
    assert(static_cast<std::size_t>(is.gcount()) == size() * sizeof(attribute_t));

    is.read(reinterpret_cast<char*>(label_buf), label_units() * sizeof(label_unit_t));
    assert(static_cast<std::size_t>(is.gcount()) == label_units() * sizeof(label_unit_t));

    if(!is)
    {
        std::cerr << "Deserialization of super k-mer chunk of size " << size() << " failed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Chunk<false>;
template class cuttlefish::Super_Kmer_Chunk<true>;
