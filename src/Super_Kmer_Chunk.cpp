
#include "Super_Kmer_Chunk.hpp"


namespace cuttlefish
{

template <bool Colored_>
Super_Kmer_Chunk<Colored_>::Super_Kmer_Chunk(const uint16_t k, const uint16_t l):
      max_sup_kmer_len(2 * k - l + 2)
    , sup_kmer_word_c((max_sup_kmer_len + 31) / 32)
{}


template <bool Colored_>
void Super_Kmer_Chunk<Colored_>::clear()
{
    att_buf.clear();
    label_buf.clear();
}


template <bool Colored_>
void Super_Kmer_Chunk<Colored_>::serialize(std::ofstream& os) const
{
    os.write(reinterpret_cast<const char*>(att_buf.data()), att_buf.size() * sizeof(attribute_t));
    os.write(reinterpret_cast<const char*>(label_buf.data()), label_buf.size() * sizeof(label_unit_t));
}


template <bool Colored_>
void Super_Kmer_Chunk<Colored_>::deserialize(std::ifstream& is, const std::size_t sz)
{
    att_buf.resize(sz), label_buf.resize(sz * sup_kmer_word_c);

    is.read(reinterpret_cast<char*>(att_buf.data()), sz * sizeof(attribute_t));
    is.read(reinterpret_cast<char*>(label_buf.data()), label_buf.size() * sizeof(label_unit_t));
}

}



// Template-instantiations for the required instances.
template class cuttlefish::Super_Kmer_Chunk<false>;
