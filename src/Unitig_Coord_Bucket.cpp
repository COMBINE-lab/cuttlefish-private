
#include "Unitig_Coord_Bucket.hpp"
#include "globals.hpp"


namespace cuttlefish
{

template <uint16_t k>
Unitig_Coord_Bucket<k>::Unitig_Coord_Bucket(const std::string& path_pref):
      path_pref(path_pref)
    , coord_bucket(path_pref + ".coord")
    , label_bucket(path_pref + ".label")
    , size_(0)
    , label_len_(0)
{}


template <uint16_t k>
std::size_t Unitig_Coord_Bucket<k>::load_coords(Unitig_Coord<k>* const buf) const
{
    const auto b_sz = coord_bucket.load(buf);
    assert(b_sz == size_);

    return b_sz;
}


template <uint16_t k>
std::size_t Unitig_Coord_Bucket<k>::load_labels(char* const buf) const
{
    const auto len = label_bucket.load(buf);
    assert(len == label_len_);

    return len;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Unitig_Coord_Bucket)
