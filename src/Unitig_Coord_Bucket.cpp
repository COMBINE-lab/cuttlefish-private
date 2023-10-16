
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

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Unitig_Coord_Bucket)
