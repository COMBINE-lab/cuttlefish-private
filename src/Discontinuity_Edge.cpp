
#include "Discontinuity_Edge.hpp"
#include "globals.hpp"


namespace cuttlefish
{

template <uint16_t k> const Kmer<k> Discontinuity_Edge<k>::phi_(phi_label);

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Discontinuity_Edge)
