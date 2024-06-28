
#include "Unitig_Scratch.hpp"


template <uint16_t k>
Unitig_Scratch<k>::Unitig_Scratch()
{
    label_.reserve(BUFF_SZ + k - 1),
    hash_.reserve(BUFF_SZ);
}


template <uint16_t k>
void Unitig_Scratch<k>::mark_cycle()
{
    is_cycle_ = true;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Unitig_Scratch)
