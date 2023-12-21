
#include "Discontinuity_Graph.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{

template <uint16_t k>
Discontinuity_Graph<k>::Discontinuity_Graph(const std::size_t part_count, const std::size_t lmtig_bucket_count, const std::string& work_path):
      E_(part_count, work_path + "E_")
    , lmtigs(work_path + "lmtig", lmtig_bucket_count, parlay::num_workers())
{}


template <uint16_t k>
void Discontinuity_Graph<k>::close_lmtig_stream()
{
    lmtigs.close();
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Discontinuity_Graph)
