
#include "Subgraphs_Processor.hpp"
#include "Subgraph.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{

template <uint16_t k>
Subgraphs_Processor<k>::Subgraphs_Processor(const std::string& bin_path_pref, const std::size_t bin_count, Edge_Matrix<k>& E):
      bin_path_pref(bin_path_pref)
    , bin_count(bin_count)
    , E(E)
{}


template <uint16_t k>
void Subgraphs_Processor<k>::process()
{
    const auto process_subgraph =
      [&](const std::size_t bin_id)
      {
          Subgraph<k> G(bin_path_pref, bin_id, E);
          G.construct();
          G.contract();
      };

    parlay::parallel_for(0, bin_count, process_subgraph);
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraphs_Processor)
