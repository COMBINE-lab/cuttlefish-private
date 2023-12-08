
#include "Subgraphs_Processor.hpp"
#include "Subgraph.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{

template <uint16_t k>
Subgraphs_Processor<k>::Subgraphs_Processor(const std::string& bin_path_pref, const std::size_t bin_count, Edge_Matrix<k>& E, op_buf_list_t& op_buf):
      bin_path_pref(bin_path_pref)
    , bin_count(bin_count)
    , E(E)
    , op_buf(op_buf)
{}


template <uint16_t k>
void Subgraphs_Processor<k>::process()
{
    const auto process_subgraph =
      [&](const std::size_t bin_id)
      {
          Subgraph<k> G(bin_path_pref, bin_id, E, op_buf[parlay::worker_id()].data());
          G.construct();
          G.contract();
      };

    parlay::parallel_for(0, bin_count, process_subgraph);
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraphs_Processor)
