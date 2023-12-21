
#include "Subgraphs_Processor.hpp"
#include "Subgraph.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{

template <uint16_t k>
Subgraphs_Processor<k>::Subgraphs_Processor(const std::string& bin_path_pref, const std::size_t bin_count, Discontinuity_Graph<k>& G, op_buf_list_t& op_buf):
      bin_path_pref(bin_path_pref)
    , bin_count(bin_count)
    , work_path(bin_path_pref)
    , G(G)
    , op_buf(op_buf)
{}


template <uint16_t k>
void Subgraphs_Processor<k>::process()
{
    std::atomic_uint64_t trivial_mtig_count = 0;

    const auto process_subgraph =
        [&](const std::size_t bin_id)
        {
            Subgraph<k> sub_dBG(bin_path_pref, bin_id, G, op_buf[parlay::worker_id()].data());
            sub_dBG.construct();
            sub_dBG.contract();

            trivial_mtig_count += sub_dBG.trivial_mtig_count();
        };

    parlay::parallel_for(0, bin_count, process_subgraph);

    std::cerr << "Trivial maximal unitig count: " << trivial_mtig_count << "\n";
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraphs_Processor)
