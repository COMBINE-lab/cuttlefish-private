
#include "Subgraphs_Processor.hpp"
#include "Subgraph.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <atomic>


namespace cuttlefish
{

template <uint16_t k>
Subgraphs_Processor<k>::Subgraphs_Processor(const Data_Logistics& logistics, const std::size_t bin_count, Discontinuity_Graph<k>& G, op_buf_list_t& op_buf):
      bin_path_pref(logistics.subgraphs_path())
    , bin_count(bin_count)
    , G(G)
    , trivial_mtig_count_(0)
    , icc_count_(0)
    , op_buf(op_buf)
{}


template <uint16_t k>
void Subgraphs_Processor<k>::process()
{
    Subgraph<k>::init_maps();

    std::atomic_uint64_t solved = 0;

    std::vector<Padded_Data<double>> t_construction(parlay::num_workers(), 0);  // Timer-per-worker for construction work.
    std::vector<Padded_Data<double>> t_contraction(parlay::num_workers(), 0);   // Timer-per-worker for contraction work.

    const auto process_subgraph =
        [&](const std::size_t bin_id)
        {
            Subgraph<k> sub_dBG(bin_path_pref, bin_id, G, op_buf[parlay::worker_id()].data());

            const auto t_0 = timer::now();
            sub_dBG.construct();
            // sub_dBG.construct_loop_filtered();
            const auto t_1 = timer::now();
            sub_dBG.contract();
            const auto t_2 = timer::now();

            trivial_mtig_count_ += sub_dBG.trivial_mtig_count();
            icc_count_ += sub_dBG.icc_count();

            if(++solved % 8 == 0)
                std::cerr << "\rSolved " << solved << " subgraphs.";

            t_construction[parlay::worker_id()].data() += timer::duration(t_1 - t_0);
            t_contraction[parlay::worker_id()].data()  += timer::duration(t_2 - t_0);
        };

    parlay::parallel_for(0, bin_count, process_subgraph, 1);
    std::cerr << "\n";

    Subgraph<k>::free_maps();

    const auto sum_time = [&](const std::vector<Padded_Data<double>>& T)
        { double t = 0; std::for_each(T.cbegin(), T.cend(), [&t](const auto v){ t += v.data(); }); return t; };
    std::cerr << "Total work in graph construction: " << sum_time(t_construction) << " (s).\n";
    std::cerr << "Total work in graph contraction:  " << sum_time(t_contraction) << " (s).\n";
}


template <uint16_t k>
uint64_t Subgraphs_Processor<k>::trivial_mtig_count() const
{
    return trivial_mtig_count_;
}


template <uint16_t k>
uint64_t Subgraphs_Processor<k>::icc_count() const
{
    return icc_count_;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraphs_Processor)
