
#include "Subgraphs_Manager.hpp"
#include "Subgraph.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <atomic>
#include <iostream>
#include <cstdlib>


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Subgraphs_Manager<k, Colored_>::Subgraphs_Manager(const Data_Logistics& logistics, const std::size_t graph_count, const uint16_t l, Discontinuity_Graph<k>& G, op_buf_list_t& op_buf):
      path_pref(logistics.subgraphs_path())
    , graph_count(graph_count)
    , l(l)
    , G(G)
    , trivial_mtig_count_(0)
    , icc_count_(0)
    , op_buf(op_buf)
{
    if((graph_count & (graph_count - 1)) != 0)
    {
        std::cerr << "Subgraph count needs to be a power of 2. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    subgraph_bucket.reserve(graph_count);
    for(std::size_t g_id = 0; g_id < graph_count; ++g_id)
        subgraph_bucket.emplace_back(bucket_t(k, l, path_pref + "_" + std::to_string(g_id)));
}


template <uint16_t k, bool Colored_>
void Subgraphs_Manager<k, Colored_>::process()
{
    Subgraph<k>::init_maps();

    std::atomic_uint64_t solved = 0;

    std::vector<Padded_Data<double>> t_construction(parlay::num_workers(), 0);  // Timer-per-worker for construction work.
    std::vector<Padded_Data<double>> t_contraction(parlay::num_workers(), 0);   // Timer-per-worker for contraction work.
    std::vector<Padded_Data<std::size_t>> max_size(parlay::num_workers(), 0);   // Largest graph size processed per worker.

    const auto process_subgraph =
        [&](const std::size_t graph_id)
        {
            Subgraph<k> sub_dBG(path_pref, graph_id, G, op_buf[parlay::worker_id()].data());

            const auto t_0 = timer::now();
            sub_dBG.construct();
            // sub_dBG.construct_loop_filtered();
            const auto t_1 = timer::now();
            sub_dBG.contract();
            const auto t_2 = timer::now();

            auto& max_sz = max_size[parlay::worker_id()].data();
            max_sz = std::max(max_sz, sub_dBG.size());

            trivial_mtig_count_ += sub_dBG.trivial_mtig_count();
            icc_count_ += sub_dBG.icc_count();

            if(++solved % 8 == 0)
                std::cerr << "\rSolved " << solved << " subgraphs.";

            t_construction[parlay::worker_id()].data() += timer::duration(t_1 - t_0);
            t_contraction[parlay::worker_id()].data()  += timer::duration(t_2 - t_1);
        };

    parlay::parallel_for(0, graph_count, process_subgraph, 1);
    std::cerr << "\n";

    Subgraph<k>::free_maps();

    const auto sum_time = [&](const std::vector<Padded_Data<double>>& T)
        { double t = 0; std::for_each(T.cbegin(), T.cend(), [&t](const auto& v){ t += v.data(); }); return t; };
    std::cerr << "Total work in graph construction: " << sum_time(t_construction) << " (s).\n";
    std::cerr << "Total work in graph contraction:  " << sum_time(t_contraction) << " (s).\n";
    std::cerr << "Largest graph size: " <<
        [&](){  std::size_t max_sz = 0;
                std::for_each(max_size.cbegin(), max_size.cend(), [&](const auto& v){ max_sz = std::max(max_sz, v.data()); });
                return max_sz; }() << ".\n";
}


template <uint16_t k, bool Colored_>
uint64_t Subgraphs_Manager<k, Colored_>::trivial_mtig_count() const
{
    return trivial_mtig_count_;
}


template <uint16_t k, bool Colored_>
uint64_t Subgraphs_Manager<k, Colored_>::icc_count() const
{
    return icc_count_;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Subgraphs_Manager)
