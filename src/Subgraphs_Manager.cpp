
#include "Subgraphs_Manager.hpp"
#include "Subgraph.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <atomic>
#include <limits>
#include <iostream>
#include <cstdlib>
#include <algorithm>


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Subgraphs_Manager<k, Colored_>::Subgraphs_Manager(const Data_Logistics& logistics, const std::size_t graph_count_, const uint16_t l, Discontinuity_Graph<k>& G, op_buf_list_t& op_buf):
      path_pref(logistics.subgraphs_path())
    , graph_count_(graph_count_)
    , l(l)
    , HLL(graph_count_)
    , G(G)
    , trivial_mtig_count_(0)
    , icc_count_(0)
    , op_buf(op_buf)
{
    if((graph_count_ & (graph_count_ - 1)) != 0)
    {
        std::cerr << "Subgraph count needs to be a power of 2. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    subgraph_bucket.reserve(graph_count_);
    for(std::size_t g_id = 0; g_id < graph_count_; ++g_id)
        subgraph_bucket.emplace_back(bucket_t(k, l, path_pref + "_" + std::to_string(g_id)));
}


template <uint16_t k, bool Colored_>
void Subgraphs_Manager<k, Colored_>::finalize()
{
    const auto close_bucket = [&B = subgraph_bucket](const std::size_t g_id) { B[g_id].unwrap().close(); };
    parlay::parallel_for(0, graph_count_, close_bucket);
}


template <uint16_t k, bool Colored_>
uint64_t Subgraphs_Manager<k, Colored_>::estimate_size_max() const
{
    uint64_t max_est = 0;
    for(std::size_t g_id = 0; g_id < graph_count_; ++g_id)
        max_est = std::max(max_est, HLL[g_id].estimate());

    return max_est;
}


template <uint16_t k, bool Colored_>
void Subgraphs_Manager<k, Colored_>::process()
{
    uint64_t max_sub_g_sz = 0;
    std::for_each(HLL.cbegin(), HLL.cend(), [&](const auto& hll){ max_sub_g_sz = std::max(max_sub_g_sz, hll.estimate()); });
    force_free(HLL);

    Subgraphs_Scratch_Space<k> subgraphs_space(max_sub_g_sz * 1.10);
    std::atomic_uint64_t solved = 0;

    std::vector<Padded<double>> t_construction(parlay::num_workers(), 0);  // Timer-per-worker for construction work.
    std::vector<Padded<double>> t_contraction(parlay::num_workers(), 0);   // Timer-per-worker for contraction work.
    std::vector<Padded<double>> t_bucket_rm(parlay::num_workers(), 0); // Timer-per-worker for super k-mer bucket removal.
    std::vector<Padded<uint64_t>> max_kmer_count(parlay::num_workers(), 0);    // Largest k-mer bucket size per worker.
    std::vector<Padded<uint64_t>> min_kmer_count(parlay::num_workers(), std::numeric_limits<uint64_t>::max()); // Smallest k-mer bucket size per worker.
    std::vector<Padded<std::size_t>> size(parlay::num_workers(), 0);   // Sum graph size processed per worker.
    std::vector<Padded<std::size_t>> max_size(parlay::num_workers(), 0);   // Largest graph size processed per worker.
    std::vector<Padded<std::size_t>> min_size(parlay::num_workers(), std::numeric_limits<std::size_t>::max()); // Smallest graph size processed per worker.

    const auto process_subgraph =
        [&](const std::size_t graph_id)
        {
            auto& b = subgraph_bucket[graph_id].unwrap();

            const auto t_0 = timer::now();
            Subgraph<k, false> sub_dBG(b, G, op_buf[parlay::worker_id()].unwrap(), subgraphs_space);
            sub_dBG.construct();
            // sub_dBG.construct_loop_filtered();
            const auto t_1 = timer::now();
            b.remove();
            const auto t_2 = timer::now();
            sub_dBG.contract();  // Perf-diagnose.
            const auto t_3 = timer::now();

            auto& max_kmer_c = max_kmer_count[parlay::worker_id()].unwrap();
            auto& min_kmer_c = min_kmer_count[parlay::worker_id()].unwrap();
            auto& sz = size[parlay::worker_id()].unwrap();
            auto& max_sz = max_size[parlay::worker_id()].unwrap();
            auto& min_sz = min_size[parlay::worker_id()].unwrap();
            max_kmer_c = std::max(max_kmer_c, sub_dBG.kmer_count());
            min_kmer_c = std::min(min_kmer_c, sub_dBG.kmer_count());
            sz += sub_dBG.size();
            max_sz = std::max(max_sz, sub_dBG.size());
            min_sz = std::min(min_sz, sub_dBG.size());

            trivial_mtig_count_ += sub_dBG.trivial_mtig_count();
            icc_count_ += sub_dBG.icc_count();

            if(++solved % 8 == 0)
                std::cerr << "\rSolved " << solved << " subgraphs.";

            t_construction[parlay::worker_id()].unwrap() += timer::duration(t_1 - t_0);
            t_bucket_rm[parlay::worker_id()].unwrap() += timer::duration(t_2 - t_1);
            t_contraction[parlay::worker_id()].unwrap()  += timer::duration(t_3 - t_2);
        };

    parlay::parallel_for(0, graph_count_, process_subgraph, 1);
    std::cerr << "\n";

    const auto sum_time = [&](const std::vector<Padded<double>>& T)
        { double t = 0; std::for_each(T.cbegin(), T.cend(), [&t](const auto& v){ t += v.unwrap(); }); return t; };
    std::cerr << "Total work in graph construction: " << sum_time(t_construction) << " (s).\n";
    std::cerr << "Total work in graph contraction:  " << sum_time(t_contraction) << " (s).\n";
    std::cerr << "Total work in bucket removal:     " << sum_time(t_bucket_rm) << " (s).\n";

    std::cerr << "Maximum k-mer count in bucket: " <<
        [&](){  std::size_t max_sz = 0;
                std::for_each(max_kmer_count.cbegin(), max_kmer_count.cend(), [&](const auto& v){ max_sz = std::max(max_sz, v.unwrap()); });
                return max_sz; }() << ".\n";
    std::cerr << "Minimum k-mer count in bucket: " <<
        [&](){  std::size_t min_sz = std::numeric_limits<uint64_t>::max();
                std::for_each(min_kmer_count.cbegin(), min_kmer_count.cend(), [&](const auto& v){ min_sz = std::min(min_sz, v.unwrap()); });
                return min_sz; }() << ".\n";
    std::cerr << "Sum graph size: " <<
        [&](){  std::size_t sz = 0;
                std::for_each(size.cbegin(), size.cend(), [&](const auto& v){ sz += v.unwrap(); });
                return sz; }() << ".\n";
    std::cerr << "Largest graph size: " <<
        [&](){  std::size_t max_sz = 0;
                std::for_each(max_size.cbegin(), max_size.cend(), [&](const auto& v){ max_sz = std::max(max_sz, v.unwrap()); });
                return max_sz; }() << ".\n";
    std::cerr << "Smallest graph size: " <<
        [&](){  std::size_t min_sz = std::numeric_limits<std::size_t>::max();
                std::for_each(min_size.cbegin(), min_size.cend(), [&](const auto& v){ min_sz = std::min(min_sz, v.unwrap()); });
                return min_sz; }() << ".\n";
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
