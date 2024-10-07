
#include "Subgraphs_Manager.hpp"
#include "Subgraph.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <atomic>
#include <cstring>
#include <limits>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Subgraphs_Manager<k, Colored_>::Subgraphs_Manager(const Data_Logistics& logistics, const std::size_t graph_count_, const uint16_t l, Discontinuity_Graph<k, Colored_>& G, op_buf_list_t& op_buf):
      path_pref(logistics.subgraphs_path())
    , graph_count_(graph_count_)
    , l(l)
    , HLL(graph_count_)
    , G_(G)
    , trivial_mtig_count_(0)
    , icc_count_(0)
    , op_buf(op_buf)
    , color_path_pref(logistics.output_file_path())
{
    if((graph_count_ & (graph_count_ - 1)) != 0)
    {
        std::cerr << "Subgraph count needs to be a power of 2. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    subgraph_bucket.reserve(graph_count_);
    for(std::size_t g_id = 0; g_id < graph_count_; ++g_id)
        subgraph_bucket.emplace_back(bucket_t(k, l, path_pref + "_" + std::to_string(g_id)));

    // Separation between these two loops is crucial for `lz4stream` to work.

    for(std::size_t g_id = 0; g_id < graph_count_; ++g_id)
        subgraph_bucket[g_id].unwrap().open();
}


template <uint16_t k, bool Colored_>
void Subgraphs_Manager<k, Colored_>::collate_super_kmer_buffers()
{
    parlay::parallel_for(0, graph_count_, [&](const std::size_t g_id)
    {
        subgraph_bucket[g_id].unwrap().collate_buffers();
    });
}


template <uint16_t k, bool Colored_>
void Subgraphs_Manager<k, Colored_>::finalize()
{
    const auto close_bucket = [&](const std::size_t g_id) { subgraph_bucket[g_id].unwrap().close(); };
    parlay::parallel_for(0, graph_count_, close_bucket);
}


template <uint16_t k, bool Colored_>
uint64_t Subgraphs_Manager<k, Colored_>::estimate_size_max() const
{
    uint64_t max_est = 0;
    for(std::size_t g_id = 0; g_id < graph_count_; ++g_id)
        max_est = std::max(max_est, HLL[g_id].unwrap().estimate());

    return max_est;
}


template <uint16_t k, bool Colored_>
void Subgraphs_Manager<k, Colored_>::process()
{
    Subgraphs_Scratch_Space<k, Colored_> subgraphs_space(estimate_size_max() * 1.10);
    force_free(HLL);

    if constexpr(Colored_)
        subgraphs_space.color_repo().init(color_path_pref + ".col");

    std::atomic_uint64_t solved = 0;

    std::vector<Padded<double>> t_construction(parlay::num_workers(), 0);  // Timer-per-worker for construction work.
    std::vector<Padded<double>> t_contraction(parlay::num_workers(), 0);   // Timer-per-worker for contraction work.
    std::vector<Padded<double>> t_color_extract(parlay::num_workers(), 0);  // Timer-per-worker for color-extraction work.
    std::vector<Padded<double>> t_bucket_rm(parlay::num_workers(), 0); // Timer-per-worker for super k-mer bucket removal.
    std::vector<Padded<uint64_t>> max_kmer_count(parlay::num_workers(), 0);    // Largest k-mer bucket size per worker.
    std::vector<Padded<uint64_t>> min_kmer_count(parlay::num_workers(), std::numeric_limits<uint64_t>::max()); // Smallest k-mer bucket size per worker.
    std::vector<Padded<std::size_t>> size(parlay::num_workers(), 0);   // Sum graph size processed per worker.
    std::vector<Padded<std::size_t>> max_size(parlay::num_workers(), 0);   // Largest graph size processed per worker.
    std::vector<Padded<std::size_t>> min_size(parlay::num_workers(), std::numeric_limits<std::size_t>::max()); // Smallest graph size processed per worker.
    std::vector<Padded<std::size_t>> label_sz(parlay::num_workers(), 0);    // Sum label size produced per worker.
    std::vector<Padded<uint64_t>> mtig_count(parlay::num_workers(), 0); // Number of locally maximal unitigs produced per worker.
    std::vector<Padded<std::size_t>> color_shift(parlay::num_workers(), 0); // Number of color-shifting vertices per worker.
    std::vector<Padded<uint64_t>> new_colored_vertex(parlay::num_workers(), 0); // Number of vertices attempting addition to the global color-table per worker.
    std::vector<Padded<uint64_t>> old_colored_vertex(parlay::num_workers(), 0); // Number of vertices with existing colors from the global color-table per worker.
    std::vector<Padded<uint64_t>> color_rel_sorted(parlay::num_workers(), 0);   // Number of color-relationships sorted in color-extraction per worker.

    std::vector<Padded<double[4]>> color_time(parlay::num_workers());   // Time taken in various steps of coloring.
    std::for_each(color_time.begin(), color_time.end(), [&](auto& c){ std::memset(c.unwrap(), 0, 4 * sizeof(double)); });

    const auto process_subgraph =
        [&](const std::size_t graph_id)
        {
            auto& b = subgraph_bucket[graph_id].unwrap();

            const auto t_0 = timer::now();
            Subgraph<k, Colored_> sub_dBG(b, G_, op_buf[parlay::worker_id()].unwrap(), subgraphs_space);
            sub_dBG.construct();
            const auto t_1 = timer::now();
            if constexpr(!Colored_)
                b.remove();
            const auto t_2 = timer::now();
            sub_dBG.contract();  // Perf-diagnose.
            const auto t_3 = timer::now();
            if constexpr(Colored_)
                sub_dBG.extract_new_colors();
            const auto t_4 = timer::now();
            if constexpr(Colored_)
                b.remove();
            const auto t_5 = timer::now();

            auto& max_kmer_c = max_kmer_count[parlay::worker_id()].unwrap();
            auto& min_kmer_c = min_kmer_count[parlay::worker_id()].unwrap();
            auto& sz = size[parlay::worker_id()].unwrap();
            auto& l_sz = label_sz[parlay::worker_id()].unwrap();
            auto& max_sz = max_size[parlay::worker_id()].unwrap();
            auto& min_sz = min_size[parlay::worker_id()].unwrap();
            auto& mtig_c = mtig_count[parlay::worker_id()].unwrap();
            auto& color_shift_c = color_shift[parlay::worker_id()].unwrap();
            auto& v_new_col_c = new_colored_vertex[parlay::worker_id()].unwrap();
            auto& v_old_col_c = old_colored_vertex[parlay::worker_id()].unwrap();
            auto& color_rel_c = color_rel_sorted[parlay::worker_id()].unwrap();
            auto& t_color = color_time[parlay::worker_id()].unwrap();

            max_kmer_c = std::max(max_kmer_c, sub_dBG.kmer_count());
            min_kmer_c = std::min(min_kmer_c, sub_dBG.kmer_count());
            sz += sub_dBG.size();
            max_sz = std::max(max_sz, sub_dBG.size());
            min_sz = std::min(min_sz, sub_dBG.size());
            l_sz += sub_dBG.label_size();
            mtig_c += sub_dBG.mtig_count();
            color_shift_c += sub_dBG.color_shift_count();
            v_new_col_c += sub_dBG.new_colored_vertex();
            v_old_col_c += sub_dBG.old_colored_vertex();
            color_rel_c += sub_dBG.color_rel_sorted();

            t_color[0] += sub_dBG.collect_rels_time();
            t_color[1] += sub_dBG.sort_time();
            t_color[2] += sub_dBG.collect_sets_time();
            t_color[3] += sub_dBG.attach_time();

            trivial_mtig_count_ += sub_dBG.trivial_mtig_count();
            icc_count_ += sub_dBG.icc_count();

            if(++solved % 8 == 0)
                std::cerr << "\rSolved " << solved << " subgraphs.";

            t_construction[parlay::worker_id()].unwrap() += timer::duration(t_1 - t_0);
            t_bucket_rm[parlay::worker_id()].unwrap() += timer::duration(t_2 - t_1) + timer::duration(t_5 - t_4);
            t_contraction[parlay::worker_id()].unwrap()  += timer::duration(t_3 - t_2);
            t_color_extract[parlay::worker_id()].unwrap() += timer::duration(t_4 - t_3);
        };

    parlay::parallel_for(0, graph_count_, process_subgraph, 1);
    std::cerr << "\n";

    if constexpr(Colored_)
        subgraphs_space.color_repo().close();

    const auto sum_time = [&](const std::vector<Padded<double>>& T)
        { double t = 0; std::for_each(T.cbegin(), T.cend(), [&t](const auto& v){ t += v.unwrap(); }); return t; };
    std::cerr << "Total work in graph construction: " << sum_time(t_construction) << " (s).\n";
    std::cerr << "Total work in graph contraction:  " << sum_time(t_contraction) << " (s).\n";
    if constexpr(Colored_)
        std::cerr << "Total work in color extraction:  " << sum_time(t_color_extract) << " (s).\n";
    std::cerr << "Total work in bucket removal:     " << sum_time(t_bucket_rm) << " (s).\n";

    std::cerr << "Maximum k-mer count in bucket: " <<
        [&](){  std::size_t max_sz = 0;
                std::for_each(max_kmer_count.cbegin(), max_kmer_count.cend(), [&](const auto& v){ max_sz = std::max(max_sz, v.unwrap()); });
                return max_sz; }() << ".\n";
    std::cerr << "Minimum k-mer count in bucket: " <<
        [&](){  std::size_t min_sz = std::numeric_limits<uint64_t>::max();
                std::for_each(min_kmer_count.cbegin(), min_kmer_count.cend(), [&](const auto& v){ min_sz = std::min(min_sz, v.unwrap()); });
                return min_sz; }() << ".\n";
    const auto sum_g_sz = [&](){ std::size_t sz = 0;
                std::for_each(size.cbegin(), size.cend(), [&](const auto& v){ sz += v.unwrap(); });
                return sz; }();
    std::cerr << "Sum graph size: " << sum_g_sz << ".\n";
    std::cerr << "Largest graph size: " <<
        [&](){  std::size_t max_sz = 0;
                std::for_each(max_size.cbegin(), max_size.cend(), [&](const auto& v){ max_sz = std::max(max_sz, v.unwrap()); });
                return max_sz; }() << ".\n";
    std::cerr << "Smallest graph size: " <<
        [&](){  std::size_t min_sz = std::numeric_limits<std::size_t>::max();
                std::for_each(min_size.cbegin(), min_size.cend(), [&](const auto& v){ min_sz = std::min(min_sz, v.unwrap()); });
                return min_sz; }() << ".\n";
    std::cerr << "Sum label size: " <<
        [&](){  std::size_t sz = 0;
                std::for_each(label_sz.cbegin(), label_sz.cend(), [&](const auto& v){ sz += v.unwrap(); });
                return sz; }() << ".\n";
    std::cerr << "lm-tig count: " <<
        [&](){  uint64_t c = 0;
                std::for_each(mtig_count.cbegin(), mtig_count.cend(), [&](const auto& v){ c += v.unwrap(); });
                return c; }() << ".\n";
    if constexpr(Colored_)
    {
        const auto color_shift_c = [&](){ std::size_t c = 0;
                    std::for_each(color_shift.cbegin(), color_shift.cend(), [&](const auto& v){ c += v.unwrap(); });
                    return c; }();
        std::cerr << "Color-shifting vertex count: " << color_shift_c << ".\n";
        std::cerr << "Number of vertices attempting new color-extraction: " <<
            [&](){  std::size_t c = 0;
                    std::for_each(new_colored_vertex.cbegin(), new_colored_vertex.cend(), [&](const auto& v){ c += v.unwrap(); });
                    return c; }() << ".\n";
        std::cerr << "Number of vertices with colors already available: " <<
            [&](){  std::size_t c = 0;
                    std::for_each(old_colored_vertex.cbegin(), old_colored_vertex.cend(), [&](const auto& v){ c += v.unwrap(); });
                    return c; }() << ".\n";
        std::cerr << "Number of vertices skipped from color-consideration: " << (sum_g_sz - color_shift_c) << ".\n";
        std::cerr << "Number of (k-mer, source) pairs sorted: " <<
            [&](){  std::size_t c = 0;
                    std::for_each(color_rel_sorted.cbegin(), color_rel_sorted.cend(), [&](const auto& v){ c += v.unwrap(); });
                    return c; }() << ".\n";

        double t[4] = {0, 0, 0, 0};
        std::for_each(color_time.cbegin(), color_time.cend(), [&](const auto& c)
        {
            const auto& c_t = c.unwrap();
            for(std::size_t i = 0; i < 4; ++i)  t[i] += c_t[i];
        });

        std::cerr << "Color count: " << subgraphs_space.color_map().size() << ".\n";
        std::cerr << "Color-repository size: " << subgraphs_space.color_repo().bytes() << " bytes.\n";

        std::cerr << "Color timings:\n";
        std::cerr << "\t Total work in collecting color-relationships:   " << t[0] << "s.\n";
        std::cerr << "\t Total work in semi-sorting color-relationships: " << t[1] << "s.\n";
        std::cerr << "\t Total work in collecting color-sets:            " << t[2] << "s.\n";
        std::cerr << "\t Total work in attaching color-sets to vertices: " << t[3] << "s.\n";
    }
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
