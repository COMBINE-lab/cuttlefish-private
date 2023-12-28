
#include "dBG_Contractor.hpp"
#include "Discontinuity_Graph_Bootstrap.hpp"
#include "Subgraphs_Processor.hpp"
#include "Discontinuity_Graph_Contractor.hpp"
#include "Contracted_Graph_Expander.hpp"
#include "Unitig_Collator.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <chrono>


namespace cuttlefish
{

template <uint16_t k>
dBG_Contractor<k>::dBG_Contractor(const std::size_t subgraph_count, const std::size_t part_count, const std::size_t lmtig_bucket_count, const std::string& output_path, const std::string& temp_path):
      subgraph_count(subgraph_count)
    , vertex_part_count(part_count)
    , lmtig_bucket_count(lmtig_bucket_count)
    , output_path(output_path)
    , work_path(temp_path)
    , G(vertex_part_count, lmtig_bucket_count, work_path)
    , op_buf(parlay::num_workers(), op_buf_t(output_sink.sink()))
{
    // TODO: centralize all the file-naming schemes.
    P_v.reserve(part_count + 1);
    P_v.emplace_back(); // No vertex other than the ϕ-vertex belongs to partition 0.
    for(std::size_t j = 1; j <= part_count; ++j)
        P_v.emplace_back(work_path + std::string("P_v_") + std::to_string(j));

    P_e.reserve(lmtig_bucket_count + 1);
    P_e.emplace_back(); // Using edge-partition 0 with edges that do not have any associated lm-tig (i.e. has weight > 1).
    for(std::size_t b = 1; b <= lmtig_bucket_count; ++b)
        P_e.emplace_back(work_path + std::string("P_e_") + std::to_string(b));
}


template <uint16_t k>
void dBG_Contractor<k>::contract(const uint16_t l, const std::string& cdbg_path)
{
    // TODO: move these utility functionalities out.
    constexpr auto now = std::chrono::high_resolution_clock::now;
    constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };
    const auto t_s = now();


    // Clear the output file and initialize the output sink.
    clear_file(output_path);
    output_sink.init_sink(output_path);


    (void)l, (void)cdbg_path;
    // Discontinuity_Graph_Bootstrap<k> dgb(cdbg_path, l, E, work_path, unitig_bucket_count);
    // dgb.generate();

    Subgraphs_Processor<k> subgraphs(work_path, subgraph_count, G, op_buf);
    subgraphs.process();

    // TODO: fix the following assertion and `n_disc_v` in light of false-phantom edges' absence from edge matrix.
    std::cerr << "Expecting at most " << ((G.E().row_size(0) + G.phantom_edge_upper_bound()) / 2) << " more non-DCC maximal unitigs\n";
    n_disc_v = G.E().size() - G.E().row_size(0) / 2;    // Each separate chain has exactly two ϕ-adjacent edges.

    Discontinuity_Graph_Contractor<k> contractor(G, P_v, work_path);
    contractor.contract();

    G.close_lmtig_stream();

    const auto t_c = now();
    std::cerr << "Discontinuity-graph contraction completed. Time taken: " << duration(t_c - t_s) << " seconds.\n";

    Contracted_Graph_Expander<k> expander(G, P_v, P_e, work_path);
    expander.expand();

    const auto t_e = now();
    std::cerr << "Expansion of contracted graph completed. Time taken: " << duration(t_e - t_c) << " seconds.\n";

    Unitig_Collator<k> collator(P_e, work_path, op_buf);
    collator.par_collate();

    // Flush data and close the output sink.
    parlay::parallel_for(0, parlay::num_workers(),
                        [&](const std::size_t idx){ op_buf[idx].data().close(); }, 1);
    output_sink.close_sink();

    const auto t_uc = now();
    std::cerr << "Unitigs-collation completed. Time taken: " << duration(t_uc - t_e) << " seconds.\n";
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::dBG_Contractor)
