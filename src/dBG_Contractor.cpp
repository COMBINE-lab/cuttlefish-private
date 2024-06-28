
#include "dBG_Contractor.hpp"
#include "Discontinuity_Graph_Bootstrap.hpp"
#include "Graph_Partitioner.hpp"
#include "State_Config.hpp"
#include "Subgraphs_Manager.hpp"
#include "Discontinuity_Graph_Contractor.hpp"
#include "Contracted_Graph_Expander.hpp"
#include "Unitig_Collator.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "profile.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{

template <uint16_t k>
dBG_Contractor<k>::dBG_Contractor(const Build_Params& params):
      params(params)
    , logistics(params)
    , G(params.vertex_part_count(), params.lmtig_bucket_count(), logistics)
    , op_buf(parlay::num_workers(), op_buf_t(output_sink.sink()))
{
    cuttlefish::State_Config::set_edge_threshold(params.cutoff());
    std::cerr << "Edge frequency cutoff: " << params.cutoff() << ".\n";

    // TODO: no need to instantiate `P_v` and `P_e` right away—waste of working memory.

    const auto p_v_path_pref = logistics.vertex_path_info_buckets_path();
    P_v.reserve(params.vertex_part_count() + 1);
    P_v.emplace_back(); // No vertex other than the ϕ-vertex belongs to partition 0.
    for(std::size_t j = 1; j <= params.vertex_part_count(); ++j)
        P_v.emplace_back(p_v_bucket_t(p_v_path_pref + "_" + std::to_string(j), 128 * 1024));    // 128 KB for `P_v` buffers.

    const auto p_e_path_pref = logistics.edge_path_info_buckets_path();
    P_e.reserve(params.lmtig_bucket_count() + 1);
    P_e.emplace_back(); // Using edge-partition 0 with edges that do not have any associated lm-tig (i.e. has weight > 1).
    for(std::size_t b = 1; b <= params.lmtig_bucket_count(); ++b)
        P_e.emplace_back(p_e_bucket_t(p_e_path_pref + "_" + std::to_string(b), 128 * 1024));    // 128 KB for `P_e` buffers.
}


template <uint16_t k>
void dBG_Contractor<k>::construct()
{
    // Clear the output file and initialize the output sink.
    const auto op_file_path = logistics.output_file_path();
    clear_file(op_file_path);
    output_sink.init_sink(op_file_path);

    const auto t_0 = timer::now();

    Subgraphs_Manager<k, false> subgraphs(logistics, params.subgraph_count(), params.min_len(), G, op_buf);

    Graph_Partitioner<k, false> super_kmer_splitter(subgraphs, logistics, params.min_len());
    EXECUTE("partition", super_kmer_splitter.partition);
    subgraphs.finalize();

    const auto t_part = timer::now();
    std::cerr << "Sequence splitting into subgraphs completed. Time taken: " << timer::duration(t_part - t_0) << " seconds.\n";

    EXECUTE("subgraphs", subgraphs.process);

    std::cerr << "Trivial maximal unitig count: " << subgraphs.trivial_mtig_count() << ".\n";
    std::cerr << "Trivial ICC count: " << subgraphs.icc_count() << ".\n";

    const auto t_subg = timer::now();
    std::cerr << "Subgraphs construction and contraction completed. Time taken: " << timer::duration(t_subg - t_part) << " seconds.\n";

    std::cerr << "Edge-matrix size: " << G.E().size() << "\n";
    std::cerr << "Phantom edge upper-bound: " << G.phantom_edge_upper_bound() << "\n";
    std::cerr << "Expecting at most " << ((G.E().row_size(0) + G.phantom_edge_upper_bound()) / 2) << " more non-DCC maximal unitigs\n";
    return; // Perf-diagnose

    Discontinuity_Graph_Contractor<k> contractor(G, P_v, logistics);
    EXECUTE("contract", contractor.contract)

    G.close_lmtig_stream();

    const auto t_c = timer::now();
    std::cerr << "Discontinuity-graph contraction completed. Time taken: " << timer::duration(t_c - t_subg) << " seconds.\n";

    Contracted_Graph_Expander<k> expander(G, P_v, P_e, logistics);
    EXECUTE("expand", expander.expand)

    const auto t_e = timer::now();
    std::cerr << "Expansion of contracted graph completed. Time taken: " << timer::duration(t_e - t_c) << " seconds.\n";

    Unitig_Collator<k> collator(P_e, logistics, op_buf, params.gmtig_bucket_count());
    EXECUTE("collate", collator.collate);

    // Flush data and close the output sink.
    parlay::parallel_for(0, parlay::num_workers(),
                        [&](const std::size_t idx){ op_buf[idx].data().close(); }, 1);
    output_sink.close_sink();

    const auto t_uc = timer::now();
    std::cerr << "Unitigs-collation completed. Time taken: " << timer::duration(t_uc - t_e) << " seconds.\n";
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::dBG_Contractor)
