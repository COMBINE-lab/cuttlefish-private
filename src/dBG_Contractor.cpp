
#include "dBG_Contractor.hpp"
#include "Discontinuity_Graph_Bootstrap.hpp"
#include "Discontinuity_Graph_Contractor.hpp"
#include "Contracted_Graph_Expander.hpp"
#include "Unitig_Collator.hpp"
#include "globals.hpp"

#include <chrono>


namespace cuttlefish
{

template <uint16_t k>
dBG_Contractor<k>::dBG_Contractor(const std::size_t part_count, const std::size_t unitig_bucket_count, const std::string& output_path, const std::string& temp_path):
      part_count(part_count)
    , unitig_bucket_count(unitig_bucket_count)
    , output_path(output_path)
    , work_path(temp_path)
    , E(part_count, temp_path + std::string("E_"), true)
    , n_disc_v(E.size() - E.row_size(0) / 2)    // Each separate chain has exactly two ϕ-adjacent edges.
{
    P_v.reserve(part_count + 1);
    P_v.emplace_back();
    for(std::size_t j = 1; j <= part_count; ++j)
        P_v.emplace_back(work_path + std::string("P_v_") + std::to_string(j));

    P_e.reserve(part_count + 1);
    P_e.emplace_back();
    for(std::size_t b = 1; b <= unitig_bucket_count; ++b)
        P_e.emplace_back(work_path + std::string("P_e_") + std::to_string(b));
}


template <uint16_t k>
void dBG_Contractor<k>::contract(const uint16_t l, const std::string& cdbg_path)
{
    // TODO: move these utility functionalities out.
    constexpr auto now = std::chrono::high_resolution_clock::now;
    constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };
    const auto t_s = now();

    (void)l, (void)cdbg_path;
    // Discontinuity_Graph_Bootstrap<k> dgb(cdbg_path, l, E, work_path, unitig_bucket_count);
    // dgb.generate();

    Discontinuity_Graph_Contractor<k> contractor(E, n_disc_v, P_v, work_path);
    contractor.contract();

    const auto t_c = now();
    std::cerr << "Discontinuity-graph contraction completed. Time taken: " << duration(t_c - t_s) << " seconds.\n";

    Contracted_Graph_Expander<k> expander(E, n_disc_v, P_v, P_e, work_path);
    expander.expand();

    const auto t_e = now();
    std::cerr << "Expansion of contracted graph completed. Time taken: " << duration(t_e - t_c) << " seconds.\n";

    Unitig_Collator<k> collator(P_e, output_path, work_path);
    collator.collate();

    const auto t_uc = now();
    std::cerr << "Unitigs-collation completed. Time taken: " << duration(t_uc - t_e) << " seconds.\n";
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::dBG_Contractor)
