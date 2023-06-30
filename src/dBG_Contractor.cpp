
#include "dBG_Contractor.hpp"
#include "Discontinuity_Graph_Bootstrap.hpp"
#include "Discontinuity_Graph_Contractor.hpp"
#include "Contracted_Graph_Expander.hpp"
#include "globals.hpp"


namespace cuttlefish
{

template <uint16_t k>
dBG_Contractor<k>::dBG_Contractor(const std::size_t part_count, const std::size_t unitig_bucket_count, const std::string& temp_path):
      part_count(part_count)
    , unitig_bucket_count(unitig_bucket_count)
    , work_path(temp_path)
    , E(part_count, temp_path + std::string("E_"))
{
    P_v.emplace_back();
    for(std::size_t j = 1; j <= part_count; ++j)
        P_v.emplace_back(work_path + std::string("P_v_") + std::to_string(j));
}


template <uint16_t k>
void dBG_Contractor<k>::contract(const uint16_t l, const std::string& cdbg_path)
{
    Discontinuity_Graph_Bootstrap<k> dgb(cdbg_path, l, E, work_path, unitig_bucket_count);
    dgb.generate();

    Discontinuity_Graph_Contractor<k> contractor(E, P_v, work_path);
    contractor.contract();

    Contracted_Graph_Expander<k> expander(E, P_v, work_path);
    expander.expand();
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::dBG_Contractor)
