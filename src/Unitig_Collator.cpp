
#include "Unitig_Collator.hpp"
#include "globals.hpp"

#include <cassert>


namespace cuttlefish
{

template <uint16_t k>
Unitig_Collator<k>::Unitig_Collator(const std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const std::string& output_path, const std::string& temp_path):
      P_e(P_e)
    , output_path(output_path)
    , work_path(temp_path)
{}


template <uint16_t k>
void Unitig_Collator<k>::collate()
{
    const std::size_t unitig_bucket_count = P_e.size() - 1;
    std::size_t e_inf_count = 0;    // Debug

    for(std::size_t b = 1; b <= unitig_bucket_count; ++b)
    {
        M.clear();
        load_path_info(b);
        e_inf_count += p_e_buf.size();
    }


    std::cerr << "Read " << e_inf_count << " edge-info.\n";
}


template <uint16_t k>
void Unitig_Collator<k>::load_path_info(std::size_t b)
{
    P_e[b].load(p_e_buf);
    M.resize(p_e_buf.size());
    for(const auto& p_e : p_e_buf)
    {
        assert(p_e.obj() < M.size());
        M[p_e.obj()] = p_e.path_info();
    }
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Unitig_Collator)
