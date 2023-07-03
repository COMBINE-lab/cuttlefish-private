
#include "Unitig_Collator.hpp"
#include "globals.hpp"


namespace cuttlefish
{

template <uint16_t k>
Unitig_Collator<k>::Unitig_Collator(const std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const std::string& temp_path):
      P_e(P_e)
    , work_path(temp_path)
{}


template <uint16_t k>
void Unitig_Collator<k>::collate()
{
    const std::size_t unitig_bucket_count = P_e.size() - 1;

    for(std::size_t i = 1; i <= unitig_bucket_count; ++i)
    {
        M.clear();
        load_path_info(i);
    }
}


template <uint16_t k>
void Unitig_Collator<k>::load_path_info(std::size_t b)
{
    P_e[b].load(p_e_buf);
    M.resize(p_e_buf.size());
    for(const auto& p_e : p_e_buf)
        M[p_e.obj()] = p_e.path_info();
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Unitig_Collator)
