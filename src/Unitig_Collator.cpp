
#include "Unitig_Collator.hpp"
#include "Unitig_File.hpp"
#include "Kmer.hpp"
#include "dBG_Utilities.hpp"
#include "globals.hpp"

#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <algorithm>
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

    typedef Kmer<k> path_id_t;
    typedef std::tuple<weight_t, std::string, side_t> lmtig_info_t;
    typedef std::pair<path_id_t, lmtig_info_t> kv_t;    // (p, <r, u ,o>)
    std::vector<kv_t> kv_store;

    std::ofstream output(output_path);
    std::string unitig;
    for(std::size_t b = 1; b <= unitig_bucket_count; ++b)
    {
        M.clear();
        load_path_info(b);

        Unitig_File_Reader unitig_reader(work_path + std::string("lmutig_") + std::to_string(b));
        uni_idx_t uni_idx = 0;
        while(unitig_reader.read_next_unitig(unitig))
        {
            assert(uni_idx < M.size());
            const auto p_e = M[uni_idx];

            kv_store.emplace_back(p_e.p(), lmtig_info_t(p_e.r(), unitig, p_e.o()));
            uni_idx++;
        }
    }

    std::cerr << "KV-store has " << kv_store.size() << " pairs.\n";
    std::sort(kv_store.begin(), kv_store.end());


    std::string max_unitig, max_u_rc;
    std::size_t i, j;
    std::size_t max_path_sz = 0;
    std::size_t min_sz = std::numeric_limits<std::size_t>::max();
    std::size_t max_sz = 0;
    std::size_t max_lmtig_sz = 0;
    for(i = 0; i < kv_store.size(); i = j)
    {
        max_unitig.clear();

        const std::size_t s = i;
        std::size_t e;
        for(e = s + 1; e < kv_store.size() && kv_store[e].first == kv_store[s].first; ++e);

        if(e - s == 2)
        {
            auto& [r0, u0, o0] = kv_store[s].second;
            if(o0 == side_t::front)
                reverse_complement(u0);

            auto& [r1, u1, o1] = kv_store[s + 1].second;
            if(o1 == side_t::front)
                reverse_complement(u1);

            reverse_complement(u1);

            max_unitig = u0 + std::string(u1.begin() + k, u1.end());

            j = e;
            max_lmtig_sz = std::max(max_lmtig_sz, std::max(u0.size(), u1.size()));
        }
        else
            for(j = i; j < kv_store.size() && kv_store[j].first == kv_store[i].first; ++j)
            {
                auto& [r, u, o] = kv_store[j].second;
                if(o == side_t::front)
                    reverse_complement(u);

                max_unitig += (max_unitig.empty() ? u : std::string(u.begin() + k, u.end()));

                max_lmtig_sz = std::max(max_lmtig_sz, u.size());
            }

        max_u_rc = max_unitig;
        reverse_complement(max_u_rc);
        if(max_unitig > max_u_rc)
            max_unitig = max_u_rc;

        output << max_unitig << "\n";

        max_path_sz = std::max(max_path_sz, j - i);
        min_sz = std::min(min_sz, max_unitig.size());
        max_sz = std::max(max_sz, max_unitig.size());
    }

    std::size_t mu_tig = 0;
    Unitig_File_Reader mu_tig_reader(work_path + std::string("mutig"));
    while(mu_tig_reader.read_next_unitig(max_unitig))
    {
        max_u_rc = max_unitig;
        reverse_complement(max_u_rc);
        if(max_unitig > max_u_rc)
            max_unitig = max_u_rc;

        mu_tig++;
        output << max_unitig << "\n";

        min_sz = std::min(min_sz, max_unitig.size());
        max_sz = std::max(max_sz, max_unitig.size());
    }

    output.close();

    std::cerr << "Read " << mu_tig << " trivially maximal unitigs.\n";


    std::cerr << "Longest path length in discontinuity graph: " << max_path_sz << "\n";
    std::cerr << "Minimum maximal-unitig size: " << min_sz << "\n";
    std::cerr << "Maximum maximal-unitig size: " << max_sz << "\n";
    std::cerr << "Maximum locally-maximal unitig size: " << max_lmtig_sz << "\n";
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
