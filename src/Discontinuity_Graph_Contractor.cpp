
#include "Discontinuity_Graph_Contractor.hpp"
#include "globals.hpp"

#include <fstream>


namespace cuttlefish
{

template <uint16_t k>
Discontinuity_Graph_Contractor<k>::Discontinuity_Graph_Contractor(Edge_Matrix<k>& E, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v, const std::string& temp_path):
      E(E)
    , P_v(P_v)
    , work_path(temp_path)
    , n_(E.size() - E.row_size(0) / 2)  // each separate chain has exactly two ϕ-adjacent edges
    , M(static_cast<std::size_t>((n_ / E.vertex_part_count()) * 1.1))
{}


template <uint16_t k>
void Discontinuity_Graph_Contractor<k>::contract()
{
    // TODO: set `buf` to max-block size of E.
    for(auto j = E.vertex_part_count(); j >= 1; --j)
    {
        M.clear();
        D_c.clear();

        contract_diagonal_block(j);

        Other_End* p_z;
        while(E.read_column_buffered(j, buf))
            for(const auto& e : buf)
            {
                if(M.insert(e.y(), Other_End(e.x(), e.s_x(), e.x_is_phi(), e.w(), false), p_z))
                    continue;

                auto& z = *p_z;
                if(z.is_phi() && e.x_is_phi())
                {
                    form_meta_vertex(e.y(), j, e.s_y(), e.w(), z.w());
                    continue;
                }

                if(z.in_same_part())
                {
                    assert(!z.is_phi());
                    if(e.y() < z.v())
                        D_c.emplace_back(e.y(), inv_side(e.s_y()), z.v(), z.s_v(), z.w(), 0, 0, false, false, side_t::unspecified);

                    z = Other_End(e.x(), e.s_x(), e.x_is_phi(), e.w(), false);
                    continue;
                }

                E.add(e.x(), e.s_x(), z.v(), z.s_v(), e.w() + z.w(), 0, 0, e.x_is_phi(), z.is_phi());
            }


        for(const auto& e : D_c)
        {
            assert(M.find(e.x()));
            assert(M.find(e.y()));

            const auto w_xy = e.w();
            const auto& m_x = *M.find(e.x());
            const auto& m_y = *M.find(e.y());

            if(m_x.is_phi() && m_y.is_phi())
                form_meta_vertex(e.x(), j, inv_side(e.s_x()), m_x.w(), w_xy + m_y.w());
            else
                E.add(m_x.v(), m_x.s_v(), m_y.v(), m_y.s_v(), m_x.w() + w_xy + m_y.w(), 0, 0, m_x.is_phi(), m_y.is_phi());
        }
    }


    std::cerr << "Formed " << meta_v_c << " meta-vertices.\n";
}


template <uint16_t k>
void Discontinuity_Graph_Contractor<k>::contract_diagonal_block(const std::size_t j)
{
    D_j.clear();
    E.read_diagonal_block(j, buf);
    for(const auto& e : buf)
    {
        const Other_End* const end_x = M.find(e.x());
        const Other_End* const end_y = M.find(e.y());

        const auto& u  = (end_x ? end_x->v() : e.x());
        const auto s_u = (end_x ? end_x->s_v() : e.s_x());
        const auto w_u = (end_x ? end_x->w() : 0);
        const auto& v  = (end_y ? end_y->v() : e.y());
        const auto s_v = (end_y ? end_y->s_v() : e.s_y());
        const auto w_v = (end_y ? end_y->w() : 0);
        const auto w   = w_u + e.w() + w_v;

        M.insert_overwrite(u, Other_End(v, s_v, false, w, true));
        M.insert_overwrite(v, Other_End(u, s_u, false, w, true));

        D_j.emplace_back(u, s_u, v, s_v, w, w == 1 ? e.b() : 0, w == 1 ? e.b_idx() : 0, false, false, side_t::unspecified);
    }


    std::ofstream output(work_path + std::string("D_") + std::to_string(j));
    output.write(reinterpret_cast<const char*>(D_j.data()), D_j.size() * sizeof(Discontinuity_Edge<k>));
    output.close();
}


}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Discontinuity_Graph_Contractor)
