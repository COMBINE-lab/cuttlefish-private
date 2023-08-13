
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

        while(E.read_column_buffered(j, buf))
            for(const auto& e : buf)
            {
                auto it = M.find(e.y());
                if(it == M.end())
                    M.emplace(e.y(), Other_End(e.x(), e.s_x(), e.x_is_phi(), e.w(), false));
                else
                {
                    auto& w = it->second;

                    if(w.is_phi() && e.x_is_phi())
                    {
                        form_meta_vertex(e.y(), j, e.s_y(), e.w(), w.w());
                        continue;
                    }

                    if(w.in_same_part())
                    {
                        assert(!w.is_phi());
                        if(e.y() < w.v())
                            D_c.emplace_back(e.y(), inv_side(e.s_y()), w.v(), w.s_v(), w.w(), 0, 0, false, false, side_t::unspecified);

                        w = Other_End(e.x(), e.s_x(), e.x_is_phi(), e.w(), false);
                        continue;
                    }

                    E.add(e.x(), e.s_x(), w.v(), w.s_v(), e.w() + w.w(), 0, 0, e.x_is_phi(), w.is_phi());
                }
            }


        for(const auto& e : D_c)
        {
            const auto w_uv = e.w();
            const auto& m_x = M[e.x()];
            const auto& m_y = M[e.y()];

            if(m_x.is_phi() && m_y.is_phi())
                form_meta_vertex(e.x(), j, inv_side(e.s_x()), m_x.w(), w_uv + m_y.w());
            else
                E.add(m_x.v(), m_x.s_v(), m_y.v(), m_y.s_v(), m_x.w() + w_uv + m_y.w(), 0, 0, m_x.is_phi(), m_y.is_phi());
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
        auto [x, s_x, w_x] = traverse_chain(e.x());
        auto [y, s_y, w_y] = traverse_chain(e.y());

        if(s_x == side_t::unspecified)
            s_x = e.s_x();
        if(s_y == side_t::unspecified)
            s_y = e.s_y();

        const auto w = w_x + e.w() + w_y;
        M[x] = Other_End(y, s_y, false, w, true);
        M[y] = Other_End(x, s_x, false, w, true);

        D_j.emplace_back(x, s_x, y, s_y, w, w == 1 ? e.b() : 0, w == 1 ? e.b_idx() : 0, false, false, side_t::unspecified);
    }


    std::ofstream output(work_path + std::string("D_") + std::to_string(j));
    output.write(reinterpret_cast<const char*>(D_j.data()), D_j.size() * sizeof(Discontinuity_Edge<k>));
    output.close();
}


template <uint16_t k>
std::tuple<Kmer<k>, side_t, weight_t> Discontinuity_Graph_Contractor<k>::traverse_chain(Kmer<k> u) const
{
    side_t s_u = side_t::unspecified;   // Side through which `u` is connected to the chain.
    weight_t w = 0;   // Weight of the traversed chain.
    Kmer<k> p_u = u;  // Last vertex seen before `u`.

    while(true)
    {
        const auto it = M.find(u);
        if(it == M.end())
            break;

        const Other_End& other_end = it->second;
        if(other_end.v() == p_u)
            break;

        p_u = u;
        u = other_end.v();
        s_u = other_end.s_v();
        w += other_end.w();
    }

    return std::make_tuple(u, s_u, w);
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Discontinuity_Graph_Contractor)
