
#include "Discontinuity_Graph_Contractor.hpp"


namespace cuttlefish
{

template <uint16_t k>
Discontinuity_Graph_Contractor<k>::Discontinuity_Graph_Contractor(Edge_Matrix<k>& E, const std::string& temp_path):
      E(E)
    , work_path(temp_path)
{
    for(std::size_t i = 0; i <= E.vertex_part_count() + 1; ++i)
        P_v.emplace_back(work_path + std::string("P_v_") + std::to_string(i));
}


template <uint16_t k>
void Discontinuity_Graph_Contractor<k>::contract()
{
    // TODO: move these utility functionalities out.
    constexpr auto now = std::chrono::high_resolution_clock::now;
    constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };
    const auto t_b = now();

    // TODO: set `buf` to max-block size of E.
    for(auto j = E.vertex_part_count(); j >= 1; --j)
    {
        M.clear();
        D.clear();


        contract_diagonal_block(j);

        while(E.read_column_buffered(j, buf))
            for(const auto& e : buf)
            {
                auto it = M.find(e.v());
                if(it == M.end())
                    M.emplace(e.v(), Other_End(e.u(), e.s_u(), e.u_is_phi(), e.w(), false));
                else
                {
                    auto& w = it->second;

                    if(w.is_phi() && e.u_is_phi())
                    {
                        form_meta_vertex(e.v(), j, e.s_v(), e.w(), w.w());
                        continue;
                    }

                    if(w.in_same_part())
                    {
                        if(e.v() < w.v())
                            D.emplace_back(e.v(), w.v(), w.w());

                        w = Other_End(e.u(), e.s_u(), e.u_is_phi(), e.w(), false);
                        continue;
                    }

                    E.add(e.u(), e.s_u(), w.v(), w.s_v(), e.w() + w.w(), 0, e.u_is_phi(), w.is_phi());
                }
            }


        for(const auto& e : D)
        {
            const auto u = std::get<0>(e);
            const auto v = std::get<1>(e);
            const auto w_uv = std::get<2>(e);
            const auto& m_u = M[u];
            const auto& m_v = M[v];

            if(m_u.is_phi() && m_v.is_phi())
                form_meta_vertex(u, j, m_u.s_v(), m_u.w(), w_uv + m_v.w());
            else
                E.add(m_u.v(), m_u.s_v(), m_v.v(), m_v.s_v(), m_u.w() + w_uv + m_v.w(), 0, m_u.is_phi(), m_v.is_phi());
        }


        P_v[j].close(); // TODO: remove
    }


    std::cerr << "Formed " << meta_v_c << " meta-vertices.\n";


    const auto t_e = now();
    std::cerr << "Discontinuity-graph contraction completed. Time taken: " << duration(t_e - t_b) << " seconds.\n";
}


template <uint16_t k>
void Discontinuity_Graph_Contractor<k>::contract_diagonal_block(const std::size_t j)
{
    E.read_diagonal_block(j, buf);
    for(const auto& e : buf)
    {
        auto [u, s_u, w_u] = traverse_chain(e.u());
        auto [v, s_v, w_v] = traverse_chain(e.v());

        if(s_u == side_t::unspecified)
            s_u = e.s_u();
        if(s_v == side_t::unspecified)
            s_v = e.s_v();

        M[u] = Other_End(v, s_v, false, w_u + e.w() + w_v, true);
        M[v] = Other_End(u, s_u, false, w_v + e.w() + w_u, true);
      }
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
