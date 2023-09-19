
#include "Discontinuity_Graph_Contractor.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <fstream>
#include <algorithm>


namespace cuttlefish
{

template <uint16_t k>
Discontinuity_Graph_Contractor<k>::Discontinuity_Graph_Contractor(Edge_Matrix<k>& E, const std::size_t n, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v, const std::string& temp_path):
      E(E)
    , n_(n)
    , P_v(P_v)
    , P_v_local(parlay::num_workers())
    , work_path(temp_path)
    , M(static_cast<std::size_t>((n_ / E.vertex_part_count()) * 1.1))
    , D_c(parlay::num_workers())
{}


template <uint16_t k>
void Discontinuity_Graph_Contractor<k>::contract()
{
    // Debug
    double map_clr_time = 0;    // Time taken to clear the hash map.
    double edge_proc_time = 0;  // Time taken to process the non-diagonal edges.
    double diag_elim_time = 0;  // Time taken to eliminate the compressed diagonal chains.
    double diag_cont_time = 0;  // Time taken to compress the diagonal chains.
    double v_info_cp_time = 0;  // Time taken to copy worker-local vertex path-info to global repo.

    // TODO: document the phases.

    buf.resize(E.max_block_size());
    for(auto j = E.vertex_part_count(); j >= 1; --j)
    {
        auto t_s = now();
        M.clear();
        auto t_e = now();
        map_clr_time += duration(t_e - t_s);

        std::for_each(D_c.begin(), D_c.end(), [](auto& v){ v.data().clear(); });
        std::for_each(P_v_local.begin(), P_v_local.end(), [](auto& v){ v.data().clear(); });

        std::cerr << "\rPart: " << j;

        t_s = now();
        contract_diagonal_block(j);
        t_e = now();
        diag_cont_time += duration(t_e - t_s);

        const auto process_non_diagonal_edge = [&](const std::size_t idx)
        {
            const auto& e = buf[idx];
            Other_End* p_z;

            assert(E.partition(e.y()) == j);
            assert(e.x_is_phi() || E.partition(e.x()) < j);

            if(M.insert(e.y(), Other_End(e.x(), e.s_x(), e.x_is_phi(), e.w(), false), p_z))
            {
                assert(M.find(e.y()));
                return;
            }

            assert(M.find(e.y()));
            auto& z = *p_z;
            if(z.is_phi() && e.x_is_phi())
            {
                form_meta_vertex(e.y(), j, e.s_y(), e.w(), z.w());
                return;
            }

            if(z.in_same_part())
            {
                assert(!z.is_phi());
                assert(M.find(z.v()));
                assert(E.partition(z.v()) == j);
                if(e.y() < z.v())
                    D_c[parlay::worker_id()].data().emplace_back(e.y(), inv_side(e.s_y()), z.v(), z.s_v(), z.w(), 0, 0, false, false, side_t::unspecified);

                z = Other_End(e.x(), e.s_x(), e.x_is_phi(), e.w(), false);
                return;
            }

            // TODO: add edges in a lock-free manner, accumulating new edges in thread-local buffers and copying to the global repo afterwards.
            E.add(e.x(), e.s_x(), z.v(), z.s_v(), e.w() + z.w(), 0, 0, e.x_is_phi(), z.is_phi());
        };

        while(true)
        {
            t_s = now();
            if(!E.read_column_buffered(j, buf))
                break;
            t_e = now();
            edge_read_time += duration(t_e - t_s);

            t_s = now();
            parlay::parallel_for(0, buf.size(), process_non_diagonal_edge, buf.size() / parlay::num_workers());
            t_e = now();
            edge_proc_time += duration(t_e - t_s);
        }


        t_s = now();
        for(const auto& v : D_c)
            for(const auto& e : v.data())
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

        t_e = now();
        diag_elim_time += duration(t_e - t_s);

        t_s = now();
        for(const auto& v : P_v_local)
        {
            P_v[j].add(v.data().data(), v.data().size());
            meta_v_c += v.data().size();
        }

        t_e = now();
        v_info_cp_time += duration(t_e - t_s);
    }

    std::cerr << "\n";


    std::cerr << "Formed " << meta_v_c << " meta-vertices.\n";
    std::cerr << "Map clearing time: " << map_clr_time << ".\n";
    std::cerr << "Edges reading time: " << edge_read_time << ".\n";
    std::cerr << "Non-diagonal edges contraction time: " << edge_proc_time << ".\n";
    std::cerr << "Meta-vertex copy time: " << v_info_cp_time << ".\n";
    std::cerr << "Diagonal-contraction time: " << diag_cont_time << ".\n";
    std::cerr << "Diagonal-elimination time: " << diag_elim_time << ".\n";
}


template <uint16_t k>
void Discontinuity_Graph_Contractor<k>::contract_diagonal_block(const std::size_t j)
{
    D_j.clear();
    const auto t_s = now();
    E.read_diagonal_block(j, buf);
    const auto t_e = now();
    edge_read_time += duration(t_e - t_s);

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

        assert(E.partition(u) == j);
        assert(E.partition(v) == j);
        assert(u != v);

        M.insert_overwrite(u, Other_End(v, s_v, false, w, true));
        M.insert_overwrite(v, Other_End(u, s_u, false, w, true));

        assert(M.find(u)); assert(M.find(e.x()));
        assert(M.find(v)); assert(M.find(e.y()));

        D_j.emplace_back(u, s_u, v, s_v, w, w == 1 ? e.b() : 0, w == 1 ? e.b_idx() : 0, false, false, side_t::unspecified);
    }


    std::ofstream output(work_path + std::string("D_") + std::to_string(j));
    output.write(reinterpret_cast<const char*>(D_j.data()), D_j.size() * sizeof(Discontinuity_Edge<k>));
    output.close();
}


}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Discontinuity_Graph_Contractor)
