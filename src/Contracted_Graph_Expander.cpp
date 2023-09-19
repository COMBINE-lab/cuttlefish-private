
#include "Contracted_Graph_Expander.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <fstream>
#include <filesystem>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k>
Contracted_Graph_Expander<k>::Contracted_Graph_Expander(const Edge_Matrix<k>& E, const std::size_t n, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>>& P_v, std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const std::string& temp_path):
      E(E)
    , n_(n)
    , P_v(P_v)
    , P_e(P_e)
    , work_path(temp_path)
    , M(static_cast<std::size_t>((n_ / E.vertex_part_count()) * 1.25))
{}


template <uint16_t k>
void Contracted_Graph_Expander<k>::expand()
{
    // Debug
    double map_clr_time = 0;    // Time taken to clear the hash map.
    double diag_exp_time = 0;   // Time taken to expand the compressed diagonal chains.
    double edge_proc_time = 0;  // Time taken to process the non-diagonal edges.
    double spec_case_time = 0;  // Time taken to process the special edge-blocks.

    std::vector<Discontinuity_Edge<k>> buf; // Buffer to read-in edges from the edge-matrix.

    buf.resize(E.max_block_size());
    p_v_buf.reserve(static_cast<std::size_t>((n_ / E.vertex_part_count()) * 1.25)); // TODO: check the maximum size empirically, as repetitions are possible in the info-buckets.
    for(std::size_t i = 1; i <= E.vertex_part_count(); ++i)
    {
        auto t_s = now();
        M.clear();
        auto t_e = now();
        map_clr_time += duration(t_e - t_s);

        load_path_info(i);

        t_s = now();
        expand_diagonal_block(i);
        t_e = now();
        diag_exp_time += duration(t_e - t_s);

        Path_Info<k> x_inf; // Computed (earlier) path-information of the endpoint of the edge that is in the current partition.
        while(true)
        {
            t_s = now();
            if(!E.read_row_buffered(i, buf))
                break;
            t_e = now();
            edge_read_time += duration(t_e - t_s);

            t_s = now();
            for(const auto& e : buf)
            {
                const bool has_x = M.find(e.x(), x_inf);
                assert(has_x);
                (void)has_x;

                const auto y_inf = infer(x_inf, e.s_x(), e.s_y(), e.w());
                const auto j = E.partition(e.y());  // TODO: consider obtaining this info during the edge-reading process.
                assert(j > i);
                P_v[j].emplace(e.y(), y_inf.p(), y_inf.r(), y_inf.o());

                if(e.w() == 1)
                {
                    add_edge_path_info(e, x_inf, y_inf);
                    og_edge_c++;
                }
            }
            t_e = now();
            edge_proc_time += duration(t_e - t_s);
        }


        // TODO: consider making the following two blocks more efficient.

        t_s = now();
        E.read_diagonal_block(i, buf);
        t_e = now();
        edge_read_time += duration(t_e - t_s);

        t_s = now();
        for(const auto& e : buf)
            if(e.w() == 1)
            {
                assert(M.find(e.x()) && M.find(e.y()));
                add_edge_path_info(e, *M.find(e.x()), *M.find(e.y()));
                og_edge_c++;
            }
        t_e = now();
        spec_case_time += duration(t_e - t_s);

        t_s = now();
        E.read_block(0, i, buf);
        t_e = now();
        edge_read_time += duration(t_e - t_s);

        t_s = now();
        for(const auto& e : buf)
            if(e.w() == 1)
            {
                assert(M.find(e.y()));
                add_edge_path_info(e, *M.find(e.y()));
                og_edge_c++;
            }
        t_e = now();
        spec_case_time += duration(t_e - t_s);
    }


    std::cerr << "Encountered " << og_edge_c << " original edges\n";
    std::cerr << "Map clearing time: " << map_clr_time << ".\n";
    std::cerr << "Edges reading time: " << edge_read_time << ".\n";
    std::cerr << "Vertex-info load time: " << p_v_load_time << ".\n";
    std::cerr << "Map filling time: " << map_fill_time << ".\n";
    std::cerr << "Non-diagonal blocks expansion time: " << edge_proc_time << ".\n";
    std::cerr << "Diagonal block expansion time: " << diag_exp_time << ".\n";
    std::cerr << "Special case time: " << spec_case_time << ".\n";
}


template <uint16_t k>
void Contracted_Graph_Expander<k>::load_path_info(const std::size_t i)
{
    auto t_s = now();
    P_v[i].load(p_v_buf);
    auto t_e = now();
    p_v_load_time += duration(t_e - t_s);

    t_s = now();
    const auto load_vertex_info = [&](const std::size_t idx)
    {
        const auto& p_v = p_v_buf[idx];
        if(!M.insert(p_v.obj(), p_v.path_info()))
            assert(*M.find(p_v.obj()) == p_v.path_info());
    };

    parlay::parallel_for(0, p_v_buf.size(), load_vertex_info, p_v_buf.size() / parlay::num_workers());

    t_e = now();
    map_fill_time += duration(t_e - t_s);
}


template <uint16_t k>
void Contracted_Graph_Expander<k>::expand_diagonal_block(const std::size_t i)
{
    const std::string d_i_path(work_path + std::string("D_") + std::to_string(i));
    std::error_code ec;
    const auto file_sz = std::filesystem::file_size(d_i_path, ec);
    D_i.resize(file_sz / sizeof(Discontinuity_Edge<k>));

    const auto t_s = now();
    std::ifstream input(d_i_path);
    input.read(reinterpret_cast<char*>(D_i.data()), file_sz);
    input.close();
    const auto t_e = now();
    edge_read_time += duration(t_e - t_s);

    Path_Info<k> y_inf;
    // In reverse order of the newly introduced diagonal-edges to always ensure one endpoint having path-info ready.
    for(auto d_it = D_i.rbegin(); d_it != D_i.rend(); ++d_it)
    {
        const auto& e = *d_it;
        assert(M.find(e.x()) || M.find(e.y()));

        if(!M.find(e.y(), y_inf))
        {
            assert(M.find(e.x()));
            M.insert(e.y(), infer(*M.find(e.x()), e.s_x(), e.s_y(), e.w()));
        }
        else
        {
            assert(M.find(e.y()));
            M.insert(e.x(), infer(y_inf, e.s_y(), e.s_x(), e.w()));
        }
    }
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Contracted_Graph_Expander)
