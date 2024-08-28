
#include "Contracted_Graph_Expander.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Contracted_Graph_Expander<k, Colored_>::Contracted_Graph_Expander(Discontinuity_Graph<k, Colored_>& G, P_v_t& P_v, P_e_t& P_e, const Data_Logistics& logistics):
      G(G)
    , P_v(P_v)
    , P_e(P_e)
    , compressed_diagonal_path(logistics.compressed_diagonal_path())
    , M(G.vertex_part_size_upper_bound())
    , P_v_ip1(parlay::num_workers())
{
    std::cerr << "Hash table capacity during expansion: " << M.capacity() << ".\n";
}


template <uint16_t k, bool Colored_>
void Contracted_Graph_Expander<k, Colored_>::expand()
{
    // Debug
    double map_clr_time = 0;    // Time taken to clear the hash map.
    double diag_exp_time = 0;   // Time taken to expand the compressed diagonal chains.
    double edge_proc_time = 0;  // Time taken to process the non-diagonal edges.
    double spec_case_time = 0;  // Time taken to process the special edge-blocks.

    Buffer<Discontinuity_Edge<k>> buf;  // Buffer to read-in edges from the edge-matrix.

    Buffer<Obj_Path_Info_Pair<Kmer<k>, k>> p_v_buf; // Buffer to read-in path-information of vertices.
    Buffer<Obj_Path_Info_Pair<Kmer<k>, k>> swap_p_v_buf;    // Path-info buffer to be used in every other iteration.
    std::size_t max_p_v_buf_sz = 0;
    std::for_each(P_v.cbegin(), P_v.cend(), [&](const auto& P_v_i){ max_p_v_buf_sz = std::max(max_p_v_buf_sz, P_v_i.unwrap().size()); });
    p_v_buf.reserve_uninit(max_p_v_buf_sz);
    swap_p_v_buf.reserve_uninit(max_p_v_buf_sz);

    const auto max_blk_sz = G.E().max_block_size();
    buf.resize_uninit(max_blk_sz);
    Buffer<Discontinuity_Edge<k>> swap_buf; // Buffer to be used in every other iteration to hide edge-read latency.
    swap_buf.resize_uninit(buf.capacity());

    const auto v_part_c = G.E().vertex_part_count();    // Number of vertex-partitions.
    auto p_v_c_curr = load_path_info(1, p_v_buf);   // Count of vertex path-information instances to process in the current batch.
    std::size_t p_v_c_next = 0; // Count of vertex path-information instances to process in the next batch.
    for(std::size_t i = 1; i <= v_part_c; ++i)
    {
        std::cerr << "\rPart: " << i;

        parlay::par_do3(
            [&]()
            {
                P_v[i].unwrap().remove();
            },
            [&]()
            {
                if(i == G.E().vertex_part_count())
                    return;

                auto& p_v_buf_next = ((i & 1) == 1 ? swap_p_v_buf : p_v_buf);
                p_v_c_next = load_path_info(i + 1, p_v_buf_next);
            }
            ,
            [&]()
            {
                auto t_s = now();
                M.clear();
                auto t_e = now();
                map_clr_time += duration(t_e - t_s);

                const auto& p_v_buf_cur = ((i & 1) == 1 ? p_v_buf : swap_p_v_buf);
                fill_path_info(p_v_buf_cur, p_v_c_curr);
                std::for_each(P_v_ip1.begin(), P_v_ip1.end(), [](auto& v){ v.unwrap().clear(); });

                t_s = now();
                expand_diagonal_block(i);
                t_e = now();
                diag_exp_time += duration(t_e - t_s);


                std::size_t cur_col = i + 1;    // Column of the current block of edges to be processed.
                auto edge_c_curr = (cur_col <= v_part_c ?   // Count of edges to process in the current batch.
                                                    G.E().read_block(i, cur_col, ((cur_col & 1) == 0 ? buf : swap_buf)) : 0);
                std::size_t edge_c_next;    // Count of edges to process in the next batch.
                while(true)
                {
                    assert(edge_c_curr <= max_blk_sz);
                    if(edge_c_curr == 0)
                        break;

                    parlay::par_do3(
                        [&]()
                        {
                            G.E().remove_block(i, cur_col);
                        },
                        [&]()
                        {
                            t_s = now();
                            auto& buf_next = (((cur_col + 1) & 1) == 0 ? buf : swap_buf);   // Buffer for the next batch.
                            edge_c_next = (cur_col < v_part_c ? G.E().read_block(i, cur_col + 1, buf_next) : 0);
                            t_e = now();
                            edge_read_time += duration(t_e - t_s);
                        }
                        ,
                        [&]()
                        {
                            const auto& buf_cur = ((cur_col & 1) == 0 ? buf : swap_buf);    // Buffer for the current batch.
                            const auto process_non_diagonal_edge = [&](const std::size_t idx)
                            {
                                const auto& e = buf_cur[idx];

                                assert(M.find(e.x()));
                                const auto x_inf = *M.find(e.x());  // Path-information of the endpoint of the edge that is in the current partition, `i`.
                                const auto y_inf = infer(x_inf, e.s_x(), e.s_y(), e.w());
                                if(y_inf.r() > 0)   // `e` is not the back-edge to the rank-1 vertex in an ICC.
                                {
                                    const auto j = G.E().partition(e.y());  // TODO: consider obtaining this info during the edge-reading process.
                                    assert(j > i);
                                    if(j > i + 1)
                                        P_v[j].unwrap().emplace(e.y(), y_inf.p(), y_inf.r(), y_inf.o(), y_inf.is_cycle());
                                    else
                                        P_v_ip1[parlay::worker_id()].unwrap().emplace_back(e.y(), y_inf.p(), y_inf.r(), y_inf.o(), y_inf.is_cycle());
                                }

                                if(e.w() == 1)
                                    add_edge_path_info(e, x_inf, y_inf);
                            };

                            t_s = now();
                            parlay::parallel_for(0, edge_c_curr, process_non_diagonal_edge);
                            t_e = now();
                            edge_proc_time += duration(t_e - t_s);
                        }
                    );

                    cur_col++;
                    edge_c_curr = edge_c_next;
                }
            }
        );

        p_v_c_curr = p_v_c_next;

        double diag_blk_time, top_blk_time;
        parlay::par_do(
            [&]()
            {
                const auto t_0 = now();
                const auto diag_blk_sz = G.E().read_diagonal_block(i, buf);
                const auto t_1 = now();
                edge_read_time += duration(t_1 - t_0);

                parlay::par_do(
                    [&]()
                    {
                        G.E().remove_block(i, i);
                    },
                    [&]()
                    {
                        const auto t_2 = now();
                        parlay::parallel_for(0, diag_blk_sz,
                            [&](const std::size_t idx)
                            {
                                const auto& e = buf[idx];
                                if(e.w() == 1)
                                {
                                    assert(M.find(e.x()) && M.find(e.y()));
                                    add_diagonal_edge_path_info(e, *M.find(e.x()));
                                }
                            });

                        const auto t_3 = now();
                        diag_blk_time = duration(t_3 - t_2);
                    }
                );
            },
            [&]()
            {
                const auto t_0 = now();
                const auto top_blk_sz = G.E().read_block(0, i, swap_buf);
                const auto t_1 = now();
                edge_read_time += duration(t_1 - t_0);

                parlay::par_do(
                    [&]()
                    {
                        G.E().remove_block(0, i);
                    },
                    [&]()
                    {
                        const auto t_2 = now();
                        parlay::parallel_for(0, top_blk_sz,
                            [&](const std::size_t idx)
                            {
                                const auto& e = swap_buf[idx];
                                if(e.w() == 1)
                                {
                                    assert(M.find(e.y()));
                                    add_edge_path_info(e, *M.find(e.y()));
                                }
                            });

                        const auto t_3 = now();
                        top_blk_time = duration(t_3 - t_2);
                    }
                );
            }
        );

        spec_case_time += std::max(diag_blk_time, top_blk_time);
    }

    std::cerr << "\n";


    // std::cerr << "Encountered " << og_edge_c << " original edges\n";
    std::cerr << "Map clearing time: " << map_clr_time << ".\n";
    std::cerr << "Edges reading time: " << edge_read_time << ".\n";
    std::cerr << "Vertex-info load time: " << p_v_load_time << ".\n";
    std::cerr << "Map filling time: " << map_fill_time << ".\n";
    std::cerr << "Non-diagonal blocks expansion time: " << edge_proc_time << ".\n";
    std::cerr << "Diagonal block expansion time: " << diag_exp_time << ".\n";
    std::cerr << "Special case time: " << spec_case_time << ".\n";
}


template <uint16_t k, bool Colored_>
std::size_t Contracted_Graph_Expander<k, Colored_>::load_path_info(const std::size_t i, Buffer<Obj_Path_Info_Pair<Kmer<k>, k>>& p_v_buf)
{
    const auto t_s = now();
    const auto sz = P_v[i].unwrap().size();
    p_v_buf.reserve_uninit(sz);
    P_v[i].unwrap().load(p_v_buf.data());
    const auto t_e = now();
    p_v_load_time += duration(t_e - t_s);

    return sz;
}


template <uint16_t k, bool Colored_>
void Contracted_Graph_Expander<k, Colored_>::fill_path_info(const Buffer<Obj_Path_Info_Pair<Kmer<k>, k>>& p_v_buf, const std::size_t buf_sz)
{
    const auto t_s = now();
    const auto load_vertex_info = [&](const std::size_t idx)
    {
        const auto& p_v = p_v_buf[idx];
        if(!M.insert(p_v.obj(), p_v.path_info()))
            assert(*M.find(p_v.obj()) == p_v.path_info());
    };

    parlay::parallel_for(0, buf_sz, load_vertex_info);

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            for(const auto& p_v : P_v_ip1[w_id].unwrap())
                if(!M.insert(p_v.obj(), p_v.path_info()))
                    assert(*M.find(p_v.obj()) == p_v.path_info());
        },
        1);

    const auto t_e = now();
    map_fill_time += duration(t_e - t_s);
}


template <uint16_t k, bool Colored_>
void Contracted_Graph_Expander<k, Colored_>::expand_diagonal_block(const std::size_t i)
{
    const std::string d_i_path(compressed_diagonal_path + "_" + std::to_string(i));
    const auto file_sz = file_size(d_i_path);
    const auto edge_c = file_sz / sizeof(Discontinuity_Edge<k>);
    D_i.reserve_uninit(edge_c);

    const auto t_s = now();
    std::ifstream input(d_i_path);
    input.read(reinterpret_cast<char*>(D_i.data()), file_sz);
    assert(static_cast<std::size_t>(input.gcount()) == file_sz);
    input.close();
    const auto t_e = now();
    edge_read_time += duration(t_e - t_s);


    parlay::par_do(
        [&]()
        {
            if(!input || !remove_file(d_i_path))
            {
                std::cerr << "Error reading / removing of contracted edge block from " << d_i_path << ". Aborting.\n";
                std::exit(EXIT_FAILURE);
            }
        },
        [&]()
        {
            Path_Info<k> x_inf, y_inf;
            for(int64_t idx = static_cast<int64_t>(edge_c) - 1; idx >= 0; --idx)    // In reverse order of the newly introduced diagonal-edges to always ensure one endpoint having path-info ready.
            {
                const auto& e = D_i[idx];
                assert(M.find(e.x()) || M.find(e.y()));

                if(!M.find(e.y(), y_inf))
                {
                    assert(M.find(e.x()));
                    x_inf = *M.find(e.x()); // M.find(e.x(), x_inf);
                    y_inf = infer(x_inf, e.s_x(), e.s_y(), e.w());
                    if(y_inf.r() > 0)
                        M.insert(e.y(), y_inf);
                }
                else
                {
                    assert(M.find(e.y()));
                    x_inf = infer(y_inf, e.s_y(), e.s_x(), e.w());
                    if(x_inf.r() > 0)
                        M.insert(e.x(), x_inf);
                }
            }
        }
    );
}


/*
template <uint16_t k>
template <typename T_s_, typename T_d_>
uint64_t Contracted_Graph_Expander<k>::collate_w_local_bufs(T_s_& source, const std::size_t beg, const size_t end, T_d_& dest)
{
    std::vector<Padded_Data<uint64_t>> C(parlay::num_workers(), 0);

    parlay::parallel_for(beg, end,
        [&](const std::size_t idx)
        {
            for(auto& w_local_buf_list : source)
            {
                auto& buf = w_local_buf_list.data()[idx];
                dest[idx].add(buf.data(), buf.size());
                C[parlay::worker_id()].data() += buf.size();

                buf.clear();
            }
        }, 1);

    uint64_t c = 0;
    std::for_each(C.cbegin(), C.cend(), [&](auto& v){ c += v.data(); });
    return c;
}
*/

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Contracted_Graph_Expander)
