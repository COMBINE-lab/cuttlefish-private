
#include "Contracted_Graph_Expander.hpp"
#include "globals.hpp"

#include <fstream>
#include <filesystem>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k>
Contracted_Graph_Expander<k>::Contracted_Graph_Expander(Edge_Matrix<k>& E, std::vector<Ext_Mem_Bucket<Vertex_Path_Info_Pair<k>>>& P_v, const std::string& temp_path):
      E(E)
    , P_v(P_v)
    , work_path(temp_path)
{}


template <uint16_t k>
void Contracted_Graph_Expander<k>::expand_diagonal_block(const std::size_t i)
{
    const std::string d_i_path(work_path + std::string("D_") + std::to_string(i));
    std::error_code ec;
    const auto file_sz = std::filesystem::file_size(d_i_path, ec);
    D_i.resize(file_sz / sizeof(Discontinuity_Edge<k>));

    std::ifstream input(d_i_path);
    input.read(reinterpret_cast<char*>(D_i.data()), file_sz);
    input.close();

    // In reverse order of the newly introduced diagonal-edges to always ensure one endpoint having path-info ready.
    for(auto d_it = D_i.rbegin(); d_it != D_i.rend(); ++d_it)
    {
        const auto& e = *d_it;

        assert(M.find(e.u()) != M.end() || M.find(e.v()) != M.end);

        auto it = M.find(e.v());
        if(it == M.end())
            M.emplace(e.v(), infer(M.find(e.u())->second, e.s_u(), e.s_v(), e.w()));
        else
            M.emplace(e.u(), infer(it->second, e.s_v(), e.s_u(), e.w()));


        if(e.w() == 1)
        {
            // TODO: add edge path-info appropriately.
        }
    }
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Contracted_Graph_Expander)
