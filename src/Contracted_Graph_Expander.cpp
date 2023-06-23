
#include "Contracted_Graph_Expander.hpp"


namespace cuttlefish
{

template <uint16_t k>
Contracted_Graph_Expander<k>::Contracted_Graph_Expander(Edge_Matrix<k>& E, std::vector<Ext_Mem_Bucket<Vertex_Path_Info_Pair<k>>>& P_v, const std::string& temp_path):
      E(E)
    , P_v(P_v)
    , work_path(temp_path)
{}

}
