
#include "Contracted_Graph_Expander.hpp"


namespace cuttlefish
{

template <uint16_t k>
Contracted_Graph_Expander<k>::Contracted_Graph_Expander(Edge_Matrix<k>& E, const std::string& temp_path):
      E(E)
    , work_path(temp_path)
{}

}
