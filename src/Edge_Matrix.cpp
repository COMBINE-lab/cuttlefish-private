
#include "Edge_Matrix.hpp"
#include "globals.hpp"

#include <cstddef>


namespace cuttlefish
{

template <uint16_t k>
Edge_Matrix<k>::Edge_Matrix(std::size_t part_count, const std::string& path):
      vertex_part_count(part_count)
    , path(path)
{
    // TODO: fix better policy.
    if((part_count & (part_count - 1)) != 0)
    {
        std::cerr << "Vertex-partition count needs to be a power of 2. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    edge_matrix.resize(part_count + 1);
    for(std::size_t i = 0; i <= part_count; ++i)
        for(std::size_t j = 0; j <= part_count; ++j)
            j < i ? edge_matrix[i].emplace_back() :
                    edge_matrix[i].emplace_back(path + std::string(".") + std::to_string(i) + "-" + std::to_string(j) + std::string(".blk"));
}


template <uint16_t k>
void Edge_Matrix<k>::serialize()
{
    for(std::size_t i = 0; i <= vertex_part_count; ++i)
        for(std::size_t j = 0; j <= vertex_part_count; ++j)
            edge_matrix[i][j].close();
}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Edge_Matrix)
