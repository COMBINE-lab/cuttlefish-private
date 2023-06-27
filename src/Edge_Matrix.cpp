
#include "Edge_Matrix.hpp"
#include "globals.hpp"

#include <cstddef>


namespace cuttlefish
{

template <uint16_t k> const std::string Edge_Matrix<k>::edge_block_ext(".blk");

template <uint16_t k>
Edge_Matrix<k>::Edge_Matrix(std::size_t part_count, const std::string& path):
      vertex_part_count_(part_count)
    , path(path)
    , row_to_read(part_count + 1, 0)
    , col_to_read(part_count + 1)
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
                    edge_matrix[i].emplace_back(bucket_file_path(i, j));

    for(std::size_t i = 1; i <= part_count; ++i)
        col_to_read[i] = i + 1;
}


template <uint16_t k>
const std::string Edge_Matrix<k>::bucket_file_path(const std::size_t i, const std::size_t j) const
{
    return path + std::to_string(i) + "-" + std::to_string(j) + edge_block_ext;
}


template <uint16_t k>
void Edge_Matrix<k>::serialize()
{
    for(std::size_t i = 0; i <= vertex_part_count_; ++i)
        for(std::size_t j = 0; j <= vertex_part_count_; ++j)
            edge_matrix[i][j].close();
}


template <uint16_t k>
void Edge_Matrix<k>::read_diagonal_block(const std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const
{
    edge_matrix[j][j].load(buf);
}


template <uint16_t k>
bool Edge_Matrix<k>::read_column_buffered(const std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const
{
    if(row_to_read[j] >= j)
        return false;

    edge_matrix[row_to_read[j]][j].load(buf);
    row_to_read[j]++;

    return true;
}


template <uint16_t k>
bool Edge_Matrix<k>::read_row_buffered(const std::size_t i, std::vector<Discontinuity_Edge<k>>& buf) const
{
    if(col_to_read[i] > vertex_part_count_)
        return false;

    edge_matrix[i][col_to_read[i]].load(buf);
    col_to_read[i]++;

    return true;
}


template <uint16_t k>
void Edge_Matrix<k>::read_block(std::size_t i, std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const
{
    assert(i <= vertex_part_count_ && j <= vertex_part_count_);
    edge_matrix[i][j].load(buf);
}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Edge_Matrix)
