
#include "Edge_Matrix.hpp"
#include "globals.hpp"

#include <cstddef>
#include <algorithm>


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
    lock.resize(part_count + 1);
    for(std::size_t i = 0; i <= part_count; ++i)
    {
        for(std::size_t j = 0; j <= part_count; ++j)
            j < i ? edge_matrix[i].emplace_back() :
                    edge_matrix[i].emplace_back(bucket_file_path(i, j));

        lock[i] = new Spin_Lock[part_count + 1];
    }

    for(std::size_t i = 1; i <= part_count; ++i)
        col_to_read[i] = i + 1;
}


template <uint16_t k>
Edge_Matrix<k>::~Edge_Matrix()
{
    for(auto lck : lock)
        delete[] lck;
}


template <uint16_t k>
const std::string Edge_Matrix<k>::bucket_file_path(const std::size_t i, const std::size_t j) const
{
    return path +  "_" + std::to_string(i) + "-" + std::to_string(j) + edge_block_ext;
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


template <uint16_t k>
std::size_t Edge_Matrix<k>::row_size(const std::size_t i) const
{
    assert(i <= vertex_part_count_);

    std::size_t sz = 0;
    for(std::size_t j = std::max(1lu, i); j <= vertex_part_count_; ++j)
        sz += edge_matrix[i][j].size();

    return sz;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::col_size(const std::size_t j) const
{
    assert(j <= vertex_part_count_);

    std::size_t sz = 0;
    for(std::size_t i = 0; i <= j; ++i)
        sz += edge_matrix[i][j].size();

    return sz;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::block_size(const std::size_t i, const std::size_t j) const
{
    return edge_matrix[i][j].size();
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::size() const
{
    std::size_t sz = 0;
    for(std::size_t i = 0; i <= vertex_part_count_; ++i)
        sz += row_size(i);

    return sz;
}


template <uint16_t k>
std::size_t Edge_Matrix<k>::max_block_size() const
{
    std::size_t max_sz = 0;
    for(std::size_t i = 0; i <= vertex_part_count_; ++i)
        for(std::size_t j = std::max(1lu, i); j <= vertex_part_count_; ++j)
            max_sz = std::max(max_sz, edge_matrix[i][j].size());

    return max_sz;
}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Edge_Matrix)
