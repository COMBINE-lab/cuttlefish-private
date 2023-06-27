
#ifndef EDGE_MATRIX_HPP
#define EDGE_MATRIX_HPP



#include "Discontinuity_Edge.hpp"
#include "Ext_Mem_Bucket.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>


namespace cuttlefish
{

// =============================================================================
// Blocked edge-matrix of a discontinuity-graph of `k`-mers.
template <uint16_t k>
class Edge_Matrix
{
private:

    static const std::string edge_block_ext;

    const std::size_t vertex_part_count_;   // Number of vertex-partitions in the graph; it needs to be a power of 2.
    const std::string path; // File-path prefix to the external-memory blocks of the matrix.
    std::vector<std::vector<Ext_Mem_Bucket<Discontinuity_Edge<k>>>> edge_matrix;    // Blocked edge matrix.

    mutable std::vector<std::size_t> row_to_read;   // `j`'th entry contains the row of the next block to read from column `j`.
    mutable std::vector<std::size_t> col_to_read;   // `i`'th entry contains the column of the next block to read from row `i`.


    // Returns the path to the file storing the `[i, j]`'th block in external
    // memory.
    const std::string bucket_file_path(std::size_t i, std::size_t j) const;


public:

    // Constructs a blocked edge-matrix for `part_count` vertex-partitions. The
    // partition-count needs to be a power of 2.
    Edge_Matrix(std::size_t part_count, const std::string& path);

    // Returns the number of vertex-partitions in the graph.
    std::size_t vertex_part_count() const { return vertex_part_count_; }

    // Returns the partition ID for the k-mer `kmer`.
    std::size_t partition(const Kmer<k>& kmer) const;

    // Adds a discontinuity-edge `e` to the matrix.
    void add(Discontinuity_Edge<k> e);

    // Adds the discontinuity-edge `({(u, s_u), (v, s_v)}, w, b)` to the matrix.
    // `u_is_phi` and `v_is_phi` denote whether the `u` and the `v` endpoints
    // are ϕ, respectively.
    void add(Kmer<k> u, side_t s_u, Kmer<k> v, side_t s_v, uint16_t w, uint16_t b, bool u_is_phi = false, bool v_is_phi = false);

    // Serializes the matrix to external-memory. Edges should not be added
    // anymore once this has been invoked. This method is required only if the
    // entirety of the edge-matrix needs to live in external-memory after the
    // parent process finishes.
    void serialize();

    // Reads the edges from the `[j, j]`'th block into `buf`.
    void read_diagonal_block(std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const;

    // Reads a chunk of edges from the column `j` into `buf`. Returns `true` iff
    // some edges are read, i.e. the column had remaining edges to be read off.
    // NB: this does not read the blocks in the diagonal.
    bool read_column_buffered(std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const;

    // Reads a chunk of edges from the row `i` into `buf`. Returns `true` iff
    // some edges are read, i.e. the row had remaining edges to be read off.
    // NB: this does not read the blocks in the diagonal.
    bool read_row_buffered(std::size_t i, std::vector<Discontinuity_Edge<k>>& buf) const;

    // Reads the edges from the `[i, j]`'th block into `buf`.
    void read_block(std::size_t i, std::size_t j, std::vector<Discontinuity_Edge<k>>& buf) const;
};


template <uint16_t k>
inline std::size_t Edge_Matrix<k>::partition(const Kmer<k>& kmer) const
{
    return (kmer.to_u64() & (vertex_part_count_ - 1)) + 1;
}


template <uint16_t k>
inline void Edge_Matrix<k>::add(Discontinuity_Edge<k> e)
{
    auto p = e.u_is_phi() ? 0 : partition(e.u());
    auto q = e.v_is_phi() ? 0 : partition(e.v());

    if(p > q)
        e.invert(),
        std::swap(p, q);

    edge_matrix[p][q].add(e);
}


template <uint16_t k>
inline void Edge_Matrix<k>::add(const Kmer<k> u, const side_t s_u, const Kmer<k> v, const side_t s_v, const uint16_t w, const uint16_t b, const bool u_is_phi, const bool v_is_phi)
{
    auto p = u_is_phi ? 0 : partition(u);
    auto q = v_is_phi ? 0 : partition(v);

    if(p <= q)
        edge_matrix[p][q].emplace(u, s_u, v, s_v, w, b, u_is_phi, v_is_phi);
    else    // p and q needs to be swapped along with the edge endpoints
        edge_matrix[q][p].emplace(v, s_v, u, s_u, w, b, v_is_phi, u_is_phi);
}

}



#endif
