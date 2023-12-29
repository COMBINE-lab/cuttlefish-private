
#include "Discontinuity_Graph.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{

template <uint16_t k> const Kmer<k> Discontinuity_Graph<k>::phi_(phi_label);


template <uint16_t k>
Discontinuity_Graph<k>::Discontinuity_Graph(const std::size_t part_count, const std::size_t lmtig_bucket_count, const std::string& work_path):
      E_(part_count, work_path + "E_")
    , lmtigs(work_path + "lmtig", lmtig_bucket_count, parlay::num_workers())
    , phantom_edge_count_(0)
{}


template <uint16_t k>
uint64_t Discontinuity_Graph<k>::phantom_edge_upper_bound() const
{
    return phantom_edge_count_;
}


template <uint16_t k>
void Discontinuity_Graph<k>::close_lmtig_stream()
{
    lmtigs.close();
}


template <uint16_t k>
std::size_t Discontinuity_Graph<k>::vertex_part_size_upper_bound() const
{
    std::size_t bound = 0;
    for(std::size_t j = 1; j <= E_.vertex_part_count(); ++j)
        // Each *original* edge from the non-diagonal blocks of column `j` corresponds to a
        // unique vertex of partition `j`, and in the worst-case, each original edge from the
        // diagonal block corresponds to two unique vertices.
        bound = std::max(bound, E_.col_size(j) + E_.block_size(j, j));

    return bound;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Discontinuity_Graph)
