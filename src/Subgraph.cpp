
#include "Subgraph.hpp"
#include "globals.hpp"
#include "kmc-super-kmers-iterator/iterate_super_kmers.h"



namespace cuttlefish
{

template <uint16_t k>
Subgraph<k>::Subgraph(const std::string& bin_dir_path, const std::size_t bin_id):
      graph_bin_dir_path(bin_dir_path)
    , bin_id(bin_id)
{}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraph)
