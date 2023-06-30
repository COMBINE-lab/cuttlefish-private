
#ifndef DBG_CONTRACTOR_HPP
#define DBG_CONTRACTOR_HPP



#include "Edge_Matrix.hpp"
#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>


namespace cuttlefish
{

// =============================================================================
// Compacted de Bruijn graph constructor.
template <uint16_t k>
class dBG_Contractor
{
private:

    const std::size_t part_count;   // Number of vertex-partitions in the discontinuity graph; needs to be a power of 2.
    const std::string work_path;    // Path-prefix for temporary working files.

    Edge_Matrix<k> E;  // Edge-matrix of the discontinuity graph.

    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>> P_v;    // `P_v[j]` contains path-info for vertices in partition `j`.


public:

    dBG_Contractor(const dBG_Contractor&) = delete;

    // Constructs a compacted de Bruijn constructor that operates with
    // `part_count` vertex-partitions in its discontinuity graph and stores
    // temporary working files at path-prefix `temp_path`.
    dBG_Contractor(std::size_t part_count, const std::string& temp_path);

    // Contracts the bootstrapped discontinuity graph generated from the
    // compacted dBG at path `cdbg_path`. `l`-minimizers are used in generating
    // the discontinuity graph.
    void contract(uint16_t l, const std::string& cdbg_path);
};

}



#endif
