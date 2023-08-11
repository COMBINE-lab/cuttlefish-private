
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
    typedef uint32_t uni_idx_t; // Type of the index of a unitig in a bucket.

private:

    // TODO: wrap the following primitive fields in `Build_Params`.
    const std::size_t part_count;   // Number of vertex-partitions in the discontinuity graph; needs to be a power of 2.
    const std::size_t unitig_bucket_count;  // Number of buckets storing literal unitigs.

    const std::string output_path;  // Path-prefix to output files.
    const std::string work_path;    // Path-prefix for temporary working files.

    Edge_Matrix<k> E;  // Edge-matrix of the discontinuity graph.

    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>> P_v;    // `P_v[j]` contains path-info for vertices in partition `j`.
    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>> P_e;  // `P_e[b]` contains path-info for edges induced by unitigs in bucket `b`.


public:

    dBG_Contractor(const dBG_Contractor&) = delete;
    dBG_Contractor& operator=(const dBG_Contractor&) = delete;

    // Constructs a compacted de Bruijn graph constructor that operates with
    // `part_count` vertex-partitions in its discontinuity graph, and stores the
    // locally-maximal unitigs from the partitioned subgraphs of the dBG into
    // `unitig_bucket_count` buckets. Output files are stored at the path-prefix
    // `output_path`, and temporary working files at path-prefix `temp_path`.
    dBG_Contractor(std::size_t part_count, std::size_t unitig_bucket_count, const std::string& output_path, const std::string& temp_path);

    // Contracts the bootstrapped discontinuity graph generated from the
    // compacted dBG at path `cdbg_path`. `l`-minimizers are used in generating
    // the discontinuity graph.
    void contract(uint16_t l, const std::string& cdbg_path);
};

}



#endif
