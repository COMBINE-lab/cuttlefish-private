
#ifndef UNITIG_COLLATOR_HPP
#define UNITIG_COLLATOR_HPP



#include "Ext_Mem_Bucket.hpp"
#include "Path_Info.hpp"

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>


namespace cuttlefish
{

// =============================================================================
// Collates locally-maximal unitigs from different unitig-buckets as per their
// path-information in a discontinuity graph of `k`-mers.
template <uint16_t k>
class Unitig_Collator
{
    typedef uint32_t uni_idx_t; // Type of the index of a unitig in a bucket.

private:

    const std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e;   // `P_e[b]` contains path-info for edges in bucket `b`.

    const std::string output_path;  // Path-prefix to output paths.
    const std::string work_path;    // Path-prefix to temporary working files.

    std::vector<Path_Info<k>> M;    // `M[idx]` is the path-info for the `idx`'th edge in some bucket.

    std::vector<Obj_Path_Info_Pair<uni_idx_t, k>> p_e_buf;  // Buffer to read-in path-information of edges.


    // Loads the path-info of edges from bucket `b` into the table `M`.
    void load_path_info(std::size_t b);


public:

    // Constructs a unitig-collator for unitigs with their associated path-info
    // at `P_e`, i.e. `P_e[b]` contains path-information of the unitigs'
    // corresponding edges at bucket `b`. Output files and temporary files are
    // stored at path-prefixes `output_path` and`temp_path`, respectively.
    Unitig_Collator(const std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const std::string& output_path, const std::string& temp_path);

    void collate();
};

}



#endif
