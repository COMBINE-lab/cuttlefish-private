
#ifndef UNITIG_COLLATOR_HPP
#define UNITIG_COLLATOR_HPP



#include "Ext_Mem_Bucket.hpp"
#include "Path_Info.hpp"
#include "Unitig_Coord_Bucket.hpp"
#include "Output_Sink.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "globals.hpp"

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
private:

    typedef Obj_Path_Info_Pair<uni_idx_t, k> unitig_path_info_t;    // A locally-maximal unitig's index in its bucket and its path-info.

    const std::vector<Ext_Mem_Bucket<unitig_path_info_t>>& P_e; // `P_e[b]` contains path-info for edges in bucket `b`.

    const std::string output_path;  // Path-prefix to output paths.
    const std::string work_path;    // Path-prefix to temporary working files.

    std::size_t max_bucket_sz;  // Maximum size of the locally-maximal unitigs' buckets.

    static constexpr std::size_t max_unitig_bucket_count = 1024;    // Must be a power-of-2.
    std::vector<Unitig_Coord_Bucket<k>> max_unitig_bucket;  // Key-value collation buckets for lm-unitigs.

    // TODO: remove? This is for the naive-collator.
    Path_Info<k>* M;    // `M[idx]` is the path-info for the `idx`'th edge in some bucket.

    // TODO: remove?  This is for the naive-collator.
    unitig_path_info_t* p_e_buf;    // Buffer to read-in path-information of edges.

    typedef Async_Logger_Wrapper sink_t;
    Output_Sink<sink_t> output_sink;    // Sink for the output maximal unitigs.


    // Maps each locally-maximal unitig to its maximal unitig's corresponding
    // bucket.
    void map();

    // Reduces each maximal unitig bucket to its contained maximal unitigs.
    void reduce();

    // Loads the path-info of edges from bucket `b` into the table `M`, and
    // returns the size of the bucket.
    std::size_t load_path_info(std::size_t b);

    // Loads the path-info of edges from bucket `b` into the table `M`, and
    // returns the size of the bucket. Uses the buffer `buf` to transfer the
    // information from the bucket to the table.
    std::size_t load_path_info(std::size_t b, Path_Info<k>* M, unitig_path_info_t* buf);


public:

    // Constructs a unitig-collator for unitigs with their associated path-info
    // at `P_e`, i.e. `P_e[b]` contains path-information of the unitigs'
    // corresponding edges at bucket `b`. Output files and temporary files are
    // stored at path-prefixes `output_path` and`temp_path`, respectively.
    Unitig_Collator(const std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const std::string& output_path, const std::string& temp_path);

    // Collates the locally-maximal unitigs into global ones.
    void collate();

    // Collates the locally-maximal unitigs into global ones.
    void par_collate();
};

}



#endif
