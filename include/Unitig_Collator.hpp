
#ifndef UNITIG_COLLATOR_HPP
#define UNITIG_COLLATOR_HPP



#include "dBG_Contractor.hpp"
#include "Path_Info.hpp"
#include "Unitig_Coord_Bucket.hpp"
#include "Output_Sink.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Character_Buffer.hpp"
#include "utility"
#include "globals.hpp"

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Collates locally-maximal unitigs from different unitig-buckets as per their
// path-information in a discontinuity graph of `k`-mers.
template <uint16_t k>
class Unitig_Collator
{
private:

    typedef typename dBG_Contractor<k>::unitig_path_info_t unitig_path_info_t;

    typedef typename dBG_Contractor<k>::P_e_t P_e_t;
    const P_e_t& P_e;   // `P_e[b]` contains path-info for edges in bucket `b`.

    const std::string lmtig_buckets_path;   // Path-prefix to the lm-tig buckets.
    const std::string unitig_coord_buckets_path;    // Path-prefix to the unitig-coordinate buckets produced in map-reduce.

    std::size_t max_bucket_sz;  // Maximum size of the locally-maximal unitigs' buckets.

    static constexpr std::size_t max_unitig_bucket_count = 1024;    // Must be a power-of-2.
    std::vector<Padded_Data<Unitig_Coord_Bucket<k>>> max_unitig_bucket; // Key-value collation buckets for lm-unitigs.

    // TODO: remove? This is for the naive-collator.
    Path_Info<k>* M;    // `M[idx]` is the path-info for the `idx`'th edge in some bucket.

    // TODO: remove?  This is for the naive-collator.
    unitig_path_info_t* p_e_buf;    // Buffer to read-in path-information of edges.

    // TODO: move the following to some CF3-centralized location.

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded_Data<op_buf_t>> op_buf_list_t;
    op_buf_list_t& op_buf;  // Worker-specific output buffers.


    // Maps each locally-maximal unitig to its maximal unitig's corresponding
    // bucket.
    void map();

    // Reduces each maximal unitig bucket to its contained maximal unitigs.
    void reduce();

    // Loads the path-info of edges from bucket `b` into the table `M`, and
    // returns the size of the bucket. Uses the buffer `buf` to transfer the
    // information from the bucket to the table.
    std::size_t load_path_info(std::size_t b, Path_Info<k>* M, Buffer<unitig_path_info_t>& buf);

    // Returns the index of the lexicographically minimum k-mer in the `s`, and
    // puts the k-mer in `min`.
    static std::size_t min_kmer(const std::string& s, Kmer<k>& min);


public:

    // Constructs a unitig-collator for unitigs with their associated path-info
    // at `P_e`, i.e. `P_e[b]` contains path-information of the unitigs'
    // corresponding edges at bucket `b`. `logistics` is the data logistics
    // manager for the algorithm execution. Worker-specific maximal unitigs are
    // written to the buffers in `op_buf`.
    Unitig_Collator(const P_e_t& P_e, const Data_Logistics& logistics, op_buf_list_t& op_buf);

    // Collates the locally-maximal unitigs into global ones.
    // void collate();

    // Collates the locally-maximal unitigs into global ones.
    void par_collate();
};

}



#endif
