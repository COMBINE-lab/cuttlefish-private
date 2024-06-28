
#ifndef UNITIG_COLLATOR_HPP
#define UNITIG_COLLATOR_HPP



#include "dBG_Contractor.hpp"
#include "Path_Info.hpp"
#include "Unitig_Coord_Bucket.hpp"
#include "Character_Buffer.hpp"
#include "utility"
#include "globals.hpp"

#include <cstdint>
#include <vector>
#include <string>


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

    std::size_t max_bucket_sz;  // Maximum size of the edge-buckets.

    const std::size_t max_unitig_bucket_count;  // Number of buckets storing literal globally-maximal unitigs.
    std::vector<Padded_Data<Unitig_Coord_Bucket<k>>> max_unitig_bucket; // Key-value collation buckets for lm-unitigs.

    typedef typename dBG_Contractor<k>::op_buf_list_t op_buf_list_t;
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
    // written to the buffers in `op_buf`. `gmtig_bucket_count` many buckets
    // are used to partition the lm-tigs to their maximal unitigs.
    Unitig_Collator(const P_e_t& P_e, const Data_Logistics& logistics, op_buf_list_t& op_buf, std::size_t gmtig_bucket_count);

    // Collates the locally-maximal unitigs into global ones.
    void collate();
};

}



#endif
