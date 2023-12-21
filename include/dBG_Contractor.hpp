
#ifndef DBG_CONTRACTOR_HPP
#define DBG_CONTRACTOR_HPP



#include "Discontinuity_Graph.hpp"
#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Output_Sink.hpp"
#include "Character_Buffer.hpp"
#include "utility.hpp"
#include "globals.hpp"

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

    // TODO: wrap the following primitive fields in `Build_Params`.
    const std::size_t subgraph_count;   // Number of subgraphs the original de Bruijn graph is broken into.
    const std::size_t vertex_part_count;    // Number of vertex-partitions in the discontinuity graph; needs to be a power of 2.
    const std::size_t lmtig_bucket_count;   // Number of buckets storing literal locally-maximal unitigs.

    const std::string output_path;  // Path-prefix to output files.
    const std::string work_path;    // Path-prefix for temporary working files.

    Discontinuity_Graph<k> G;   // The discontinuity graph.

    std::size_t n_disc_v;   // Number of discontinuity-vertices.

    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<Kmer<k>, k>>> P_v;    // `P_v[j]` contains path-info for vertices in partition `j`.
    std::vector<Ext_Mem_Bucket<Obj_Path_Info_Pair<uni_idx_t, k>>> P_e;  // `P_e[b]` contains path-info for edges induced by unitigs in bucket `b`.

    typedef Async_Logger_Wrapper sink_t;
    Output_Sink<sink_t> output_sink;    // Sink for the output maximal unitigs.

    // 100 KB (soft limit) worth of maximal unitig records (FASTA) can be retained in memory per worker, at most, before flushes.
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded_Data<op_buf_t>> op_buf_list_t;
    op_buf_list_t op_buf;   // Worker-specific output buffers.


public:

    dBG_Contractor(const dBG_Contractor&) = delete;
    dBG_Contractor& operator=(const dBG_Contractor&) = delete;

    // Constructs a compacted de Bruijn graph constructor that operates with
    // `subgraph_count` subgraphs of the underlying de Bruijn graph,
    // `part_count` vertex-partitions in its discontinuity graph, and stores the
    // locally-maximal unitigs from the partitioned subgraphs of the dBG into
    // `lmtig_bucket_count` buckets. Output files are stored at the path-prefix
    // `output_path`, and temporary working files at path-prefix `temp_path`.
    dBG_Contractor(std::size_t subgraph_count, std::size_t part_count, std::size_t lmtig_bucket_count, const std::string& output_path, const std::string& temp_path);

    // Contracts the bootstrapped discontinuity graph generated from the
    // compacted dBG at path `cdbg_path`. `l`-minimizers are used in generating
    // the discontinuity graph.
    // TODO: parameters were there for bootstrapping purposes. Useless now.
    void contract(uint16_t l, const std::string& cdbg_path);
};

}



#endif
