
#ifndef DBG_CONTRACTOR_HPP
#define DBG_CONTRACTOR_HPP



#include "Discontinuity_Graph.hpp"
#include "Path_Info.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Output_Sink.hpp"
#include "Character_Buffer.hpp"
#include "Build_Params.hpp"
#include "Data_Logistics.hpp"
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
public:

    typedef Obj_Path_Info_Pair<Kmer<k>, k> vertex_path_info_t;  // A vertex and its path-info.
    typedef Obj_Path_Info_Pair<uni_idx_t, k> unitig_path_info_t;    // A locally-maximal unitig's index in its bucket and its path-info.

    typedef std::vector<Ext_Mem_Bucket_Concurrent<vertex_path_info_t>> P_v_t;
    typedef std::vector<Ext_Mem_Bucket_Concurrent<unitig_path_info_t>> P_e_t;

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded_Data<op_buf_t>> op_buf_list_t;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    const Data_Logistics logistics; // Data logistics manager for the algorithm execution.

    Discontinuity_Graph<k> G;   // The discontinuity graph.

    std::size_t n_disc_v;   // Number of discontinuity-vertices.

    // TODO: consider using padding.
    P_v_t P_v;  // `P_v[j]` contains path-info for vertices in partition `j`.
    P_e_t P_e;  // `P_e[b]` contains path-info for edges induced by unitigs in bucket `b`.

    Output_Sink<sink_t> output_sink;    // Sink for the output maximal unitigs.

    // 100 KB (soft limit) worth of maximal unitig records (FASTA) can be retained in memory per worker, at most, before flushes.
    op_buf_list_t op_buf;   // Worker-specific output buffers.


public:

    dBG_Contractor(const dBG_Contractor&) = delete;
    dBG_Contractor& operator=(const dBG_Contractor&) = delete;

    // Constructs a compacted de Bruijn graph constructor with the parameters
    // required for the construction wrapped in `params`.
    dBG_Contractor(const Build_Params& params);

    // Contracts the bootstrapped discontinuity graph generated from the
    // compacted dBG at path `cdbg_path`. `l`-minimizers are used in generating
    // the discontinuity graph.
    void construct();
};

}



#endif
