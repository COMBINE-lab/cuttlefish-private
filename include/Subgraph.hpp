
#ifndef SUBGRAPH_HPP
#define SUBGRAPH_HPP



#include "State_Config.hpp"
#include "Directed_Vertex.hpp"
#include "Kmer.hpp"
#include "DNA.hpp"
#include "DNA_Utility.hpp"
#include "Unitig_Scratch.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "Discontinuity_Graph.hpp"
#include "dBG_Utilities.hpp"
#include "Character_Buffer.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "utility.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <unordered_map>


namespace cuttlefish
{

// =============================================================================
// A subgraph of the de Bruijn graph, induced by a KMC bin.
template <uint16_t k>
class Subgraph
{
private:

    enum class Walk_Termination;    // Type of scenarios how a unitig-walk terminates in the subgraph.
    typedef Walk_Termination termination_t;

    const std::string graph_bin_dir_path;   // Path to the directory with all the graph KMC-bins.
    const std::size_t bin_id; // ID of the graph KMC-bin.

    typedef std::unordered_map<Kmer<k>, State_Config, Kmer_Hasher<k>> map_t;
    static std::vector<Padded_Data<map_t>> map; // Map collection for different workers.

    map_t& M;   // Map to be used for this subgraph.

    uint64_t edge_c;    // Number of edges in the graph.
    uint64_t label_sz;  // Total number of characters in the literal representations of all the maximal unitigs.
    uint64_t disc_edge_c;   // Number of edges of the discontinuity graph induced from this subgraph.
    uint64_t isolated;  // Count of isolated vertices—not part of any edge.

    Discontinuity_Graph<k>& G;  // The discontinuity graph.

    uint64_t trivial_mtig_c;    // Number of trivial maximal unitigs in the graph (i.e. also maximal unitigs in the supergraph).
    uint64_t icc_count_;    // Number of trivial maximal unitigs in the graph that are ICCs.


    // TODO: move the following out to a central location.

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    op_buf_t& op_buf;   // Output buffer for trivially maximal unitigs of the underlying dBG.


    // Extracts the maximal unitig containing the vertex `v_hat`, and
    // `maximal_unitig` is used as the working scratch for the extraction, i.e.
    // to build and store two unitigs connecting to the two sides of `v_hat`.
    // Returns `true` iff the containing maximal unitig has not been outputted
    // earlier.
    bool extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig);

    // Traverses a unitig starting from the vertex `v_hat`, exiting it through
    // the side `s_v_hat`. `st_v` is the state of `s_v_hat`. `unitig` is used as
    // the scratch space to build the unitig. Returns `true` iff the walk tried
    // to exit the subgraph through a discontinuous side; in which case that
    // vertex is stored in `exit_v`.
    termination_t walk_unitig(const Kmer<k>& v_hat, const State_Config& v_inf, side_t s_v_hat, Unitig_Scratch<k>& unitig, Directed_Vertex<k>& exit_v);


public:

    // Constructs a subgraph object for the `bin_id`'th bin in the graph bin
    // directory `bin_dir_path`. Updates the discontinuity graph `G` with its
    // edges observed from this subgraph and writes the trivially maximal
    // unitigs to `op_buf`.
    Subgraph(const std::string& bin_dir_path, std::size_t bin_id, Discontinuity_Graph<k>& d_graph, op_buf_t& op_buf);

    Subgraph(const Subgraph&) = delete;
    Subgraph(Subgraph&&) = delete;

    // Initialize the map collection of different workers.
    static void init_maps();

    // Free the map collection from memory of different workers.
    static void free_maps();

    // Constructs the subgraph from the KMC bin into an internal navigable and
    // membership data structure.
    void construct();

    // Constructs the subgraph from the KMC bin into an internal navigable and
    // membership data structure. Addresses "exact" loop-filtering opposed to
    // `construct`.
    void construct_loop_filtered();

    // Builds the compacted graph from the original graph.
    void contract();

    // Returns the size of the graph.
    std::size_t size() const;

    // Returns the count of isolated vertices—not part of any edge.
    uint64_t isolated_vertex_count() const;

    // Returns the number of (multi-)edges in the graph.
    uint64_t edge_count() const;

    // Returns the number of edges of the discontinuity graph produced from this
    // subgraph.
    uint64_t discontinuity_edge_count() const;

    // Returns the number of trivial maximal unitigs in the graph (i.e. also
    // maximal unitigs in the supergraph).
    uint64_t trivial_mtig_count() const;

    // Returns the number of trivial maximal unitigs in the graph that are ICCs.
    uint64_t icc_count() const;

    // Returns the total number of characters in the literal representations of
    // all the maximal unitigs.
    uint64_t label_size() const;
};


// Type of scenarios how a unitig-walk terminates in the subgraph.
template <uint16_t k>
enum class Subgraph<k>::Walk_Termination
{
    null,       // non-existent walk
    branched,   // branched off
    crossed,    // crossed to a different unitig, or looped / cycled back to the same unitig
    dead_ended, // no extension existed
    exitted,    // exitted the subgraph
};


template <uint16_t k>
inline bool Subgraph<k>::extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig)
{
    constexpr auto back = side_t::back;
    constexpr auto front = side_t::front;
    constexpr auto exitted = termination_t::exitted;

    assert(M.find(v_hat) != M.end());
    const auto& state = M[v_hat];
    if(state.is_visited())  // The containing maximal unitig has already been outputted.
        return false;


    maximal_unitig.mark_linear();

    Directed_Vertex<k> v_l, v_r;    // Possible discontinuity ends of the maximal unitig at the left and the right extensions.
    termination_t walk_end_l(termination_t::null), walk_end_r(termination_t::null); // Whether the maximal unitig tried to exit the subgraph through the left and the right extensions.

    walk_end_r = walk_unitig(v_hat, state, back, maximal_unitig.unitig(back), v_r);
    if(maximal_unitig.unitig(back).is_cycle())
    {
        assert(walk_end_r == termination_t::crossed);
        maximal_unitig.mark_cycle(back);
    }
    else
    {
        walk_end_l = walk_unitig(v_hat, state, front, maximal_unitig.unitig(front), v_l);
        assert(!maximal_unitig.unitig(front).is_cycle());
    }

    if(walk_end_l == exitted || walk_end_r == exitted)  // The maximal unitig containing `v_hat` spans multiple subgraphs.
    {
        maximal_unitig.finalize_weak();
        G.add_edge( walk_end_l == exitted ? v_l.canonical() : Discontinuity_Graph<k>::phi(), walk_end_l == exitted ? v_l.entrance_side() : side_t::back,
                    walk_end_r == exitted ? v_r.canonical() : Discontinuity_Graph<k>::phi(), walk_end_r == exitted ? v_r.entrance_side() : side_t::back,
                    walk_end_l != exitted, walk_end_r != exitted, maximal_unitig);
        disc_edge_c++;
    }
    else    // Extracted a trivial maximal unitig.
    {
        trivial_mtig_c++;
        if(maximal_unitig.is_cycle())
            icc_count_++;

        maximal_unitig.finalize();
        maximal_unitig.add_fasta_rec_to_buffer(op_buf);
    }

    return true;
}


template <uint16_t k>
inline typename Subgraph<k>::termination_t Subgraph<k>::walk_unitig(const Kmer<k>& v_hat, const State_Config& st_v, const side_t s_v_hat, Unitig_Scratch<k>& unitig, Directed_Vertex<k>& exit_v)
{
    const auto s_icc_return = inv_side(s_v_hat);    // The side through which to return to `v_hat` if it's contained in an ICC.
    Directed_Vertex<k> v(s_v_hat == side_t::back ? v_hat : v_hat.reverse_complement()); // Current vertex being added to the unitig.
    side_t s_v = s_v_hat;   // The side of the current vertex `v_hat` through which to extend the unitig, i.e. to exit `v`.
    State_Config state = st_v;  // State of `v`.
    base_t b_ext;   // The nucleobase encoding the edge(s) incident to the side `s_v` of `v`.

    unitig.init(v);

    auto it = M.find(v.canonical());
    // auto it = M.find_positive(v.canonical());
    while(true)
    {
        it->second.mark_visited();

        b_ext = state.edge_at(s_v);
        assert(!state.is_discontinuous(s_v) || b_ext == base_t::E); // If a side is discontinuous, it must be empty.
        if(b_ext == base_t::N)  // Reached a branching endpoint.
            return termination_t::branched;

        if(b_ext == base_t::E)
        {
            if(!state.is_discontinuous(s_v))    // Reached a truly empty side.
                return termination_t::dead_ended;

            // Trying to exit the subgraph through a discontinuity vertex.
            exit_v = v;
            return termination_t::exitted;
        }


        if(s_v == side_t::front)
            b_ext = DNA_Utility::complement(b_ext);

        v.roll_forward(b_ext);  // Walk to the next vertex.
        assert(M.find(v.canonical()) != M.end());
        it = M.find(v.canonical());
        // it = M.find_positive(v.canonical());
        state = it->second;

        s_v = v.entrance_side();
        assert(!state.is_empty_side(s_v));
        if(state.is_branching_side(s_v))    // Crossed an endpoint and reached a different unitig.
            return termination_t::crossed;

        if(state.is_visited())  // Hit the same unitig.
        {
            // The unitig is an ICC; crossed back to the same unitig.
            if(v.canonical() == v_hat && s_v == s_icc_return)
                unitig.mark_cycle();
            else
                // Otherwise, hit a looping edge—visiting the immediate predecessor vertex.
                // A special case; crossed to the same unitig from a different orientation.
                assert(std::as_const(v)
                        .roll_backward(s_v == side_t::front ? state.edge_at(s_v) : DNA_Utility::complement(state.edge_at(s_v)))
                        .is_same_vertex(v));

            return termination_t::crossed;
        }


        // Still within the unitig.
        const bool e = unitig.extend(v, DNA_Utility::map_char(b_ext));
        assert(e); (void)e;

        s_v = opposite_side(s_v);
    }

    return termination_t::null;
}

}




#endif
