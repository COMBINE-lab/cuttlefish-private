
#ifndef SUBGRAPH_HPP
#define SUBGRAPH_HPP



#include "State_Config.hpp"
#include "Directed_Vertex.hpp"
#include "Kmer.hpp"
#include "DNA.hpp"
#include "DNA_Utility.hpp"
#include "Unitig_Scratch.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "Edge_Matrix.hpp"
#include "dBG_Utilities.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
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

    const std::string graph_bin_dir_path;   // Path to the directory with all the graph KMC-bins.
    const std::size_t bin_id; // ID of the graph KMC-bin.

    typedef std::unordered_map<Kmer<k>, State_Config, Kmer_Hasher<k>> map_t;
    map_t M;

    uint64_t edge_c;    // Number of edges in the graph.
    uint64_t label_sz;  // Total number of characters in the literal representations of all the maximal unitigs.
    uint64_t isolated;  // Count of isolated vertices—not part of any edge.

    Edge_Matrix<k>& E;  // Edge-matrix of the discontinuity graph.


    // Extracts the maximal unitig containing the vertex `v_hat`, and
    // `maximal_unitig` is used as the working scratch for the extraction, i.e.
    // to build and store two unitigs connecting to the two sides of `v_hat`.
    // Returns `true` iff the containing maximal unitig has not been outputted
    // earlier.
    bool extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig);

    // Traverses a unitig starting from the vertex `v_hat`, exiting it through
    // the side `s_v_hat`. `st_v` is the state of `s_v_hat`. `unitig` is used as
    // the scratch space to build the unitig.
    void walk_unitig(const Kmer<k>& v_hat, const State_Config& v_inf, side_t s_v_hat, Unitig_Scratch<k>& unitig);


public:

    // Constructs a subgraph object for the `bin_id`'th bin in the graph bin
    // directory `bin_dir_path`. Updates the edge-matrix `E` of the
    // discontinuity graph with its edges observed from this subgraph.
    Subgraph(const std::string& bin_dir_path, std::size_t bin_id, Edge_Matrix<k>& E);

    Subgraph(const Subgraph&) = delete;
    Subgraph(Subgraph&&) = delete;

    // Constructs the subgraph from the KMC bin into an internal navigable and
    // membership data structure.
    void construct();

    // Builds the compacted graph from the original graph.
    void contract();

    // Returns the size of the graph.
    std::size_t size() const;

    // Returns the count of isolated vertices—not part of any edge.
    uint64_t isolated_vertex_count() const;

    // Returns the number of (multi-)edges in the graph.
    uint64_t edge_count() const;

    // Returns the total number of characters in the literal representations of
    // all the maximal unitigs.
    uint64_t label_size() const;
};


template <uint16_t k>
inline bool Subgraph<k>::extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig)
{
    constexpr auto back = side_t::back;
    constexpr auto front = side_t::front;

    assert(M.find(v_hat) != M.end());

    const auto& state = M[v_hat];
    if(state.is_visited())  // The containing maximal unitig has already been outputted.
        return false;

    maximal_unitig.mark_linear();

    walk_unitig(v_hat, state, back, maximal_unitig.unitig(back));
    if(maximal_unitig.unitig(back).is_cycle())
        maximal_unitig.mark_cycle(back);
    else
        walk_unitig(v_hat, state, front, maximal_unitig.unitig(front));

    maximal_unitig.finalize();

    return true;
}


template <uint16_t k>
inline void Subgraph<k>::walk_unitig(const Kmer<k>& v_hat, const State_Config& st_v, const side_t s_v_hat, Unitig_Scratch<k>& unitig)
{
    Directed_Vertex<k> v(s_v_hat == side_t::back ? v_hat : v_hat.reverse_complement()); // Current vertex being added to the unitig.
    side_t s_v = s_v_hat;   // The side of the current vertex `v_hat` through which to extend the unitig, i.e. to exit `v`.
    State_Config state = st_v;  // State of `v`.
    base_t b_ext;   // The nucleobase encoding the edge(s) incident to the side `s_v` of `v`.

    unitig.init(v);

    while(true)
    {
        M[v.canonical()].mark_visited();

        b_ext = state.edge_at(s_v);
        if(b_ext == base_t::N)  // Reached an endpoint.
            break;

        if(s_v == side_t::front)
            b_ext = DNA_Utility::complement(b_ext);

        v.roll_forward(b_ext);  // Walk to the next vertex.
        assert(M.find(v.canonical()) != M.end());
        state = M[v.canonical()];

        if(state.is_visited())  // Hit a looping edge—visiting the immediate predecessor vertex.
            break;

        s_v = v.entrance_side();
        assert(!state.is_empty_side(s_v));
        if(state.is_branching_side(s_v))    // Crossed an endpoint and reached a different unitig.
            break;

        // Still within the unitig.
        if(!unitig.extend(v, DNA_Utility::map_char(b_ext)))
            break;  // The unitig is a DCC (Detached Chordless Cycle).

        s_v = opposite_side(s_v);
    }
}

}




#endif
