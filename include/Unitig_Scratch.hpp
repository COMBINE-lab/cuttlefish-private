
#ifndef UNITIG_SCRATCH_HPP
#define UNITIG_SCRATCH_HPP



#include "Directed_Vertex.hpp"
#include "dBG_Utilities.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>


// =============================================================================
// A class to keep scratch data (i.e. working space) for unitigs.
template <uint16_t k>
class Unitig_Scratch
{
private:

    Directed_Vertex<k> anchor;      // The anchor vertex of the unitig traversal.
    Directed_Vertex<k> endpoint_;   // The current end of the unitig through which farther extensions can be done.
                                    // (The side for the extension is to be handled by the client code, although can
                                    // also be inferred from the "directed" vertex.)
    Directed_Vertex<k> min_vertex_; // The lexicographically minimum vertex in the unitig.
    std::size_t vertex_idx;         // Index of the vertex in the path being traversed.
    std::size_t min_v_idx;          // Index of the lexicographically minimum vertex in the path.

    std::vector<char> label_;       // Literal label of the unitig.
    std::vector<uint64_t> hash_;    // Hashes of the constituent vertices of the unitig.
    std::vector<Kmer<k>> V;         // Set of the vertices (in their canonical-form) in the unitig.
    bool is_cycle_;                 // Whether the unitig is cyclical or not.


    // Clears the scratch data.
    void clear();

public:

    // Initializes the unitig scratch with the vertex `v`. `Keep_Vertices`
    // denote whether the vertices in their canonical form are to be kept or
    // not.
    template <bool Keep_Vertices_ = false>
    void init(const Directed_Vertex<k>& v);

    // Extends the unitig scratch with the vertex `v`, and its literal form
    // with the symbol `b`. Returns `true` iff adding `v` to the unitig does
    // not render itself a cycle. `Keep_Vertices` denote whether the vertices
    // in their canonical form are to be kept or not.
    template <bool Keep_Vertices_ = false>
    bool extend(const Directed_Vertex<k>& v, char b);

    // Marks the unitig as a cycle.
    void mark_cycle() { is_cycle_ = true; }

    // Reverse complements the unitig.
    void reverse_complement();

    // Returns the literal label of the unitig.
    const std::vector<char>& label() const { return label_; }

    // Returns the hash collection of the unitig vertices.
    const std::vector<uint64_t>& hash() const { return hash_; }

    // Returns the vertices (in their canonical-form) in the unitig, in the
    // order of the label.
    const std::vector<Kmer<k>>& vertices() const { return V; }

    // Returns the current extension-end vertex of the unitig.
    const Directed_Vertex<k>& endpoint() const { return endpoint_; }

    // Returns the count of vertices in this unitig.
    std::size_t size() const { return hash_.size(); }

    // Returns `true` iff unitig is a cycle.
    bool is_cycle() const { return is_cycle_; }

    // Returns the lexicographically minimum vertex in the unitig.
    const Directed_Vertex<k>& min_vertex() const { return min_vertex_; }

    // Returns the index of the lexicographically minimum vertex in the unitig.
    std::size_t min_vertex_idx() const { return min_v_idx; }
};


template <uint16_t k>
inline void Unitig_Scratch<k>::clear()
{
    label_.clear();
    hash_.clear();
    V.clear();
}


template <uint16_t k>
template <bool Keep_Vertices_>
inline void Unitig_Scratch<k>::init(const Directed_Vertex<k>& v)
{
    clear();

    min_vertex_ = endpoint_ = anchor = v;
    min_v_idx = vertex_idx = 0;

    endpoint_.kmer().get_label(label_);
    hash_.emplace_back(endpoint_.hash());
    if constexpr(Keep_Vertices_) V.emplace_back(v.canonical());
    is_cycle_ = false;
}


template <uint16_t k>
template <bool Keep_Vertices_>
inline bool Unitig_Scratch<k>::extend(const Directed_Vertex<k>& v, const char b)
{
    if(v.is_same_vertex(anchor))
    {
        is_cycle_ = true;
        return false;
    }


    endpoint_ = v;
    vertex_idx++;

    if(min_vertex_.canonical() > endpoint_.canonical())
        min_vertex_ = endpoint_,
        min_v_idx = vertex_idx;

    label_.emplace_back(b);
    hash_.emplace_back(endpoint_.hash());
    if constexpr(Keep_Vertices_) V.emplace_back(v.canonical());

    return true;
}


template <uint16_t k>
inline void Unitig_Scratch<k>::reverse_complement()
{
    cuttlefish::reverse_complement(label_);
    std::reverse(hash_.begin(), hash_.end());
    std::reverse(V.begin(), V.end());
    min_v_idx = (size() - 1 - min_v_idx);
}



#endif
