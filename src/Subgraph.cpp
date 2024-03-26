
#include "Subgraph.hpp"
#include "Super_Kmer_Bucket.hpp"
#include "DNA.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include "fstream"


namespace cuttlefish
{

template <uint16_t k, bool Colored_> std::vector<Padded_Data<typename Subgraph<k, Colored_>::map_t>> Subgraph<k, Colored_>::map;


template <uint16_t k, bool Colored_>
Subgraph<k, Colored_>::Subgraph(const Super_Kmer_Bucket<Colored_>& B, Discontinuity_Graph<k>& G, op_buf_t& op_buf):
      B(B)
    , M(map[parlay::worker_id()].data())
    , edge_c(0)
    , label_sz(0)
    , disc_edge_c(0)
    , isolated(0)
    , G(G)
    , trivial_mtig_c(0)
    , icc_count_(0)
    , op_buf(op_buf)
{
    M.clear();
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::init_maps()
{
    map.resize(parlay::num_workers());
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::free_maps()
{
    for(auto& M : map)
        force_free(M.data());

    map.clear();
    map.shrink_to_fit();
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::construct()
{
    auto super_kmer_it = B.iterator();  // Iterator over the weak super k-mers inducing this graph.

    typedef typename decltype(super_kmer_it)::label_unit_t label_unit_t;
    const auto word_count = super_kmer_it.super_kmer_word_count();  // Fixed number of words in a super k-mer label.

    Directed_Vertex<k> v;

    Super_Kmer_Attributes<Colored_> att;
    label_unit_t* label;
    while(super_kmer_it.next(att, label))
    {
        const auto len = att.len();
        assert(len >= k);
        assert(len < 2 * (k - 1));

        v.from_super_kmer(label, word_count);
        std::size_t kmer_idx = 0;
        while(true)
        {
            assert(kmer_idx + k - 1 < len);

            const auto is_canonical = v.in_canonical_form();
            const auto pred_base = (kmer_idx == 0 ? base_t::E : get_base(label, word_count, kmer_idx - 1));
            const auto succ_base = (kmer_idx + k == len ? base_t::E : get_base(label, word_count, kmer_idx + k));
            const auto front = (is_canonical ? pred_base : DNA_Utility::complement(succ_base));
            const auto back  = (is_canonical ? succ_base : DNA_Utility::complement(pred_base));

            edge_c += (succ_base != base_t::E);

            // Update hash table with the neighborhood info.
            auto& st = M[v.canonical()];
            st.update_edges(front, back);
            if(kmer_idx == 0 && att.left_discontinuous())
                st.mark_discontinuous(v.entrance_side());

            if(kmer_idx + k == len)
            {
                if(att.right_discontinuous())
                    st.mark_discontinuous(v.exit_side());

                break;
            }

            if((kmer_idx > 0 || !att.left_discontinuous()) && (kmer_idx + k < len || !att.right_discontinuous()))
                assert(!st.is_discontinuity());


            v.roll_forward(succ_base);
            kmer_idx++;
        }
    }
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::construct_loop_filtered()
{
    auto super_kmer_it = B.iterator();  // Iterator over the weak super k-mers inducing this graph.

    const auto word_count = super_kmer_it.super_kmer_word_count();  // Fixed number of words in a super k-mer label.

    Directed_Vertex<k> v_pre, v_suf;

    Super_Kmer_Attributes<Colored_> att;
    label_unit_t* label;
    while(super_kmer_it.next(att, label))
    {
        const auto len = att.len();
        assert(len >= k);
        assert(len < 2 * (k - 1));

        if(len == k)
            continue;

        std::size_t edge_idx = 0;
        v_pre.from_super_kmer(label, word_count);
        auto succ_base_pre = get_base(label, word_count, edge_idx + k);
        v_suf = std::as_const(v_pre).roll_forward(succ_base_pre);
        auto pred_base_suf = get_base(label, word_count, edge_idx);

        while(true)
        {
            assert(edge_idx + k < len);

            auto& st_pre = M[v_pre.canonical()];
            auto& st_suf = M[v_suf.canonical()];

            v_pre.in_canonical_form() ? st_pre.update_edge(side_t::back, succ_base_pre) : st_pre.update_edge(side_t::front, DNA_Utility::complement(succ_base_pre));
            if(!v_suf.is_same_vertex(v_pre))    // Avoid double-counting of self-loops.
                v_suf.in_canonical_form() ? st_suf.update_edge(side_t::front, pred_base_suf) : st_suf.update_edge(side_t::back, DNA_Utility::complement(pred_base_suf));
            edge_c++;

            if(edge_idx == 0 && att.left_discontinuous())
                st_pre.mark_discontinuous(v_pre.entrance_side());

            if(edge_idx + 1 + k == len)
            {
                if(att.right_discontinuous())
                    st_suf.mark_discontinuous(v_suf.exit_side());

                break;
            }


            edge_idx++;
            v_pre = v_suf;
            succ_base_pre = get_base(label, word_count, edge_idx + k);
            pred_base_suf = get_base(label, word_count, edge_idx);
            v_suf.roll_forward(succ_base_pre);
        }
    }
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::contract()
{
    Maximal_Unitig_Scratch<k> maximal_unitig;   // Scratch space to be used to construct maximal unitigs.
    uint64_t vertex_count = 0;  // Count of vertices processed.
    uint64_t unitig_count = 0;  // Count of maximal unitigs.
    uint64_t non_isolated = 0;  // Count of non-isolated vertices.

    for(const auto& p : M)
    {
        const auto& v = p.first;
        const auto& v_st = p.second;
        assert(!v_st.is_discontinuous(side_t::front) || !v_st.is_discontinuous(side_t::back));

        if(v_st.is_visited())   // The containing maximal unitig has already been outputted.
            continue;

        if(v_st.is_isolated())
        {
            if(v_st.is_discontinuity()) // A potential phantom-edge for the discontinuity graph is incident to `v`.
                G.inc_potential_phantom_edge();

            isolated++;
            continue;
        }

        non_isolated++;


        if(extract_maximal_unitig(v, maximal_unitig))
        {
            vertex_count += maximal_unitig.size();
            unitig_count++;
            label_sz += maximal_unitig.size() + k - 1;
        }
    }

    assert(vertex_count + isolated == M.size());
}


template <uint16_t k, bool Colored_>
std::size_t Subgraph<k, Colored_>::size() const
{
    return M.size();
}


template <uint16_t k, bool Colored_>
uint64_t Subgraph<k, Colored_>::edge_count() const
{
    return edge_c;
}


template <uint16_t k, bool Colored_>
uint64_t Subgraph<k, Colored_>::discontinuity_edge_count() const
{
    return disc_edge_c;
}


template <uint16_t k, bool Colored_>
uint64_t Subgraph<k, Colored_>::trivial_mtig_count() const
{
    return trivial_mtig_c;
}


template <uint16_t k, bool Colored_>
uint64_t Subgraph<k, Colored_>::icc_count() const
{
    return icc_count_;
}


template <uint16_t k, bool Colored_>
uint64_t Subgraph<k, Colored_>::label_size() const
{
    return label_sz;
}


template <uint16_t k, bool Colored_>
uint64_t Subgraph<k, Colored_>::isolated_vertex_count() const
{
    return isolated;
}

}

/*
```
for b in B:    // B: collection of buckets
    for s in b:    // s: super k-mer
        i = 0
        for x in s:    // x: k-mer
            is_can = x.is_canonical()
            pn = pred_nuc(x) // or, pn = s[i - 1]
            sn = succ_nuc(x) // or, sn = s[i + k]
            front = (is_can ? pn : bar(sn))
            back  = (is_can ? sn : bar(pn))
            update HT for x.canonical() with front and back
            i = i + 1
```
*/



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Subgraph)
