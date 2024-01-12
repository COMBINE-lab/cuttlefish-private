
#include "Subgraph.hpp"
#include "DNA.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"
#include "kmc-super-kmers-iterator/iterate_super_kmers.h"

#include "fstream"


namespace cuttlefish
{

template<uint16_t k> std::vector<Padded_Data<typename Subgraph<k>::map_t>> Subgraph<k>::map;


template <uint16_t k>
Subgraph<k>::Subgraph(const std::string& bin_dir_path, const std::size_t bin_id, Discontinuity_Graph<k>& G, op_buf_t& op_buf):
      graph_bin_dir_path(bin_dir_path)
    , bin_id(bin_id)
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


template <uint16_t k>
void Subgraph<k>::init_maps()
{
    map.resize(parlay::num_workers());
}


template <uint16_t k>
void Subgraph<k>::free_maps()
{
    for(auto& M : map)
        force_free(M.data());

    map.clear();
    map.shrink_to_fit();
}


template <uint16_t k>
void Subgraph<k>::construct()
{
    const std::size_t q_sz = 1; // Number of internal queues for the iterator; set to the count of workers using the iterator.
    IterateSuperKmers super_kmer_it(graph_bin_dir_path, bin_id, q_sz);  // The iterator over the super k-mers in this graph-bin.

    typedef unsigned long long super_kmer_data_t;   // KMC's super k-mer's representation type.
    const auto word_count = super_kmer_it.GetSuperKmerDataLen();    // Fixed number of words in KMC super k-mer.

    // Returns the `idx`'th base in a KMC super k-mer `super_kmer`.
    const auto get_base =
        [word_count](const uint64_t* const super_kmer, const std::size_t idx)
        {
            assert(idx / 32 < word_count);
            return DNA::Base((super_kmer[word_count - 1 - (idx >> 5)] >> ((31 - (idx & 31)) << 1)) & uint64_t(0b11));
        };

    Directed_Vertex<k> v;

    // Extracts and processes each k-mer (vertex) from a given super k-mer
    // `super_kmer` of length `len` bases. `disc_l` and `disc_r` denote whether
    // the left- and the right-end k-mers are discontinuity vertices or not.
    const auto extract_kmers =
        [&](const super_kmer_data_t* const super_kmer, const std::size_t len, const bool disc_l, const bool disc_r)
        {
            assert(len >= k);

            auto const kmc_data = reinterpret_cast<const uint64_t*>(super_kmer);
            v.from_KMC_super_kmer(kmc_data, word_count);
            std::size_t kmer_idx = 0;

            while(true)
            {
                assert(kmer_idx + k - 1 < len);

                const auto is_canonical = v.in_canonical_form();
                const auto pred_base = (kmer_idx == 0 ? base_t::E : get_base(kmc_data, kmer_idx - 1));
                const auto succ_base = (kmer_idx + k == len ? base_t::E : get_base(kmc_data, kmer_idx + k));
                const auto front = (is_canonical ? pred_base : DNA_Utility::complement(succ_base));
                const auto back  = (is_canonical ? succ_base : DNA_Utility::complement(pred_base));

                edge_c += (succ_base != base_t::E);

                // Update hash table with the neighborhood info.
                auto& st = M[v.canonical()];
                st.update_edges(front, back);
                if(kmer_idx == 0 && disc_l)
                    st.mark_discontinuous(v.entrance_side());

                if(kmer_idx + k == len)
                {
                    if(disc_r)
                        st.mark_discontinuous(v.exit_side());

                    break;
                }


                if((kmer_idx > 0 || !disc_l) && (kmer_idx + k < len || !disc_r))
                    assert(!st.is_discontinuity());


                v.roll_forward(succ_base);
                kmer_idx++;
            }
        };

    super_kmer_it.AddConsumer(extract_kmers);
    super_kmer_it.WaitForAll();
}


template <uint16_t k>
void Subgraph<k>::construct_loop_filtered()
{
    const std::size_t q_sz = 1; // Number of internal queues for the iterator; set to the count of workers using the iterator.
    IterateSuperKmers super_kmer_it(graph_bin_dir_path, bin_id, q_sz);  // The iterator over the super k-mers in this graph-bin.

    typedef unsigned long long super_kmer_data_t;   // KMC's super k-mer's representation type.
    const auto word_count = super_kmer_it.GetSuperKmerDataLen();    // Fixed number of words in KMC super k-mer.

    // Returns the `idx`'th base in a KMC super k-mer `super_kmer`.
    const auto get_base =
        [word_count](const uint64_t* const super_kmer, const std::size_t idx)
        {
            assert(idx / 32 < word_count);
            return DNA::Base((super_kmer[word_count - 1 - (idx >> 5)] >> ((31 - (idx & 31)) << 1)) & uint64_t(0b11));
        };

    Directed_Vertex<k> v_pre, v_suf;

    // Extracts and processes each k-mer (vertex) from a given super k-mer
    // `super_kmer` of length `len` bases. `disc_l` and `disc_r` denote whether
    // the left- and the right-end k-mers are discontinuity vertices or not.
    const auto extract_kmers =
        [&](const super_kmer_data_t* const super_kmer, const std::size_t len, const bool disc_l, const bool disc_r)
        {
            assert(len >= k);

            if(len == k)
                return;

            auto const kmc_data = reinterpret_cast<const uint64_t*>(super_kmer);
            std::size_t edge_idx = 0;
            v_pre.from_KMC_super_kmer(kmc_data, word_count);
            auto succ_base_pre = get_base(kmc_data, edge_idx + k);
            v_suf = std::as_const(v_pre).roll_forward(succ_base_pre);
            auto pred_base_suf = get_base(kmc_data, edge_idx);

            while(true)
            {
                assert(edge_idx + k < len);

                auto it_l = M.find(v_pre.canonical());
                if(it_l == M.end())
                    it_l = M.emplace(v_pre.canonical(), State_Config()).first;
                auto& st_pre = it_l->second;

                auto it_r = M.find(v_suf.canonical());
                if(it_r == M.end())
                    it_r = M.emplace(v_suf.canonical(), State_Config()).first;
                auto& st_suf = it_r->second;

                v_pre.in_canonical_form() ? st_pre.update_edge(side_t::back, succ_base_pre) : st_pre.update_edge(side_t::front, DNA_Utility::complement(succ_base_pre));
                if(!v_suf.is_same_vertex(v_pre))    // Avoid double-counting of self-loops.
                    v_suf.in_canonical_form() ? st_suf.update_edge(side_t::front, pred_base_suf) : st_suf.update_edge(side_t::back, DNA_Utility::complement(pred_base_suf));
                edge_c++;

                if(edge_idx == 0 && disc_l)
                    st_pre.mark_discontinuous(v_pre.entrance_side());

                if(edge_idx + 1 + k == len)
                {
                    if(disc_r)
                        st_suf.mark_discontinuous(v_suf.exit_side());

                    break;
                }


                edge_idx++;
                v_pre = v_suf;
                succ_base_pre = get_base(kmc_data, edge_idx + k);
                pred_base_suf = get_base(kmc_data, edge_idx);
                v_suf.roll_forward(succ_base_pre);
            }
        };

    super_kmer_it.AddConsumer(extract_kmers);
    super_kmer_it.WaitForAll();
}


template <uint16_t k>
void Subgraph<k>::contract()
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


template <uint16_t k>
std::size_t Subgraph<k>::size() const
{
    return M.size();
}


template <uint16_t k>
uint64_t Subgraph<k>::edge_count() const
{
    return edge_c;
}


template <uint16_t k>
uint64_t Subgraph<k>::discontinuity_edge_count() const
{
    return disc_edge_c;
}


template <uint16_t k>
uint64_t Subgraph<k>::trivial_mtig_count() const
{
    return trivial_mtig_c;
}


template <uint16_t k>
uint64_t Subgraph<k>::icc_count() const
{
    return icc_count_;
}


template <uint16_t k>
uint64_t Subgraph<k>::label_size() const
{
    return label_sz;
}


template <uint16_t k>
uint64_t Subgraph<k>::isolated_vertex_count() const
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
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraph)
