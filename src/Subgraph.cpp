
#include "Subgraph.hpp"
#include "DNA.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"
#include "kmc-super-kmers-iterator/iterate_super_kmers.h"

#include "fstream"


namespace cuttlefish
{

template <uint16_t k>
Subgraph<k>::Subgraph(const std::string& bin_dir_path, const std::size_t bin_id):
      graph_bin_dir_path(bin_dir_path)
    , bin_id(bin_id)
    , edge_c(0)
{}


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

    // std::ofstream output(graph_bin_dir_path + "kmers." + std::to_string(parlay::worker_id()), std::ios::app);
    Directed_Vertex<k> v;

    // Extracts and processes each k-mer (vertex) from a given super k-mer
    // `super_kmer` of length `len`.
    const auto extract_kmers =
        [&](const super_kmer_data_t* const super_kmer, const std::size_t len)
        {
            auto const kmc_data = reinterpret_cast<const uint64_t*>(super_kmer);
            v.from_KMC_super_kmer(kmc_data, word_count);
            std::size_t kmer_idx = 0;

            while(true)
            {
                assert(kmer_idx + k - 1 < len);
                // output << ">" << kmer_idx << "\n" << v.kmer() << "\n";

                const auto is_canonical = v.in_canonical_form();
                const auto pred_base = (kmer_idx == 0 ? DNA::Base::N : get_base(kmc_data, kmer_idx - 1));
                const auto succ_base = (kmer_idx + k == len ? DNA::Base::N : get_base(kmc_data, kmer_idx + k));
                const auto front = (is_canonical ? pred_base : DNA_Utility::complement(succ_base));
                const auto back  = (is_canonical ? succ_base : DNA_Utility::complement(pred_base));

                edge_c += (succ_base != DNA::Base::N);

                // Update hash table with the neighborhood info.
                const auto it = M.find(v.canonical());
                if(it == M.end())
                    M.emplace(v.canonical(), State_Config());

                M[v.canonical()].update_edges(front, back);

                if(kmer_idx + k == len)
                    break;

                v.roll_forward(succ_base);
                kmer_idx++;
            }
        };

    super_kmer_it.AddConsumer(extract_kmers);
    super_kmer_it.WaitForAll();

    // output.close();
}


template <uint16_t k>
void Subgraph<k>::compact()
{
    Maximal_Unitig_Scratch<k> maximal_unitig;   // Scratch space to be used to construct maximal unitigs.
    uint64_t vertex_count = 0;  // Number of vertices processed.

    for(const auto& p : M)
    {
        const auto& v = p.first;
        // TODO: extract maximal unitig containing `v`.
        if(extract_maximal_unitig(v, maximal_unitig))
            vertex_count += maximal_unitig.size();
    }

    assert(vertex_count == M.size());
}


template <uint16_t k>
std::size_t Subgraph<k>::size() const
{
    return M.size();
}


template <uint16_t k>
std::size_t Subgraph<k>::edge_count() const
{
    return edge_c;
}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraph)
