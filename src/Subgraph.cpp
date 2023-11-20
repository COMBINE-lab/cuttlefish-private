
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
{}


template <uint16_t k>
void Subgraph<k>::load()
{
    const std::size_t q_sz = 1; // Number of internal queues for the iterator; set to the count of workers using the iterator.
    IterateSuperKmers super_kmer_it(graph_bin_dir_path, bin_id, q_sz);  // The iterator.

    typedef unsigned long long super_kmer_data_t;   // KMC's super k-mer's representation type.
    const auto word_count = super_kmer_it.GetSuperKmerDataLen();    // Fixed number of words in KMC super k-mer.

    // Returns the `idx`'th base in a KMC super k-mer `super_kmer`.
    const auto get_base =
        [word_count](const uint64_t* const super_kmer, const std::size_t idx)
        {
            assert(idx / 32 < word_count);
            return DNA::Base((super_kmer[word_count - 1 - (idx >> 5)] >> ((31 - (idx & 31)) << 1)) & uint64_t(0b11));
        };

    std::ofstream output(graph_bin_dir_path + "kmers." + std::to_string(parlay::worker_id()), std::ios::app);
    Directed_Vertex<k> v;
    super_kmer_it.AddConsumer(
        [&](const super_kmer_data_t* const super_kmer, const std::size_t len)
        {
            auto const kmc_data = reinterpret_cast<const uint64_t*>(super_kmer);
            v.from_KMC_super_kmer(kmc_data, word_count);
            std::size_t kmer_idx = 0;

            while(true)
            {
                assert(kmer_idx + k - 1 < len);
                output << ">" << kmer_idx << "\n" << v.kmer() << "\n";

                const auto is_canonical = v.in_canonical_form();
                const auto pred_base = (kmer_idx == 0 ? DNA::Base::N : get_base(kmc_data, kmer_idx - 1));
                const auto succ_base = (kmer_idx + k == len ? DNA::Base::N : get_base(kmc_data, kmer_idx + k));
                const auto front = (is_canonical ? pred_base : succ_base);
                const auto back  = (is_canonical ? succ_base : pred_base);

                // Update hash table with the neighborhood info.
                auto const it = M.find(v.canonical());
                if(it == M.end())
                    M.emplace(v.canonical(), Vertex_Info());

                M[v.canonical()].add_neighbor(front, back);

                if(kmer_idx + k == len)
                    break;

                v.roll_forward(succ_base);
                kmer_idx++;
            }
        }
    );

    super_kmer_it.WaitForAll();

    output.close();
}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Subgraph)
