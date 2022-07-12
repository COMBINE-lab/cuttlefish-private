
#include "Kmer_Index.hpp"
#include "Ref_Parser.hpp"

#include <vector>
#include <unordered_map>


template <uint16_t k>
template <uint16_t l>
inline bool Kmer_Index<k>::validate(const std::string& file_path)
{
    // Build a naive and a cuttlefish index.
    Ref_Parser parser(file_path);

    std::string paths;  // The concatenated paths.
    std::vector<std::size_t> ends;  // Endpoints of the paths in the concatenated sequence.

    Kmer_Index<k> kmer_index(l, 1, true);   // The cuttlefish index.
    const auto token = kmer_index.get_token();

    std::unordered_map<cuttlefish::minimizer_t, std::vector<std::size_t>> M;    // The minimizer table.

    while(parser.read_next_seq())
    {
        if(parser.seq_len() < k)
            continue;

        const char* const seq = parser.seq();
        const std::size_t len = parser.seq_len();
        std::size_t last_min_idx = len;

        kmer_index.deposit(token, seq, len);

        // Collect the minimizer instances.
        for(std::size_t kmer_idx = 0; kmer_idx + k <= len; ++kmer_idx)
        {
            Kmer<l> min_lmer(seq, kmer_idx);
            std::size_t min_idx = kmer_idx;
            uint64_t min_hash = Minimizer_Iterator::hash(min_lmer.as_int());

            for(std::size_t i = kmer_idx + 1; i + l <= kmer_idx + k; ++i)
            {
                const Kmer<l> lmer(seq, i);
                const uint64_t lmer_hash = Minimizer_Iterator::hash(lmer.as_int());

                if(min_hash < lmer_hash)
                    continue;

                if(min_hash > lmer_hash || min_lmer > lmer)
                    min_lmer = lmer, min_idx = i, min_hash = lmer_hash;
            }

            const std::size_t abs_min_idx = paths.size() + min_idx;
            if(abs_min_idx != last_min_idx)
            {
                M[min_lmer.as_int()].emplace_back(abs_min_idx);
                last_min_idx = abs_min_idx;
            }
        }

        paths += seq;
        ends.emplace_back(paths.size());
    }


    kmer_index.index();


    return true;
}
