
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
    std::size_t inst_count = 0;

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

            if(min_idx != last_min_idx)
            {
                const std::size_t abs_min_idx = paths.size() + min_idx;
                M[min_lmer.as_int()].emplace_back(abs_min_idx);
                inst_count++;
                last_min_idx = min_idx;
            }
        }

        paths += seq;
        ends.emplace_back(paths.size());
    }


    kmer_index.index();

    std::cout << "Constructed the naive and the Cuttlefish index.\n";

    std::cout <<    "\n\nCross-checking the indices.\n"
                    "===============================\n\n";

    // Check paths' validity.

    std::cout << "Path counts:\n\tNaive idx: " << ends.size() << ", Cuttlefish idx: " << kmer_index.path_count << ".\n";
    if(ends.size() != kmer_index.path_count)
        return false;

    for(std::size_t i = 0; i < ends.size(); ++i)
        if(ends[i] != kmer_index.path_ends->at(i))
            return false;

    std::cout << "Path endpoint indices matched.\n";

    std::cout << "Path sequence lengths:\n\tNaive idx: " << paths.size() << ", Cuttlefish idx: " << kmer_index.sum_paths_len << ".\n";
    if(paths.size() != kmer_index.sum_paths_len)
        return false;

    for(std::size_t i = 0; i < paths.size(); ++i)
        if(paths[i] != DNA_Utility::map_char(DNA::Base(kmer_index.paths[i].get())))
            return false;

    std::cout << "Path sequences aligned.\n";


    // Check minimizer instance counts and their offsets in the paths.

    std::cout << "Unique minimizer count:\n\tNaive idx: " << M.size() << ", Cuttlefish idx: " << kmer_index.min_count << ".\n";
    std::cout << "Minimizer instance count:\n\tNaive idx: " << inst_count << ", Cuttlefish idx: " << kmer_index.num_instances << ".\n";
    if(M.size() != kmer_index.min_count || inst_count != kmer_index.num_instances)
        return false;

    const auto& mi_count = *kmer_index.min_instance_count;
    const auto& m_offset = *kmer_index.min_offset;
    for(const auto p : M)
    {
        const cuttlefish::minimizer_t min = p.first;
        const uint64_t hash = kmer_index.hash(min) - 1;
        if(hash >= M.size())
        {
            std::cout << "Alien minimizer encountered.\n";
            return false;
        }

        const std::size_t count = mi_count[hash] - (hash > 0 ? mi_count[hash - 1] : 0);
        if(count != p.second.size())
        {
            std::cout << "Instance count for minimizer: " << Kmer<l>(p.first) << " with hash " << hash << ":\n\t"
                            "Naive idx: " << p.second.size() << ", Cuttlefish idx: " << count << "\n";
            std::cout << "Instance counts don't match for some minimizers.\n";
            return false;
        }


        const std::size_t min_blk_off = (hash > 0 ? mi_count[hash - 1] : 0);
        for(std::size_t i = 0; i < count; ++i)
            if(m_offset[min_blk_off + i] != p.second[i])
            {
                std::cout << "Differing instance offset for minimizer: " << Kmer<l>(p.first) << " with hash " << hash << ":\n\t"
                                "Naive idx: " << p.second[i] << ", Cuttlefish idx: " << m_offset[min_blk_off + i] << "\n";
                std::cout << "Instance offsets don't match for some minimizers.\n";
                return false;
            }
    }

    std::cout << "Instance count and offsets of individual minimizers matched.\n";


    return true;
}
