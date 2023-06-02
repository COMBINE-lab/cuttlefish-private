
#ifndef MINIMIZER_ITERATOR_HPP
#define MINIMIZER_ITERATOR_HPP



#include "Minimizer_Utility.hpp"
#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "globals.hpp"
#include "xxHash/xxh3.h"

#include <cstdint>
#include <queue>
#include <cassert>


// =============================================================================
// A class to iterate over the minimizers of the constituent k-mers of a given
// sequence of type `T_seq_`, computing the minimizer for each k-mer in
// amortized O(1) time. If `is_canonical_` is `true`, then canonical minimizers
// are computed, i.e. both strand-forms of the sequence is considered.
template <typename T_seq_, bool is_canonical_ = false>
class Minimizer_Iterator
{
    typedef cuttlefish::minimizer_t minimizer_t;

private:

    T_seq_ const seq;   // The sequence on which to iterate over.
    const std::size_t seq_len;  // Length of the given sequence.

    const uint16_t k;   // Size of the k-mers.
    const uint16_t l;   // Size of the minimizers.

    minimizer_t last_lmer;  // The last l-mer processed.
    minimizer_t last_lmer_bar;  // Reverse complement of the last l-mer processed.
    std::size_t last_lmer_idx;  // Index into the sequence of the last l-mer processed.
    const minimizer_t clear_MSN_mask;   // Bitmask to clear the most-significant nucleotide bits of l-mers.

    // Collection of l-mers that have already been seen, and cannot be ruled out
    // yet as candidate minimizers for the k-mers yet to be seen fully. `dq_f`
    // is for the forward-strand and `dq_r` is for the reverse-strand.
    std::deque<Lmer_Tuple> dq_f, dq_r;

    constexpr static uint64_t seed = 0; // Seed for hashing l-mers.

public:

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of the sequence `seq`, of length `seq_len`, in a streaming
    // manner. The iterator sits at the first k-mer after the construction.
    Minimizer_Iterator(T_seq_ seq, std::size_t seq_len, uint16_t k, uint16_t l);

    // Moves the iterator to the next k-mer in the sequence. Returns `true` iff
    // the current k-mer is not the last k-mer in the sequence.
    bool operator++();

    // Puts the minimizer of the current k-mer into `minimizer`, and stores its
    // index in the sequence at `index`. The position information over the
    // sequence is maintained internally.
    void value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index) const;
};


template <typename T_seq_, bool is_canonical_>
inline Minimizer_Iterator<T_seq_, is_canonical_>::Minimizer_Iterator(T_seq_ const seq, const std::size_t seq_len, const uint16_t k, const uint16_t l):
    seq(seq),
    seq_len(seq_len),
    k(k),
    l(l),
    clear_MSN_mask(~(static_cast<minimizer_t>(0b11) << (2 * (l - 1))))
{
    assert(l <= k);

    last_lmer = 0;
    last_lmer_bar = 0;
    last_lmer_idx = 0;

    DNA::Base base;
    for(std::size_t idx = 0; idx < l; ++idx)
    {
        base = DNA_Utility::map_base(seq[idx]);
        last_lmer |= (static_cast<minimizer_t>(base) << (2 * (l - 1 - idx)));

        if constexpr(is_canonical_)
            last_lmer_bar |= (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * idx));
    }

    dq_f.emplace_back(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer));
    if constexpr(is_canonical_)
        dq_r.emplace_back(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar));

    while(last_lmer_idx + (l - 1) < static_cast<std::size_t>(k - 1))
        operator++();
}


template <typename T_seq_, bool is_canonical_>
inline bool Minimizer_Iterator<T_seq_, is_canonical_>::operator++()
{
    if(last_lmer_idx + l == seq_len)
        return false;

    last_lmer_idx++;
    const DNA::Base base = DNA_Utility::map_base(seq[last_lmer_idx + l - 1]);
    last_lmer = ((last_lmer & clear_MSN_mask) << 2) | base;
    if constexpr(is_canonical_)
        last_lmer_bar = (last_lmer_bar >> 2) | (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * (l - 1)));


    if(last_lmer_idx + l - 1 >= k)  // The window is now passing over k-length substrings, and slid by offset 1 to the right.
    {
        // End of the current k-mer, e = last_lmer_idx + l - 1; So start of the current k-mer: e - (k - 1).
        const std::size_t curr_kmer_idx = last_lmer_idx + l - 1 - (k - 1);
        if(dq_f.front().index < curr_kmer_idx)    // This candidate l-mer falls out of the current k-mer.
            dq_f.pop_front();

        if constexpr(is_canonical_)
            if(dq_r.front().index < curr_kmer_idx)  // This candidate l-mer falls out of the current k-mer.
                dq_r.pop_front();
    }


    const auto fix_dq =
        [](std::deque<Lmer_Tuple>& dq, const Lmer_Tuple& lmer_tuple)
        {
            while(!dq.empty())
                if(dq.back() < lmer_tuple)
                    break;
                else
                    dq.pop_back();

            dq.push_back(lmer_tuple);
        };


    fix_dq(dq_f, Lmer_Tuple(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer)));
    if constexpr(is_canonical_)
        fix_dq(dq_r, Lmer_Tuple(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar)));

    return true;
}


template <typename T_seq_, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, is_canonical_>::value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index) const
{
    if constexpr(is_canonical_)
    {
        if(dq_f.front() < dq_r.front())
            minimizer = dq_f.front().lmer, index = dq_f.front().index;
        else
            minimizer = dq_r.front().lmer, index = dq_r.front().index;
    }
    else
        minimizer = dq_f.front().lmer, index = dq_f.front().index;
}



#endif
