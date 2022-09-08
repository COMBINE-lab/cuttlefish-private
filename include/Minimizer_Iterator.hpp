
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
// amortized O(1) time.
template <typename T_seq_>
class Minimizer_Iterator : Minimizer_Utility
{
    typedef cuttlefish::minimizer_t minimizer_t;

private:

    T_seq_ const seq;   // The sequence on which to iterate over.
    const std::size_t seq_len;  // Length of the given sequence.

    const uint16_t k;   // Size of the k-mers.
    const uint16_t l;   // Size of the minimizers.

    minimizer_t last_lmer;  // The last l-mer processed.
    std::size_t last_lmer_idx;  // Index into the sequence of the last l-mer processed.
    const minimizer_t clear_MSN_mask;   // Bitmask to clear the most-significant nucleotide bits of l-mers.

    // Collection of l-mers that have already been seen, and cannot be ruled out
    // yet as candidate minimizers for the k-mers yet to be seen fully.
    std::deque<Lmer_Tuple> dq;

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


template <typename T_seq_>
inline Minimizer_Iterator<T_seq_>::Minimizer_Iterator(T_seq_ const seq, const std::size_t seq_len, const uint16_t k, const uint16_t l):
    seq(seq),
    seq_len(seq_len),
    k(k),
    l(l),
    clear_MSN_mask(~(static_cast<minimizer_t>(0b11) << (2 * (l - 1))))
{
    assert(l <= k);

    last_lmer = 0;
    last_lmer_idx = 0;

    for(std::size_t idx = 0; idx < l; ++idx)
        last_lmer |= (static_cast<minimizer_t>(DNA_Utility::map_base(seq[idx])) << (2 * (l - 1 - idx)));

    dq.emplace_back(last_lmer, last_lmer_idx, hash(last_lmer));


    while(last_lmer_idx + (l - 1) < static_cast<std::size_t>(k - 1))
        operator++();
}


template <typename T_seq_>
inline bool Minimizer_Iterator<T_seq_>::operator++()
{
    if(last_lmer_idx + l == seq_len)
        return false;

    last_lmer_idx++;
    last_lmer = ((last_lmer & clear_MSN_mask) << 2) | DNA_Utility::map_base(seq[last_lmer_idx + l - 1]);


    if(last_lmer_idx + l - 1 >= k)  // The window is now passing over k-length substrings, and slid by offset 1 to the right.
    {
        // End of the current k-mer, e = last_lmer_idx + l - 1; So start of the current k-mer: e - (k - 1).
        const std::size_t curr_kmer_idx = last_lmer_idx + l - 1 - (k - 1);
        if(dq.front().index < curr_kmer_idx)    // This candidate l-mer falls out of the current k-mer.
            dq.pop_front();
    }


    const Lmer_Tuple last_lmer_tuple(last_lmer, last_lmer_idx, hash(last_lmer));
    while(!dq.empty())
        if(dq.back() < last_lmer_tuple)
            break;
        else
            dq.pop_back();

    dq.emplace_back(last_lmer_tuple);

    return true;
}


template <typename T_seq_>
inline void Minimizer_Iterator<T_seq_>::value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index) const
{
    minimizer = dq.front().lmer;
    index = dq.front().index;
}



#endif
