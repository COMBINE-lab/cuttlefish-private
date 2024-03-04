
#ifndef MINIMIZER_ITERATOR_HPP
#define MINIMIZER_ITERATOR_HPP



#include "Minimizer_Utility.hpp"
#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "globals.hpp"
#include "utility.hpp"
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

    T_seq_ seq; // The sequence on which to iterate over.
    std::size_t seq_len;    // Length of the given sequence.

    const uint16_t k;   // Size of the k-mers.
    const uint16_t l;   // Size of the minimizers.

    const uint64_t seed;    // Seed for hashing l-mers.

    minimizer_t last_lmer;  // The last l-mer processed.
    minimizer_t last_lmer_bar;  // Reverse complement of the last l-mer processed.
    std::size_t last_lmer_idx;  // Index into the sequence of the last l-mer processed.
    const minimizer_t clear_MSN_mask;   // Bitmask to clear the most-significant nucleotide bits of l-mers.

    // typedef std::deque<Lmer_Tuple> deque_t;
    typedef deque<Lmer_Tuple> deque_t;

    // Collection of l-mers that have already been seen, and cannot be ruled out
    // yet as candidate minimizers for the k-mers yet to be seen fully. `dq_f`
    // is for the forward-strand and `dq_r` is for the reverse-strand.
    deque_t dq_f, dq_r;

public:

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of a given sequence in a streaming manner. The seed-value
    // `seed` is used in hashing the `l`-mers.
    Minimizer_Iterator(uint16_t k, uint16_t l, uint64_t seed = 0);

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of the sequence `seq`, of length `seq_len`, in a streaming
    // manner. The iterator sits at the first k-mer after the construction.
    // The seed-value `seed` is used in hashing the `l`-mers.
    Minimizer_Iterator(T_seq_ seq, std::size_t seq_len, uint16_t k, uint16_t l, uint64_t seed = 0);

    // Resets the iterator to the sequence `seq` of length `seq_len`. The
    // iterator sits at the first k-mer after the construction.
    void reset(T_seq_ seq, std::size_t seq_len);

    // Moves the iterator to a next k-mer by extending it with the character
    // `ch`.
    void advance(char ch);

    // Moves the iterator to the next k-mer in the sequence. Returns `true` iff
    // the current k-mer is not the last k-mer in the sequence.
    bool operator++();

    // Puts the minimizer of the current k-mer into `minimizer`, and stores its
    // index in the sequence at `index`. The position information over the
    // sequence is maintained internally.
    void value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index) const;
};


template <typename T_seq_, bool is_canonical_>
inline Minimizer_Iterator<T_seq_, is_canonical_>::Minimizer_Iterator(const uint16_t k, const uint16_t l, const uint64_t seed):
      k(k)
    , l(l)
    , seed(seed)
    , clear_MSN_mask(~(static_cast<minimizer_t>(0b11) << (2 * (l - 1))))
    , dq_f(2 * k - l)
    , dq_r(2 * k - l)
{
    assert(l <= k);
}


template <typename T_seq_, bool is_canonical_>
inline Minimizer_Iterator<T_seq_, is_canonical_>::Minimizer_Iterator(T_seq_ const seq, const std::size_t seq_len, const uint16_t k, const uint16_t l, const uint64_t seed):
    Minimizer_Iterator(k, l, seed)
{
    reset(seq, seq_len);
}


template <typename T_seq_, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, is_canonical_>::reset(T_seq_ const seq, const std::size_t seq_len)
{
    this->seq = seq;
    this->seq_len = seq_len;
    dq_f.clear(), dq_r.clear();

    last_lmer = 0;
    last_lmer_bar = 0;
    last_lmer_idx = 0;

    DNA::Base base;
    for(std::size_t idx = 0; idx < l; ++idx)
    {
        base = DNA_Utility::map_base(seq[idx]);
        assert(base != DNA::Base::N);
        last_lmer |= (static_cast<minimizer_t>(base) << (2 * (l - 1 - idx)));

        if constexpr(is_canonical_)
            last_lmer_bar |= (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * idx));
    }

    dq_f.emplace_back(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer, seed));
    if constexpr(is_canonical_)
        dq_r.emplace_back(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar, seed));

    while(last_lmer_idx + (l - 1) < static_cast<std::size_t>(k - 1))
        operator++();
}


template <typename T_seq_, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, is_canonical_>::advance(const char ch)
{
    assert(DNA_Utility::is_DNA_base(ch));

    last_lmer_idx++;
    const DNA::Base base = DNA_Utility::map_base(ch);
    last_lmer = ((last_lmer & clear_MSN_mask) << 2) | base;
    if constexpr(is_canonical_)
        last_lmer_bar = (last_lmer_bar >> 2) | (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * (l - 1)));


    // End of the current k-mer, e = last_lmer_idx + l - 1; So start of the current k-mer: e - (k - 1).
    const std::size_t curr_kmer_idx = (last_lmer_idx + l - 1 >= k - 1lu ? last_lmer_idx + l - 1 - (k - 1) : 0);
    if(dq_f.front().index < curr_kmer_idx)    // This candidate l-mer falls out of the current k-mer.
        dq_f.pop_front();

    if constexpr(is_canonical_)
        if(dq_r.front().index < curr_kmer_idx)  // This candidate l-mer falls out of the current k-mer.
            dq_r.pop_front();


    const auto fix_dq =
        [](deque_t& dq, const Lmer_Tuple& lmer_tuple)
        {
            while(!dq.empty())
                if(dq.back() < lmer_tuple)
                    break;
                else
                    dq.pop_back();

            dq.push_back(lmer_tuple);
        };


    fix_dq(dq_f, Lmer_Tuple(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer, seed)));
    if constexpr(is_canonical_)
        fix_dq(dq_r, Lmer_Tuple(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar, seed)));
}


template <typename T_seq_, bool is_canonical_>
inline bool Minimizer_Iterator<T_seq_, is_canonical_>::operator++()
{
    if(last_lmer_idx + l == seq_len)
        return false;

    advance(seq[last_lmer_idx + l]);
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
