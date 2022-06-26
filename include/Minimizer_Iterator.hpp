
#ifndef MINIMIZER_ITERATOR_HPP
#define MINIMIZER_ITERATOR_HPP



#include "DNA_Utility.hpp"
#include "globals.hpp"
#include "xxHash/xxh3.h"

#include <cstdint>
#include <queue>
#include <cassert>


// =============================================================================


// A class to iterate over the minimizers of the constituent k-mers of a given
// sequence, computing the minimizer for each k-mer in amortized O(1) time.
class Minimizer_Iterator
{
    typedef cuttlefish::minimizer_t minimizer_t;

    class Lmer_Tuple;

private:

    const char* const seq;  // The sequence on which to iterate over.
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
    Minimizer_Iterator(const char* seq, std::size_t seq_len, uint16_t k, uint16_t l);

    // Moves the iterator to the next k-mer in the sequence. Returns `true` iff
    // the current k-mer is not the last k-mer in the sequence.
    bool operator++();

    // Puts the minimizer of the current k-mer into `minimizer`, and stores its
    // index in the sequence at `index`. The position information over the
    // sequence is maintained internally.
    void value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index) const;

    // Returns the hash value of the l-mer `lmer`.
    static uint64_t hash(cuttlefish::minimizer_t lmer);

};


// A class to pack information regarding l-mers, to aid in computing
// l-minimizers of k-mers over sequences by `Minimizer_Iterator`.
class Minimizer_Iterator::Lmer_Tuple
{
    friend class Minimizer_Iterator;

    typedef cuttlefish::minimizer_t minimizer_t;

private:

    minimizer_t lmer;   // The l-mer.
    std::size_t index;  // Its index in the underlying sequence.
    uint64_t hash;      // Its hash value, determining the l-mer ordering.


    // Returns `true` iff this l-mer tuple is to be ordered as lesser to the
    // tuple `rhs`. The primary ordering is based on the tuples' hashes. If
    // equal, then the ordering is based on their literal form. If also equal,
    // then the tuple occurring earlier in the underlying sequence is lesser.
    bool operator<(const Lmer_Tuple& rhs) const;


public:

    // Constructs a tuple for an l-mer `lmer`, positioned at index `index` of
    // the underlying sequence, and having a hash value `hash`.
    Lmer_Tuple(minimizer_t lmer, std::size_t index, uint64_t hash);
};


inline Minimizer_Iterator::Lmer_Tuple::Lmer_Tuple(const minimizer_t lmer, const std::size_t index, const uint64_t hash):
    lmer(lmer),
    index(index),
    hash(hash)
{}


inline Minimizer_Iterator::Minimizer_Iterator(const char* const seq, const std::size_t seq_len, const uint16_t k, const uint16_t l):
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


inline bool Minimizer_Iterator::operator++()
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


inline void Minimizer_Iterator::value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index) const
{
    minimizer = dq.front().lmer;
    index = dq.front().index;
}


inline uint64_t Minimizer_Iterator::hash(const cuttlefish::minimizer_t lmer)
{
#ifdef CF_DEVELOP_MODE
    return lmer;    // TODO: add as debug option for developer.
#endif

    return XXH3_64bits_withSeed(&lmer, sizeof(lmer), seed);
}


inline bool Minimizer_Iterator::Lmer_Tuple::operator<(const Lmer_Tuple& rhs) const
{
    if(hash < rhs.hash)
        return true;

    if(hash > rhs.hash)
        return false;

    if(lmer < rhs.lmer)
        return true;

    if(lmer > rhs.lmer)
        return false;

    return index < rhs.index;
}



#endif
