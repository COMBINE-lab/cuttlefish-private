
#ifndef KMER_INDEX_HPP
#define KMER_INDEX_HPP



#include "Kmer.hpp"
#include "Spin_Lock.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>


// =============================================================================


template <uint16_t k>
class Kmer_Index
{
    typedef cuttlefish::minimizer_t minimizer_t;

    const uint16_t worker_count;    // Number of worker threads producing the paths.
    std::vector<std::size_t> worker_string_size;    // Total size of strings (paths) deposited by each worker to the indexer, that has not been flushed yet.

    class Minimizer_Offset_Pair;
    std::vector<std::vector<Minimizer_Offset_Pair>> worker_minimizer_buf;   // Separate buffer for each worker, to contain its minimizers and their offsets into the deposited paths.

    std::size_t curr_token; // Number of tokens produced for the workers so far.

    Spin_Lock lock; // Mutually-exclusive access lock for different workers.


public:

    // Constructs a k-mer indexer that will index sequences produced from at
    // most `worker_count` workers.
    Kmer_Index(uint16_t worker_count);

    class Worker_Token;

    // Returns a unique token object.
    const Worker_Token get_token();
};


// A `(minimizer, offset)` pair, where `minimizer` is the l-minimizer of some
// k-mer `x` in some sequence `seq`, where `minimizer` is at the index `offset`
// in `seq`.
template <uint16_t k>
class Kmer_Index<k>::Minimizer_Offset_Pair
{
private:

    typedef cuttlefish::minimizer_t minimizer_t;

    minimizer_t minimizer;
    std::size_t offset;
};


// A token class to provide distinct token objects to different workers, to
// distinguish them for the indexer.
template <uint16_t k>
class Kmer_Index<k>::Worker_Token
{
    friend class Kmer_Index;

private:

    std::size_t id;

    Worker_Token(const std::size_t id): id(id) {}

    std::size_t get_id() const
    { return id; }
};



#endif
