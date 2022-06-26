
#ifndef KMER_INDEX_HPP
#define KMER_INDEX_HPP



#include "Kmer.hpp"
#include "DNA_Utility.hpp"
#include "Spin_Lock.hpp"
#include "globals.hpp"
#include "compact_vector/compact_vector.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>


// =============================================================================


template <uint16_t k>
class Kmer_Index
{
    typedef cuttlefish::minimizer_t minimizer_t;
    typedef compact::vector<uint8_t, 2> bitvector_t;

    bitvector_t paths;  // The concatenated paths sequence.

    const uint16_t worker_count;    // Number of worker threads producing the paths.

    std::vector<bitvector_t> worker_path_buf;   // Separate buffer for each worker, to contain their deposited paths.

    class Minimizer_Offset_Pair;

    constexpr static std::size_t buf_sz_th = 5 * 1024 * 1024;   // Threshold for the total size (in bytes) of the buffers per worker: 5 MB.

    std::size_t curr_token; // Number of tokens produced for the workers so far.

    Spin_Lock lock; // Mutually-exclusive access lock for different workers.


    // Flushes the buffers of the worker with ID `worker_id`.
    void flush(std::size_t worker_id);

public:

    // Constructs a k-mer indexer that will index sequences produced from at
    // most `worker_count` workers.
    Kmer_Index(uint16_t worker_count);

    class Worker_Token;

    // Returns a unique token object.
    const Worker_Token get_token();

    // Deposits the sequence `seq` of length `len` to the index, from a worker
    // having the token `token`.
    void deposit(const Worker_Token& token, const char* seq, std::size_t len);
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


template <uint16_t k>
inline void Kmer_Index<k>::deposit(const Worker_Token& token, const char* const seq, const std::size_t len)
{
    const std::size_t id = token.get_id();  // The worker's id.
    auto& path_buf = worker_path_buf[id];   // The worker's concatenated-paths sequence buffer.


    // Get the path sequence.
    for(std::size_t i = 0; i < len; ++i)
        path_buf.push_back(DNA_Utility::map_base(seq[i]));


    const std::size_t buf_size = ((path_buf.size() * 2) / 8);   // Total buffer size (in bytes) of the worker.
    if(buf_size >= buf_sz_th)
        flush(id);
}


template <uint16_t k>
inline void Kmer_Index<k>::flush(const std::size_t worker_id)
{
    auto& path_buf = worker_path_buf[worker_id];    // The worker's concatenated-paths sequence buffer.

    // Dump the worker-specific path buffer to the global concatenated-paths sequence.
    lock.lock();

    // TODO: think of faster copying, possibly with chunks.
    for(auto p = path_buf.cbegin(); p != path_buf.cend(); ++p)
        paths.push_back(*p);

    lock.unlock();

    path_buf.clear();
}



#endif
