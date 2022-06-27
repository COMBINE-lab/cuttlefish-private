
#ifndef KMER_INDEX_HPP
#define KMER_INDEX_HPP



#include "Kmer.hpp"
#include "DNA_Utility.hpp"
#include "Spin_Lock.hpp"
#include "Minimizer_Iterator.hpp"
#include "globals.hpp"
#include "compact_vector/compact_vector.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <fstream>
#include <cstdlib>


// =============================================================================


// A class to index the k-mers of provided sequences, potentially from many
// producer threads, based on the k-mers' minimizers.
template <uint16_t k>
class Kmer_Index
{
    typedef cuttlefish::minimizer_t minimizer_t;
    typedef compact::vector<uint8_t, 2> bitvector_t;
    // typedef std::vector<char> bitvector_t;  // For testing.

    bitvector_t paths;  // The concatenated paths sequence.

    const uint16_t l;   // Size of the l-minimizers.    // TODO: consider templatizing.

    const uint16_t producer_count;  // Number of producer threads supplying the paths to the indexer.

    std::vector<bitvector_t> producer_path_buf; // Separate buffer for each producer, to contain their deposited paths.

    class Minimizer_Offset_Pair;
    std::vector<std::vector<Minimizer_Offset_Pair>> producer_minimizer_buf; // Separate buffer for each producer, to contain their minimizers and their offsets into the deposited paths.

    std::vector<std::ofstream> producer_minimizer_file; // Separate file for each producer, to store their minimizers' information.

    constexpr static std::size_t buf_sz_th = 5 * 1024 * 1024;   // Threshold for the total size (in bytes) of the buffers per producer: 5 MB.

    std::size_t curr_token; // Number of tokens produced for the producers so far.

    Spin_Lock lock; // Mutually-exclusive access lock for different producers.


    // Flushes the buffers of the producer with ID `producer_id`.
    void flush(std::size_t producer_id);

public:

    // Constructs a k-mer indexer that will index sequences produced from at
    // most `producer_count` producers, based on the k-mers' `l`-minimizers.
    Kmer_Index(uint16_t l, uint16_t producer_count);

    class Producer_Token;

    // Returns a unique token object.
    const Producer_Token get_token();

    // Deposits the sequence `seq` of length `len` to the index, from a producer
    // having the token `token`.
    void deposit(const Producer_Token &token, const char* seq, std::size_t len);

    // Flushes the remaining content from the producers.
    void finalize_production();
};


// A `(minimizer, offset)` pair, where `minimizer` is the l-minimizer (for a
// pre-defined `l`) of some k-mer `x` in some sequence `seq`, where `minimizer`
// is at the index `offset` in `seq`.
template <uint16_t k>
class Kmer_Index<k>::Minimizer_Offset_Pair
{
private:

    typedef cuttlefish::minimizer_t minimizer_t;

    minimizer_t minimizer;  // The minimizer.
    std::size_t offset;     // Offset of the minimizer in the underlying sequence.

public:

    // Constructs a tuple with the minimizer `minimizer` and its offset `offset`.
    Minimizer_Offset_Pair(const minimizer_t minimizer, const std::size_t offset):
        minimizer(minimizer),
        offset(offset)
    {}

    // Shifts (to the right) the offset of the minimizer by `offset_shit`.
    void shift(const std::size_t offset_shift)  { offset += offset_shift; }
};


// A token class to provide distinct token objects to different producers, to
// distinguish them for the indexer.
template <uint16_t k>
class Kmer_Index<k>::Producer_Token
{
    friend class Kmer_Index;

private:

    std::size_t id;

    Producer_Token(const std::size_t id): id(id) {}

    std::size_t get_id() const
    { return id; }
};


template <uint16_t k>
inline void Kmer_Index<k>::deposit(const Producer_Token& token, const char* const seq, const std::size_t len)
{
    const std::size_t id = token.get_id();      // The producer's id.
    auto& path_buf = producer_path_buf[id];     // The producer's concatenated-paths sequence buffer.
    auto& min_buf = producer_minimizer_buf[id]; // The producer's minimizer-information buffer.

    // Get the minimizers and the path sequence.

    Minimizer_Iterator minimizer_iterator(seq, len, k, l); // To iterate over each minimizer in `seq`.
    minimizer_t minimizer;                                 // The minimizer itself.
    std::size_t min_idx;                                   // Index of the minimizer within `seq`.
    std::size_t last_min_idx = len;                        // To track minimizer shifts.

    const std::size_t rel_offset = path_buf.size(); // Relative offset of each minimizer in this producer's path buffer.
    for (std::size_t i = 0; i + (k - 1) < len; ++i)
    {
        minimizer_iterator.value_at(minimizer, min_idx);
        if(min_idx != last_min_idx)
            min_buf.emplace_back(minimizer, rel_offset + min_idx),
            last_min_idx = min_idx;

        path_buf.push_back(DNA_Utility::map_base(seq[i]));
        // path_buf.push_back(seq[i]); // For testing.

        ++minimizer_iterator;
    }

    // Get the (k - 1)-length tail of the last k-mer.
    for(std::size_t i = len - (k - 1); i < len; ++i)
        path_buf.push_back(DNA_Utility::map_base(seq[i]));
        // path_buf.push_back(seq[i]); // For testing

    const std::size_t buf_size = ((path_buf.size() * 2) / 8) + (min_buf.size() * sizeof(Minimizer_Offset_Pair)); // Total buffer size (in bytes) of the producer.
    if(buf_size >= buf_sz_th)
        flush(id);
}


template <uint16_t k>
inline void Kmer_Index<k>::flush(const std::size_t producer_id)
{
    auto& path_buf = producer_path_buf[producer_id];    // The producer's concatenated paths sequence buffer.

    // Dump the producer-specific path buffer to the global concatenated-paths sequence.
    lock.lock();

    const std::size_t offset_shift = paths.size();

    // TODO: think of faster copying, possibly with chunks.
    for(auto p = path_buf.cbegin(); p != path_buf.cend(); ++p)
        paths.push_back(*p);

    lock.unlock();

    path_buf.clear();


    // Shift the producer-specific relative indices of the minimizers to their absolute indices into the concatenated paths.
    auto& min_buf = producer_minimizer_buf[producer_id];
    std::for_each(
                    min_buf.begin(), min_buf.end(),
                    [offset_shift](auto& p){ p.shift(offset_shift); }
                );

    // Dump the producer-specific minimizer information to disk.
    auto& minimizer_file = producer_minimizer_file[producer_id];
    minimizer_file.write(reinterpret_cast<const char*>(min_buf.data()), min_buf.size() * sizeof(Minimizer_Offset_Pair));
    if(!minimizer_file)
    {
        std::cerr << "Error writing to the minimizer files. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    min_buf.clear();
}



#endif
