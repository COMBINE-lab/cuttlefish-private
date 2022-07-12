
#ifndef KMER_INDEX_HPP
#define KMER_INDEX_HPP



#include "Kmer.hpp"
#include "DNA_Utility.hpp"
#include "Spin_Lock.hpp"
#include "Minimizer_Iterator.hpp"
#include "Minimizer_Instance.hpp"
#include "globals.hpp"
#include "compact_vector/compact_vector.hpp"
#include "BBHash/BooPHF.h"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <thread>


// =============================================================================


// A class to index the k-mers of provided de Bruijn graph path-sequences,
// potentially from many producer threads, based on the k-mers' minimizers.
template <uint16_t k>
class Kmer_Index
{
    typedef cuttlefish::minimizer_t minimizer_t;    // Minimizers can be represented using 64-bit integers.
    typedef compact::vector<uint8_t, 2> path_vector_t;      // Type of the bitvector storing the path sequences.
    typedef compact::vector<std::size_t> path_end_vector_t; // Type of the bitvector storing the path endpoints.
    typedef compact::ts_vector<std::size_t> min_vector_t;   // Type of the bitvectors storing the minimizer instance-counts and -indices.
    // typedef std::vector<char> path_vector_t; // For testing.

    path_vector_t paths;    // The concatenated path-sequences.

    std::vector<std::size_t> path_ends_vec; // Vector for the ending indices, possibly with many unused bits, of the paths in the concatenated sequence.
    path_end_vector_t* path_ends;   // The ending indices of the paths in the concatenated sequence.

    const uint16_t l;   // Size of the l-minimizers.    // TODO: consider templatizing.

    const uint16_t producer_count;  // Number of producer threads supplying the paths to the indexer.

    std::size_t path_count; // Number of paths deposited to the sequence.
    std::size_t sum_paths_len;  // Sum length of all the path-sequences.
    uint64_t num_instances; // Number of minimizer instances in the paths.
    uint64_t min_count; // Number of unique minimizers in the paths.
    uint64_t max_inst_count;    // Maximum count of instances for some minimizer.

    std::vector<path_vector_t> producer_path_buf;   // Separate buffer for each producer, to contain their deposited paths.
    std::vector<std::vector<std::size_t>> producer_path_end_buf;    // Separate buffer for each producer, to contain the (exclusive) indices of the path endpoints in their deposited paths.

    std::vector<std::vector<Minimizer_Instance>> producer_minimizer_buf; // Separate buffer for each producer, to contain their minimizers and their offsets into the deposited paths.

    constexpr static std::size_t buf_sz_th = 5 * 1024 * 1024;   // Threshold for the total size (in bytes) of the buffers per producer: 5 MB.

    // TODO: replace with a file-manager.
    std::vector<std::ofstream> producer_minimizer_file; // Separate file for each producer, to store their minimizers' information.

    std::vector<Minimizer_Instance*> min_group; // Separate buffer to read in each minimizer file.
    std::vector<std::size_t> min_group_size;    // Size of each minimizer file (in element count).

    typedef boomphf::SingleHashFunctor<minimizer_t> minimizer_hasher_t; // The seeded hasher class for minimizers.
                                                                        // TODO: placeholder for now; think it through.
    typedef boomphf::mphf<minimizer_t, minimizer_hasher_t, false> minimizer_mphf_t; // The minimizer-MPHF type.
    constexpr static double gamma = 2.0;    // The gamma parameter of the BBHash algorithm.
    minimizer_mphf_t* min_mphf; // MPHF of the minimizers.

    min_vector_t* min_instance_count;   // Count of instances per each unique minimizer.

    min_vector_t* min_offset;   // Offsets of the instances for the unique minimizers, laid flat all together.

    std::size_t curr_token; // Number of tokens generated for the producers so far.

    Spin_Lock lock; // Mutually-exclusive access lock for different producers.
    constexpr static std::size_t idx_lock_count = 65536;    // Number of locks in the sparse-locks used in various steps.

    std::vector<std::thread> worker;    // Worker threads.


    // Returns the path to the minimizer-information file of producer with ID
    // `producer_id`.
    static const std::string minimizer_file_path(uint16_t producer_id);

    // Dumps the data from `container` to the stream `output`, clearing
    // `container`.
    template <typename T_container_>
    static void dump(std::vector<T_container_>& container, std::ofstream& output);

    // Flushes the buffers of the producer with ID `producer_id`.
    void flush(std::size_t producer_id);

    // Reads in the minimizer instances from each separate minimizer file to
    // memory, and sorts them in parallel.
    void read_and_sort_minimizers();

    // Merges the minimizer instances of each separate producer.
    void merge_minimizers();

    // Constructs the MPHF `min_mphf` over the unique minimizers found.
    void construct_minimizer_mphf();

    // Counts the number of instances for each distinct minimizer.
    void count_minimizer_instances();

    // Enumerates the offsets of the instances for each distinct minimizer.
    void get_minimizer_offsets();

    // Looks up the minimizer `min` in the MPHF and returns its hash value + 1.
    uint64_t hash(minimizer_t min) const;

    // Closes the deposit stream incoming from the producers and flushes the
    // remaining content to disk.
    void close_deposit_stream();

    // Tries to align the k-mer `kmer` to the concatenated paths sequence at the
    // index `idx`. Returns `true` iff the alignment succeeds.
    bool align(const Kmer<k>& kmer, std::size_t idx) const;

    // Binary searches for the maximum rightmost value in the container
    // `container` within the index range `[left, right]` (both ends inclusive)
    // that is at most as `val`. If such a value exists, returns its index.
    // Otherwise, returns `left - 1`.
    template <typename T_container_>
    static int64_t lower_bound(const T_container_& container, int64_t left, int64_t right, typename T_container_::value_type val);

    // Binary searches for the minimum leftmost value in the container
    // `container` within the index range `[left, right]` (both ends inclusive)
    // that is larger than `val`. If such a value exists, returns its index.
    // Otherwise, returns `right + 1`.
    template <typename T_container_>
    static int64_t upper_bound(const T_container_& container, int64_t left, int64_t right, typename T_container_::value_type val);

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

    // Indexes the deposited path-sequences. Should only be invoked after all the
    // sequences to be indexed have been deposited.
    void index();
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

    std::size_t get_id() const { return id; }
};


template <uint16_t k>
inline void Kmer_Index<k>::deposit(const Producer_Token& token, const char* const seq, const std::size_t len)
{
    if(len < k)
        return;

    const std::size_t id = token.get_id();          // The producer's id.
    auto& path_buf = producer_path_buf[id];         // The producer's concatenated-paths sequence buffer.
    auto& path_end_buf = producer_path_end_buf[id]; // The producer's path-endpoints buffer.
    auto& min_buf = producer_minimizer_buf[id];     // The producer's minimizer-information buffer.

    // Get the minimizers and the path sequence.

    Minimizer_Iterator minimizer_iterator(seq, len, k, l); // To iterate over each minimizer in `seq`.
    minimizer_t minimizer;                                 // The minimizer itself.
    std::size_t min_idx;                                   // Index of the minimizer within `seq`.
    std::size_t last_min_idx = len;                        // To track minimizer shifts.

    const std::size_t rel_offset = path_buf.size(); // Relative offset of each minimizer in this producer's path buffer.
    for(std::size_t i = 0; i + (k - 1) < len; ++i)
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

    // Get the ending index (exclusive) of this path in the concatenated paths.
    path_end_buf.emplace_back(path_buf.size());

    const std::size_t buf_size = ((path_buf.size() * 2) / 8) + (path_end_buf.size() * sizeof(std::size_t)) + (min_buf.size() * sizeof(Minimizer_Instance)); // Total buffer size (in bytes) of the producer.
    if(buf_size >= buf_sz_th)
        flush(id);
}


template <uint16_t k>
inline void Kmer_Index<k>::flush(const std::size_t producer_id)
{
    auto& path_buf = producer_path_buf[producer_id];            // The producer's concatenated paths sequence buffer.
    auto& path_end_buf = producer_path_end_buf[producer_id];    // The producer's path-endpoints buffer.
    auto& min_buf = producer_minimizer_buf[producer_id];        // The producer's minimizer-information buffer.

    // Dump the producer-specific path buffer and the endpoints buffer to the global concatenated-paths and -endpoints.
    lock.lock();

    const std::size_t offset_shift = paths.size();

    // TODO: think of faster copying, possibly with chunks.
    for(auto p = path_buf.cbegin(); p != path_buf.cend(); ++p)
        paths.push_back(*p);

    // Shift the producer-specific relative indices of the path-endpoints to their absolute indices into the concatenated paths.
    std::for_each(path_end_buf.begin(), path_end_buf.end(),
                    [offset_shift](auto& p) { p += offset_shift; }
                );

    path_ends_vec.reserve(path_ends_vec.size() + path_end_buf.size());
    path_ends_vec.insert(path_ends_vec.end(), path_end_buf.cbegin(), path_end_buf.cend());

    num_instances += min_buf.size();

    lock.unlock();

    path_buf.clear();
    path_end_buf.clear();


    // Shift the producer-specific relative indices of the minimizers to their absolute indices into the concatenated paths.
    std::for_each(min_buf.begin(), min_buf.end(),
                    [offset_shift](auto& p){ p.shift(offset_shift); }
                );

    // Dump the producer-specific minimizer information to disk.
    auto& minimizer_file = producer_minimizer_file[producer_id];
    dump(min_buf, minimizer_file);
}


template <uint16_t k>
template <typename T_container_>
inline void Kmer_Index<k>::dump(std::vector<T_container_>& container, std::ofstream& output)
{
    if(!output.write(reinterpret_cast<const char*>(container.data()), container.size() * sizeof(T_container_)))
    {
        std::cerr << "Error writing to file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    container.clear();
}


template <uint16_t k>
inline uint64_t Kmer_Index<k>::hash(const cuttlefish::minimizer_t min) const
{
    return min_mphf->lookup(min) + 1;
}


template <uint16_t k>
template <typename T_container_>
inline int64_t Kmer_Index<k>::lower_bound(const T_container_& container, int64_t left, int64_t right, const typename T_container_::value_type val)
{
    int64_t mid, result = left - 1;

    while(left <= right)
    {
        mid = (left + right) >> 1;
        if(container[mid] > val)
            right = mid - 1;
        else
        {
            result = mid;
            left = mid + 1;
        }
    }

    return result;
}


template <uint16_t k>
template <typename T_container_>
inline int64_t Kmer_Index<k>::upper_bound(const T_container_& container, int64_t left, int64_t right, const typename T_container_::value_type val)
{
    int64_t mid, result = right + 1;

    while(left <= right)
    {
        mid = (left + right) >> 1;
        if(container[mid] <= val)
            left = mid + 1;
        else
        {
            result = mid;
            right = mid - 1;
        }
    }

    return result;
}


template <uint16_t k>
__attribute__((optimize("unroll-loops")))
inline bool Kmer_Index<k>::align(const Kmer<k>& kmer, const std::size_t idx) const
{
    constexpr uint16_t word_count = Kmer<k>::num_words();
    const uint64_t* const kmer_data = kmer.data();

    // Align the completely-packed words, i.e. except for possibly the highest-indexed word.
    for(uint16_t word_num = 1; word_num < word_count; ++word_num)
        if(paths.get_int<uint64_t, 32>((idx + k) - word_num * k) != kmer_data[word_num - 1])
            return false;

    // Align the (only) partially-packed word.
    if constexpr(k & 31)
        if(paths.get_int<uint64_t, k & 31>(idx) != kmer_data[word_count - 1])
            return false;

    return true;
}



#endif
