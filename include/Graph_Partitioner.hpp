
#ifndef SUPER_KMER_PARTITIONER_HPP
#define SUPER_KMER_PARTITIONER_HPP



#include "Data_Logistics.hpp"
#include "RabbitFX/io/FastxChunk.h"
#include "RabbitFX/io/DataQueue.h"

#include <cstdint>
#include <atomic>
#include <string>
#include <vector>


namespace cuttlefish
{

// =============================================================================
// Class to partition de Bruijn graphs into subgraphs based on minimizers of the
// `k`-mers. Effectively, splits input sequences to maximal weak super k-mers
// and distributes those to appropriate subgraphs.
// TODO: only supports FASTQ for now; extend to FASTA also.
template <uint16_t k>
class Graph_Partitioner
{

private:

    typedef rabbit::fq::FastqDataChunk chunk_t; // Type of chunks containing parsed sequences.
    typedef rabbit::fq::FastqDataPool chunk_pool_t; // Type of memory pools for chunks.
    typedef rabbit::core::TDataQueue<chunk_t> chunk_q_t;    // Type of queue of parsed chunks.

    const std::vector<std::string> seqs;    // Input sequence collection.

    const uint16_t l_;  // Size of minimizers for the super k-mers.
    static constexpr uint64_t min_seed = 0; // Seed for `l`-minimizer hashing.

    const std::string subgraphs_path_pref;  // Path prefix for the subgraphs' super k-mer buckets.

    std::atomic_uint64_t record_count_; // Number of records in the sequences.
    std::atomic_uint64_t super_kmer_count_; // Number of super k-mers in the sequences.
    std::atomic_uint64_t super_kmers_len_;  // Total length of the super k-mers in the sequences.


    // Parses the provided sequences into chunks from the memory pool
    // `chunk_pool` and puts the parsed chunks into the queue `chunk_q`.
    void parse(chunk_pool_t& chunk_pool, chunk_q_t& chunk_q);

    // Processes the parsed chunks from the queue `chunk_q` and returns the
    // processed chunks to the memory pool `chunk_pool`.
    void process(chunk_q_t& chunk_q, chunk_pool_t& chunk_pool);

public:

    // Constructs a de Bruijn graph partitioner with `l`-minimizers for the
    // sequences from the data logistics manager `logistics`.
    Graph_Partitioner(const Data_Logistics& logistics, uint16_t l);

    // Returns the size of minimizers for the super k-mers.
    auto l() const { return l_; }

    // Returns the number of records in the sequences.
    uint64_t record_count() const { return record_count_; }

    // Returns the number of super k-mers in the sequences.
    uint64_t super_kmer_count() const { return super_kmer_count_; }

    // Returns the total length of the super k-mers in the sequences.
    uint64_t super_kmers_len() const { return super_kmers_len_; }

    // Partitions the passed sequences into maximal weak super k-mers and
    // deposits those to the super k-mer buckets.
    void partition();
};

}



#endif
