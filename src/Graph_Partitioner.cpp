
#include "Graph_Partitioner.hpp"
#include "Subgraphs_Manager.hpp"
#include "Minimizer_Iterator.hpp"
#include "DNA_Utility.hpp"
#include "globals.hpp"
#include "RabbitFX/io/FastxStream.h"
#include "RabbitFX/io/Reference.h"
#include "RabbitFX/io/Formater.h"
#include "parlay/parallel.h"

#include <cmath>
#include <functional>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Graph_Partitioner<k, Colored_>::Graph_Partitioner(Subgraphs_Manager<k, Colored_>& subgraphs, const Data_Logistics& logistics, const uint16_t l):
      subgraphs(subgraphs)
    , seqs(logistics.input_paths_collection())
    , l_(l)
    , subgraphs_path_pref(logistics.subgraphs_path())
    , chunk_count_(0)
    , record_count_(0)
    , super_kmer_count_(0)
    , super_kmers_len_(0)
{}


template <uint16_t k, bool Colored_>
void Graph_Partitioner<k, Colored_>::partition()
{
    const auto chunk_count = std::ceil(parlay::num_workers() * 1.1);    // Maximum number of chunks. TODO: make a more informed choice.
    chunk_pool_t chunk_pool(chunk_count);   // Memory pool for chunks of sequences.
    chunk_q_t chunk_q(chunk_count); // Parsed chunks.

    std::thread parser(&Graph_Partitioner::parse, this, std::ref(chunk_pool), std::ref(chunk_q));

    const auto process = [&](std::size_t){ this->process(chunk_q, chunk_pool); };
    parlay::parallel_for(0, parlay::num_workers(), process, 1);

    parser.join();

    std::cerr << "Number of processed chunks: " << chunk_count_ << ".\n";
    std::cerr << "Number of records: " << record_count_ << ".\n";
    std::cerr << "Number of super k-mers: " << super_kmer_count_ << ".\n";
    std::cerr << "Total length of the super k-mers: " << super_kmers_len_ << ".\n";
}


template <uint16_t k, bool Colored_>
void Graph_Partitioner<k, Colored_>::parse(chunk_pool_t& chunk_pool, chunk_q_t& chunk_q)
{
    uint64_t chunk_count = 0;   // Number of parsed chunks.
    rabbit::int64 source_id = 0;    // Source (i.e. file) ID of a chunk.

    for(const auto& file_path : seqs)
    {
        // TODO: address memory-reuse issue within a reader for every new instance.
        rabbit::fq::FastqFileReader reader(file_path, chunk_pool);
        const chunk_t* chunk;
        while((chunk = reader.readNextChunk()) != NULL)
        {
            chunk_q.Push(source_id, chunk);
            chunk_count++;
        }

        source_id++;
    }

    chunk_q.SetCompleted();
    std::cerr << "Parsed " << chunk_count << " chunks in total from " << seqs.size() << " files.\n";
}


template <uint16_t k, bool Colored_>
void Graph_Partitioner<k, Colored_>::process(chunk_q_t& chunk_q, chunk_pool_t& chunk_pool)
{
    rabbit::int64 source_id;    // Source (i.e. file) ID of a chunk.
    uint64_t chunk_count = 0;   // Number of processed chunks.
    chunk_t* chunk; // Current chunk.
    std::vector<neoReference> parsed_chunk; // Current parsed chunk.

    uint64_t rec_count = 0; // Number of parsed records.
    uint64_t super_kmer_count = 0;  // Number of parsed super k-mers.
    uint64_t super_kmers_len = 0;   // Total length of the super k-mers, in bases.

    Minimizer_Iterator<const char*, true> min_it(k, l_, min_seed);
    while(chunk_q.Pop(source_id, chunk))
    {
        parsed_chunk.clear();
        rec_count += rabbit::fq::chunkFormat(chunk, parsed_chunk);

        for(const auto& record : parsed_chunk)
        {
            const auto seq = reinterpret_cast<const char*>(record.base + record.pseq);
            const auto seq_len = record.lseq;

            std::size_t last_frag_end = 0;  // Ending index (exclusive) of the last sequence fragment.
            while(true)
            {
                std::size_t frag_beg;   // Offset of the next fragment in the sequence.
                std::size_t frag_len;   // Length of the next fragment.

                // Skip placeholder bases.
                for(frag_beg = last_frag_end; frag_beg + k <= seq_len; ++frag_beg)
                    if(DNA_Utility::is_DNA_base(seq[frag_beg]))
                        break;

                if(frag_beg + k > seq_len)  // No more sequence fragment remains with complete k-mers.
                    break;

                const auto frag = seq + frag_beg;   // Current fragment.

                // Check whether the first k-mer of the fragment has any placeholder bases.
                for(frag_len = 1; frag_len < k; ++frag_len)
                    if(!DNA_Utility::is_DNA_base(frag[frag_len]))
                        break;

                if(frag_len < k)
                {
                    last_frag_end = frag_beg + frag_len;
                    continue;
                }


                minimizer_t last_min, last_min_idx;
                minimizer_t min, min_idx;
                std::size_t curr_sup_kmer_off = 0; // Relative offset of the current super k-mer in the fragment.
                std::size_t kmer_idx = 0;   // Index of the current k-mer in the current super k-mer.

                min_it.reset(frag, seq_len - frag_beg); // The fragment length is an estimate; upper-bound to be exact.
                min_it.value_at(last_min, last_min_idx);
                frag_len = k;

                while(DNA_Utility::is_DNA_base(frag[frag_len]))
                {
                    min_it.advance(frag[frag_len]);
                    kmer_idx++, frag_len++;

                    min_it.value_at(min, min_idx);
                    if(min_idx != last_min_idx)
                    {
                        // Either the last minimizer fell out of the k-window or the new minimizer sits at the last l-mer.
                        assert(last_min_idx < curr_sup_kmer_off + kmer_idx || min_idx == curr_sup_kmer_off + kmer_idx + k - l_);

                        const auto next_sup_kmer_off = frag_len - k;
                        const auto len = (next_sup_kmer_off - curr_sup_kmer_off) + (k - 1);
                        super_kmers_len += len;
                        curr_sup_kmer_off = next_sup_kmer_off;
                        super_kmer_count++;

                        last_min = min, last_min_idx = min_idx;
                        kmer_idx = 0;
                    }
                }

                super_kmer_count++;
                super_kmers_len += frag_len - curr_sup_kmer_off;
                last_frag_end = frag_beg + frag_len;
            }
        }


        chunk_pool.Release(chunk);
        chunk_count++;
    }


    chunk_count_ += chunk_count;
    record_count_ += rec_count;
    super_kmer_count_ += super_kmer_count;
    super_kmers_len_ += super_kmers_len;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Graph_Partitioner)
