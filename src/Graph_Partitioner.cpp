
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


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Graph_Partitioner<k, Colored_>::Graph_Partitioner(Subgraphs_Manager<k, Colored_>& subgraphs, const Data_Logistics& logistics, const uint16_t l):
      subgraphs(subgraphs)
    , seqs(logistics.input_paths_collection())
    , l_(l)
    , subgraphs_path_pref(logistics.subgraphs_path())
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

    std::cerr << "Number of records: " << record_count_ << ".\n";
    std::cerr << "Number of super k-mers: " << super_kmer_count_ << ".\n";
    std::cerr << "Total length of the super k-mers: " << super_kmers_len_ << ".\n";
}


template <uint16_t k, bool Colored_>
void Graph_Partitioner<k, Colored_>::parse(chunk_pool_t& chunk_pool, chunk_q_t& chunk_q)
{
    rabbit::int64 chunk_count = 0;

    for(const auto& file_path : seqs)
    {
        // TODO: address memory-reuse issue within a reader for every new instance.
        rabbit::fq::FastqFileReader reader(file_path, chunk_pool);
        const chunk_t* chunk;
        while((chunk = reader.readNextChunk()) != NULL)
            chunk_q.Push(chunk_count++, chunk);
    }

    chunk_q.SetCompleted();
    std::cerr << "Parsed " << chunk_count << " chunks in total from " << seqs.size() << " files.\n";
}


template <uint16_t k, bool Colored_>
void Graph_Partitioner<k, Colored_>::process(chunk_q_t& chunk_q, chunk_pool_t& chunk_pool)
{
    rabbit::int64 chunk_id = 0;
    std::vector<neoReference> parsed_chunk;
    chunk_t* chunk;

    uint64_t rec_count = 0;
    uint64_t super_kmer_count = 0;
    uint64_t super_kmers_len = 0;

    Minimizer_Iterator<const char*, true> min_it(k, l_, min_seed);
    while(chunk_q.Pop(chunk_id, chunk))
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
                std::size_t frag_beg;   // Beginning index of the next fragment.
                std::size_t frag_len;   // Length of the next fragment.

                // Skip placeholder bases.
                for(frag_beg = last_frag_end; frag_beg + k <= seq_len; ++frag_beg)
                    if(DNA_Utility::is_DNA_base(seq[frag_beg]))
                        break;

                if(frag_beg + k > seq_len)  // No more sequence fragment remains with complete k-mers.
                    break;

                // Check whether the first k-mer of the fragment has any placeholder bases.
                for(frag_len = 1; frag_len < k; ++frag_len)
                    if(!DNA_Utility::is_DNA_base(seq[frag_beg + frag_len]))
                        break;

                if(frag_len < k)
                {
                    last_frag_end = frag_beg + frag_len;
                    continue;
                }


                minimizer_t last_min, last_min_idx;
                minimizer_t min, min_idx;
                std::size_t curr_sup_kmer_off = 0; // Relative offset of the current super k-mer in the fragment.

                min_it.reset(seq + frag_beg, seq_len - frag_beg);
                min_it.value_at(last_min, last_min_idx);
                frag_len = k;

                while(DNA_Utility::is_DNA_base(seq[frag_beg + frag_len]))
                {
                    min_it.advance(seq[frag_beg + frag_len]);
                    min_it.value_at(min, min_idx);
                    frag_len++;
                    if(min_idx != last_min_idx)
                    {
                        // TODO: add assert for minimizer-diff validity.
                        const auto next_sup_kmer_off = frag_len - k;
                        super_kmers_len += (next_sup_kmer_off - curr_sup_kmer_off) + k - 1;
                        curr_sup_kmer_off = next_sup_kmer_off;
                        super_kmer_count++;

                        last_min = min, last_min_idx = min_idx;
                    }
                }

                super_kmer_count++;
                super_kmers_len += frag_len - curr_sup_kmer_off;
                last_frag_end = frag_beg + frag_len;
            }
        }


        chunk_pool.Release(chunk);
    }


    record_count_ += rec_count;
    super_kmer_count_ += super_kmer_count;
    super_kmers_len_ += super_kmers_len;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Graph_Partitioner)
