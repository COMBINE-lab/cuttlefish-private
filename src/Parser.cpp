
#include "Parser.hpp"
#include "DNA_Utility.hpp"
#include "Minimizer_Iterator.hpp"
#include "RabbitFX/io/FastxStream.h"
#include "RabbitFX/io/Formater.h"
#include "RabbitFX/io/Globals.h"

#include <cstdint>
#include <vector>
#include <iostream>
#include <functional>
#include <thread>


namespace cuttlefish
{


Parser::Parser(const std::string& file_path, const std::size_t parser_count):
      file_path(file_path)
    , consumer_count(parser_count)
{}


void Parser::parse()
{
    fq_chunk_pool_t fq_chunk_pool;  // Memory pool for chunks of sequences.
    fq_chunk_queue_t fq_chunk_q;    // Read chunks ready for parse.

    std::thread producer(&Parser::produce, this, std::ref(fq_chunk_pool), std::ref(fq_chunk_q));

    std::vector<std::thread> consumer;
    // std::vector<std::atomic_uint64_t> count(5);
    std::atomic_uint64_t count(0);
    for(uint64_t t = 0; t < consumer_count; ++t)
        consumer.emplace_back(&Parser::consume_split_super_kmers, this, std::ref(fq_chunk_pool), std::ref(fq_chunk_q), std::ref(count));

    producer.join();
    std::for_each(consumer.begin(), consumer.end(), [](auto& t){ t.join(); });

    std::cerr << "Count of super k-mers: " << count << ".\n";
}


void Parser::produce(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q)
{
    rabbit::fq::FastqFileReader fq_reader(file_path, chunk_pool);
    rabbit::int64 chunk_count = 0;

    while(true)
    {
        auto fq_chunk = fq_reader.readNextChunk();
        if(fq_chunk == NULL)
            break;

        chunk_q.Push(chunk_count++, fq_chunk);
    }

    chunk_q.SetCompleted();
    std::cerr << "File " << file_path << " has " << chunk_count << " chunks.\n";
}


void Parser::consume_count_bases(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q, std::vector<std::atomic_uint64_t>& count)
{
    rabbit::int64 id = 0;
    std::vector<neoReference> parsed_chunk;
    chunk_t* fq_chunk;

    uint64_t rec_count = 0;
    uint64_t nuc_count[] = {0, 0, 0, 0, 0};

    while(chunk_q.Pop(id, fq_chunk))
    {
        parsed_chunk.clear();
        rec_count += rabbit::fq::chunkFormat(fq_chunk, parsed_chunk);

        for(const auto& rec : parsed_chunk)
        {
            auto* const seq = reinterpret_cast<const char*>(rec.base + rec.pseq);
            const auto len = rec.lseq;
            for(std::size_t i = 0; i < len; ++i)
            {
                const auto b = DNA_Utility::to_upper(seq[i]);
                nuc_count[DNA_Utility::map_base(b)]++;
            }
        }

        chunk_pool.Release(fq_chunk);
    }


    for(std::size_t b = 0; b < 5; ++b)
        count[b] += nuc_count[b];
}


void Parser::consume_split_super_kmers(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q, std::atomic_uint64_t& count)
{
    rabbit::int64 id = 0;
    std::vector<neoReference> parsed_chunk;
    chunk_t* fq_chunk;

    constexpr uint16_t k = 31;
    constexpr uint16_t l = 18;
    constexpr uint64_t min_seed = 0;

    uint64_t rec_count = 0;
    uint64_t sup_kmer_count = 0;

    while(chunk_q.Pop(id, fq_chunk))
    {
        parsed_chunk.clear();
        rec_count += rabbit::fq::chunkFormat(fq_chunk, parsed_chunk);

        for(const auto& rec : parsed_chunk)
        {
            auto* const seq = reinterpret_cast<const char*>(rec.base + rec.pseq);
            const auto seq_len = rec.lseq;

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

                // Check whether the first k-mer has any placeholder bases.
                for(frag_len = 1; frag_len < k; ++frag_len)
                    if(!DNA_Utility::is_DNA_base(seq[frag_beg + frag_len]))
                        break;

                if(frag_len < k)
                {
                    last_frag_end = frag_beg + frag_len;
                    continue;
                }


                Minimizer_Iterator<decltype(seq), true, false> min_it(seq + frag_beg, seq_len - frag_beg, k - 1, l, min_seed);
                minimizer_t last_min, last_min_idx;
                minimizer_t min, min_idx;
                min_it.value_at(last_min, last_min_idx);
                frag_len = k - 1;
                while(++min_it)
                {
                    frag_len++;
                    min_it.value_at(min, min_idx);
                    if(min_idx != last_min_idx)
                    {
                        sup_kmer_count++;
                        last_min = min, last_min_idx = min_idx;
                    }
                }

                sup_kmer_count++;

                last_frag_end = frag_beg + frag_len;
            }
        }

        chunk_pool.Release(fq_chunk);
    }


    count += sup_kmer_count;
}


}
