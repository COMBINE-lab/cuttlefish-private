
#include "Parser.hpp"
#include "DNA_Utility.hpp"
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
    std::vector<std::atomic_uint64_t> count(5);
    for(uint64_t t = 0; t < consumer_count; ++t)
        consumer.emplace_back(&Parser::consume, this, std::ref(fq_chunk_pool), std::ref(fq_chunk_q), std::ref(count));

    producer.join();
    std::for_each(consumer.begin(), consumer.end(), [](auto& t){ t.join(); });

    for(std::size_t b = 0; b < 4; ++b)
        std::cerr << DNA_Utility::map_char(DNA::Base(b)) << " : " << count[b] << "\n";
    std::cerr << "N : " << count[4] << "\n";
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


void Parser::consume(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q, std::vector<std::atomic_uint64_t>& count)
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

}
