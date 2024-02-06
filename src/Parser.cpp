
#include "Parser.hpp"
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
    for(uint64_t t = 0; t < consumer_count; ++t)
        consumer.emplace_back(&Parser::consume, this, std::ref(fq_chunk_pool), std::ref(fq_chunk_q));

    producer.join();
    std::for_each(consumer.begin(), consumer.end(), [](auto& t){ t.join(); });
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


void Parser::consume(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q)
{
    std::size_t line_sum = 0;
    rabbit::int64 id = 0;
    std::vector<neoReference> data;
    chunk_t* fq_chunk;

    data.resize(100'000);
    while(chunk_q.Pop(id, fq_chunk))
    {
        line_sum += rabbit::fq::chunkFormat(fq_chunk, data);
        chunk_pool.Release(fq_chunk);
    }

    std::cerr << "Line count: " << line_sum << ".\n";
}

}
