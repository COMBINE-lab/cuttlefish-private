
#ifndef PARSER_HPP
#define PARSER_HPP



#include "RabbitFX/io/FastxChunk.h"

#include <cstddef>
#include <string>
#include <atomic>


namespace cuttlefish
{

// A class to parse FASTX files.
class Parser
{

private:

    typedef rabbit::fq::FastqDataChunk chunk_t; // Type of chunks containing sequences.
    typedef rabbit::fq::FastqDataPool fq_chunk_pool_t;  // Type of memory pool for chunks.
    typedef rabbit::core::TDataQueue<chunk_t> fq_chunk_queue_t; // Type of queue of chunks ready for consumption.

    std::string file_path;  // File to parse.

    std::size_t consumer_count; // Number of concurrent consumers of the read sequence data.


    // Proof-of-concept production method of parsed sequences.
    void produce(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q);

    // Proof-of-concept consumption method for parsed sequences.
    void consume(fq_chunk_pool_t& chunk_pool, fq_chunk_queue_t& chunk_q, std::vector<std::atomic_uint64_t>& count);

public:

    // Constructs a parser for the file at path `file_path` for `parser_count`
    // parsers.
    Parser(const std::string& file_path, std::size_t parser_count = 1);

    // Proof-of-concept parse method.
    void parse();
};

}



#endif
