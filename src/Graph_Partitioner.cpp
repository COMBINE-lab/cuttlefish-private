
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
    , weak_super_kmer_count_(0)
    , weak_super_kmers_len_(0)
    , super_km1_mers_len_(0)
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
    std::cerr << "Number of super (k - 1)-mers: " << weak_super_kmer_count_ << ".\n";
    std::cerr << "Total length of the weak super k-mers:  " << weak_super_kmers_len_ << ".\n";
    std::cerr << "Total length of the super (k - 1)-mers: " << super_km1_mers_len_ << ".\n";
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
    uint64_t sup_km1_mer_count = 0; // Number of parsed super (k - 1)-mers.
    std::size_t sup_km1_mers_len = 0;   // Total length of the super (k - 1)-mers, in bases.
    std::size_t weak_sup_kmers_len = 0; // Total length of the (weak) super k-mers, in bases.

    typedef Minimizer_Iterator<const char*, true> min_it_t;
    min_it_t min_it(k - 1, l_, min_seed);   // `l`-minimizer iterator for `(k - 1)`-mers.
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
                std::size_t frag_len;   // Length of the fragment.

                // Skip placeholder bases.
                for(frag_beg = last_frag_end; frag_beg + (k + 1) <= seq_len; ++frag_beg)
                    if(DNA_Utility::is_DNA_base(seq[frag_beg]))
                        break;

                if(frag_beg + (k + 1) > seq_len)    // No more fragment remains with complete (k + 1)-mers, i.e. edges.
                    break;

                const auto frag = seq + frag_beg;   // Current fragment.

                // Check whether the first (k + 1)-mer of the fragment has any placeholder bases.
                for(frag_len = 1; frag_len < (k + 1); ++frag_len)
                    if(!DNA_Utility::is_DNA_base(frag[frag_len]))
                        break;

                if(frag_len < (k + 1))
                {
                    last_frag_end = frag_beg + frag_len;
                    continue;
                }


                minimizer_t last_min;   // Last minimizer found in the iteration.
                minimizer_t min;    // Current minimizer in the iteration.
                std::size_t last_min_off, min_off;  // Relative offset of the minimizers in the fragment.
                std::size_t cur_sup_km1_mer_off = 0;    // Relative offset of the current super (k - 1)-mer in the fragment.
                std::size_t km1_mer_idx = 0;    // Index of the current (k - 1)-mer in the current super (k - 1)-mer.
                frag_len = k - 1;

                min_it.reset(frag, seq_len - frag_beg); // The fragment length is an estimate; upper-bound to be exact.
                min_it.value_at(min, min_off);
                last_min = min, last_min_off = min_off;

                while(DNA_Utility::is_DNA_base(frag[frag_len]))
                {
                    min_it.advance(frag[frag_len]);
                    km1_mer_idx++, frag_len++;

                    min_it.value_at(min, min_off);
                    if(min_off != last_min_off) // Encountered a discontinuity vertex—the k-mer whose suffix is the current (k - 1)-mer.
                    {
                        // Either the last minimizer just fell out of the (k - 1)-mer window or the new minimizer sits at the last l-mer.
                        assert( last_min_off == cur_sup_km1_mer_off + km1_mer_idx - 1 ||
                                min_off == cur_sup_km1_mer_off + km1_mer_idx + (k - 1) - l_);

                        // The `(k - 1)`-mers of this discontinuity k-mer have different minimizer instances.
                        assert( min_it_t::minimizer_idx(frag + cur_sup_km1_mer_off + km1_mer_idx - 1, k - 1, l_, min_seed) !=
                                1 + min_it_t::minimizer_idx(frag + cur_sup_km1_mer_off + km1_mer_idx, k - 1, l_, min_seed));

                        const auto next_sup_km1_mer_off = cur_sup_km1_mer_off + km1_mer_idx;
                        assert(next_sup_km1_mer_off == frag_len - (k - 1));
                        const auto len = km1_mer_idx + (k - 2); // Length of the super (k - 1)-mer that just terminated.
                        sup_km1_mers_len += len;
                        sup_km1_mer_count++;

                        const bool l_disc = (cur_sup_km1_mer_off > 0);
                        const bool r_disc = true;
                        subgraphs.add_super_kmer(last_min, frag + cur_sup_km1_mer_off - l_disc, l_disc + len + r_disc, l_disc, r_disc);
                        weak_sup_kmers_len += l_disc + len + r_disc;

                        cur_sup_km1_mer_off = next_sup_km1_mer_off;
                        last_min = min, last_min_off = min_off;
                        km1_mer_idx = 0;
                    }
                }

                const auto len = frag_len - cur_sup_km1_mer_off;
                sup_km1_mers_len += len;
                sup_km1_mer_count++;

                const bool l_disc = (cur_sup_km1_mer_off > 0);
                const bool r_disc = false;
                subgraphs.add_super_kmer(last_min, frag + cur_sup_km1_mer_off - l_disc, l_disc + len + r_disc, l_disc, r_disc);
                weak_sup_kmers_len += l_disc + len + r_disc;

                last_frag_end = frag_beg + frag_len;
            }
        }


        chunk_pool.Release(chunk);
        chunk_count++;
    }


    chunk_count_ += chunk_count;
    record_count_ += rec_count;
    weak_super_kmer_count_ += sup_km1_mer_count;
    weak_super_kmers_len_ += weak_sup_kmers_len;
    super_km1_mers_len_ += sup_km1_mers_len;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Graph_Partitioner)
