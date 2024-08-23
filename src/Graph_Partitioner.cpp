
#include "Graph_Partitioner.hpp"
#include "Subgraphs_Manager.hpp"
#include "Minimizer_Iterator.hpp"
#include "DNA_Utility.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "RabbitFX/io/Reference.h"
#include "RabbitFX/io/Formater.h"
#include "RabbitFX/io/Globals.h"
#include "parlay/parallel.h"

#include <algorithm>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k, bool Is_FASTQ_, bool Colored_>
Graph_Partitioner<k, Is_FASTQ_, Colored_>::Graph_Partitioner(Subgraphs_Manager<k, Colored_>& subgraphs, const Data_Logistics& logistics, const uint16_t l):
      subgraphs(subgraphs)
    , seqs(logistics.input_paths_collection())
    , l_(l)
    , sup_km1_mer_len_th(2 * (k - 1) - l_)
    , chunk_pool_sz(parlay::num_workers() * (Is_FASTQ_ ? 2 : 4))    // TODO: make a more informed choice.
    , chunk_pool(chunk_pool_sz)
    , chunk_q(chunk_pool_sz)
    , subgraphs_path_pref(logistics.subgraphs_path())
    , stat_w(parlay::num_workers())
{}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::partition()
{
    parlay::par_do(
        [&]()
        {
            read_chunks();
        },
        [&]()
        {
            const auto process_chunks = [&](std::size_t){ process(); };
            parlay::parallel_for(0, parlay::num_workers(), process_chunks, 1);
        });


    Worker_Stats stat;
    std::for_each(stat_w.cbegin(), stat_w.cend(), [&](const auto& s){ stat += s.unwrap(); });
    std::cerr << "Number of processed chunks: " << stat.chunk_count << ".\n";
    std::cerr << "Total size of chunks: " << stat.chunk_bytes << ".\n";
    std::cerr << "Number of records: " << stat.record_count << ".\n";
    std::cerr << "Number of super (k - 1)-mers: " << stat.weak_super_kmer_count << ".\n";
    std::cerr << "Total length of the weak super k-mers:  " << stat.weak_super_kmers_len << ".\n";
    std::cerr << "Total length of the super (k - 1)-mers: " << stat.super_km1_mers_len << ".\n";
    std::cerr << "Total work in parse: " << stat.parse_time << "s.\n";
    std::cerr << "Total work in processing records: " << stat.process_time << "s.\n";
    std::cerr << "Max work in processing records: " <<
        [&](){ double t = 0; std::for_each(stat_w.cbegin(), stat_w.cend(), [&](const auto& s){ t = std::max(t, s.unwrap().process_time); }); return t; }() << "s.\n";
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::read_chunks()
{
    const auto t_s = timer::now();
    uint64_t chunk_count = 0;   // Number of read chunks.
    rabbit::int64 source_id = 1;    // Source (i.e. file) ID of a chunk.

    for(const auto& file_path : seqs)
    {
        // TODO: address memory-reuse issue within a reader for every new instance.
        typename RabbitFX_DS_type<Is_FASTQ_>::reader_t reader(file_path, chunk_pool);
        const auto push_chunk = [&](auto const chunk)
        {
            if(chunk == NULL)
                return false;

            chunk_q.Push(source_id, chunk);
            chunk_count++;
            return true;
        };

        bool chunks_remain = true;
        while(chunks_remain)
            if constexpr(Is_FASTQ_)
                chunks_remain = push_chunk(reader.readNextChunk());
            else
                chunks_remain = push_chunk(reader.readNextChunkList());

        source_id++;
    }

    chunk_q.SetCompleted();
    const auto t_e = timer::now();
    std::cerr << "Read " << chunk_count << " chunks in total from " << seqs.size() << " files. Time elapsed: " << timer::duration(t_e - t_s) << " s.\n";
}


// TODO: remove.
template <uint16_t k, bool Is_FASTQ_, bool Colored_>
uint64_t Graph_Partitioner<k, Is_FASTQ_, Colored_>::read_chunks(const std::size_t source_id)
{
    assert(Colored_);

    uint64_t chunk_count = 0;   // Number of read chunks.

    // TODO: address memory-reuse issue within a reader for every new instance.
    typename RabbitFX_DS_type<Is_FASTQ_>::reader_t reader(seqs[source_id - 1], chunk_pool);
    const auto push_chunk = [&](auto const chunk)
    {
        if(chunk == NULL)
            return false;

        chunk_q.Push(source_id, chunk);
        chunk_count++;
        return true;
    };

    bool chunks_remain = true;
    while(chunks_remain)
        if constexpr(Is_FASTQ_)
            chunks_remain = push_chunk(reader.readNextChunk());
        else
            chunks_remain = push_chunk(reader.readNextChunkList());

    chunk_q.SetCompleted();

    return chunk_count;
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::process()
{
    rabbit::int64 source_id;    // Source (i.e. file) ID of a chunk.
    chunk_t* chunk; // Current chunk.
    std::vector<parsed_seq_t> parsed_chunk; // Current parsed chunk.
    rabbit::int64 last_source = 0;  // Source ID of the last chunk processed.
    (void)last_source;

    while(chunk_q.Pop(source_id, chunk))
    {
        assert(source_id >= last_source);

        process_chunk(chunk, source_id, parsed_chunk);
        last_source = source_id;
    }
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::process_chunk(chunk_t* chunk, const uint32_t source_id, std::vector<parsed_seq_t>& parsed_chunk)
{
    auto& w_stat = stat_w[parlay::worker_id()].unwrap();

    uint64_t sup_km1_mer_count = 0; // Number of parsed super (k - 1)-mers.
    std::size_t sup_km1_mers_len = 0;   // Total length of the super (k - 1)-mers, in bases.
    std::size_t weak_sup_kmers_len = 0; // Total length of the (weak) super k-mers, in bases.

    // Minimizer_Iterator<const char*, k - 1, true> min_it(l_, min_seed);   // `l`-minimizer iterator for `(k - 1)`-mers.
    Min_Iterator<k - 1> min_it(l_); // `l`-minimizer iterator for `(k - 1)`-mers.

    const auto t_0 = timer::now();
    parsed_chunk.clear();
    if constexpr(Is_FASTQ_)
        w_stat.chunk_bytes += chunk->size,
        w_stat.record_count += rabbit::fq::chunkFormat(chunk, parsed_chunk);
    else
    {
        w_stat.chunk_bytes += chunk->chunk->size,
        w_stat.record_count += rabbit::fa::chunkListFormat(*chunk, parsed_chunk);

        auto ptr = chunk->chunk;
        do
        {
            chunk_pool.Release(ptr);
            ptr = ptr->next;
        }
        while(ptr != NULL);
    }

    const auto t_1 = timer::now();
    w_stat.parse_time += timer::duration(t_1 - t_0);

    for(auto& record : parsed_chunk)
    {
        const char* seq;
        std::size_t seq_len;

        if constexpr(Is_FASTQ_)
            seq = reinterpret_cast<const char*>(record.base + record.pseq),
            seq_len = record.lseq;
        else
        {
            record.seq.resize(record.seq.length() + 32);    // Deal with UB in `Super_Kmer_Chunk::add_encoded_label`
            seq = record.seq.c_str(),
            seq_len = record.length;
        }

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


            // minimizer_t cur_min;    // Minimizer of the current super (k - 1)-mer in the iteration.
            // minimizer_t next_min;   // Minimizer of the next super (k - 1)-mer in the iteration.
            uint64_t cur_h; // 64-bit hash of the current super (k - 1)-mer's minimizer.
            uint64_t next_h;    // 64-bit hash of the next super (k - 1)-mer's minimizer.
            std::size_t prev_g; // Subgraph ID of the previous super (k - 1)-mer's minimizer.
            std::size_t cur_g;  // Subgraph ID of the current super (k - 1)-mer's minimizer.
            std::size_t next_g; // Subgraph ID of the next super (k - 1)-mer's minimizer.
            // std::size_t cur_min_off, next_min_off;  // Relative offsets of the minimizers in the fragment. Used for assertion checks.
            std::size_t cur_sup_km1_mer_off = 0;    // Relative offset of the current super (k - 1)-mer in the fragment.
            std::size_t km1_mer_idx = 0;    // Index of the current (k - 1)-mer in the current super (k - 1)-mer.
            frag_len = k - 1;

            // min_it.reset(frag, seq_len - frag_beg); // The fragment length is an estimate; upper-bound to be exact.
            // min_it.value_at(cur_min, cur_min_off, cur_h);
            min_it.reset(frag);
            cur_h = min_it.hash();
            cur_g = subgraphs.graph_ID(cur_h);
            prev_g = subgraphs.graph_count();   // To deal with false-positive `-Wmaybe-uninitialized` later on.

            while(DNA_Utility::is_DNA_base(frag[frag_len]))
            {
                const auto len = km1_mer_idx + (k - 1); // Length of the current super (k - 1)-mer.

                min_it.advance(frag[frag_len]);
                km1_mer_idx++, frag_len++;

                // min_it.value_at(next_min, next_min_off, next_h);
                next_h = min_it.hash();
                next_g = subgraphs.graph_ID(next_h);
/*                  assert(next_min_off >= cur_sup_km1_mer_off + km1_mer_idx);

                if(next_min_off != cur_min_off)
                    // Either the last minimizer just fell out of the (k - 1)-mer window or the new minimizer sits at the last l-mer.
                    assert( cur_min_off == cur_sup_km1_mer_off + km1_mer_idx - 1 ||
                            next_min_off == cur_sup_km1_mer_off + km1_mer_idx + (k - 1) - l_);
*/
                // Either encountered a discontinuity vertex—the k-mer whose suffix is the current (k - 1)-mer, or the super (k - 1)-mer extends too long.
                if(next_g != cur_g || len == sup_km1_mer_len_th)
                {
                    if(next_g != cur_g)
                        // The `(k - 1)`-mers of this discontinuity k-mer have minimizers mapping to different subgraphs.
                        assert(is_discontinuity(frag + cur_sup_km1_mer_off + km1_mer_idx - 1));

                    const auto next_sup_km1_mer_off = cur_sup_km1_mer_off + km1_mer_idx;
                    assert(next_sup_km1_mer_off == frag_len - (k - 1));
                    sup_km1_mers_len += len;
                    sup_km1_mer_count++;

                    const bool l_joined = (cur_sup_km1_mer_off > 0);    // Whether this weak super k-mer is joined to the one to its left.
                    const bool r_joined = true; // Whether this weak super k-mer is joined to the one to its right.
                    const bool l_disc = (l_joined && prev_g != cur_g);  // Whether it's left-discontinuous.
                    const bool r_disc = (next_g != cur_g);  // Whether it's right-discontinuous.
                    // const bool l_cont = (l_joined && prev_g == cur_g);  // Whether it's left-continuous.
                    // const bool r_cont = (next_g == cur_g);  // Whether it's right-continuous.
                    const auto len_weak = l_joined + len + r_joined;    // Length of the weak super k-mer.
                    assert(len_weak >= k);
                    // TODO: the following add, being to different subgraphs' different worker-buffers, causes lots of cache misses.
                    if constexpr(!Colored_)
                        subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, l_disc, r_disc);
                    else
                        subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, source_id, l_disc, r_disc);
                    weak_sup_kmers_len += len_weak;

                    cur_sup_km1_mer_off = next_sup_km1_mer_off;
                    prev_g = cur_g;
                    cur_g = next_g;
                    km1_mer_idx = 0;
                }

                // cur_min = next_min;
                // cur_min_off = next_min_off; // The minimizer-instance offsets are tracked only for assertion checks.
                cur_h = next_h;
            }

            const auto len = frag_len - cur_sup_km1_mer_off;
            sup_km1_mers_len += len;
            sup_km1_mer_count++;

            const bool l_joined = (cur_sup_km1_mer_off > 0);
            const bool r_joined = false;
            const bool l_disc = (l_joined && prev_g != cur_g);
            const bool r_disc = false;
            // const bool l_cont = (l_joined && prev_g == cur_g);
            // const bool r_cont = false;
            const auto len_weak = l_joined + len + r_joined;
            assert(len_weak >= k);
            // TODO: the following add, being to different subgraphs' different worker-buffers, causes lots of cache misses.
            if constexpr(!Colored_)
                subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, l_disc, r_disc);
            else
                subgraphs.add_super_kmer(cur_g, frag + cur_sup_km1_mer_off - l_joined, len_weak, source_id, l_disc, r_disc);
            weak_sup_kmers_len += len_weak;

            last_frag_end = frag_beg + frag_len;
        }
    }

    const auto t_2 = timer::now();
    w_stat.process_time += timer::duration(t_2 - t_1);

    if constexpr(Is_FASTQ_)
        chunk_pool.Release(chunk);

    w_stat.chunk_count++;
    w_stat.weak_super_kmer_count += sup_km1_mer_count;
    w_stat.weak_super_kmers_len += weak_sup_kmers_len;
    w_stat.super_km1_mers_len += sup_km1_mers_len;
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
void Graph_Partitioner<k, Is_FASTQ_, Colored_>::Worker_Stats::operator+=(const Worker_Stats& rhs)
{
    chunk_count += rhs.chunk_count;
    chunk_bytes += rhs.chunk_bytes;
    record_count += rhs.record_count;
    weak_super_kmer_count += rhs.weak_super_kmer_count;
    weak_super_kmers_len += rhs.weak_super_kmers_len;
    super_km1_mers_len += rhs.super_km1_mers_len;
    parse_time += rhs.parse_time;
    process_time += rhs.process_time;
}


template <uint16_t k, bool Is_FASTQ_, bool Colored_>
bool Graph_Partitioner<k, Is_FASTQ_, Colored_>::is_discontinuity(const char* const seq) const
{
    minimizer_t min_l, min_r;
    uint64_t h_l, h_r;
    std::size_t idx_l, idx_r;

    typedef Minimizer_Iterator<const char*, k - 1, true> min_it_t;
    min_it_t::minimizer(seq, l_, min_seed, min_l, h_l, idx_l);
    min_it_t::minimizer(seq + 1, l_, min_seed, min_r, h_r, idx_r);

    if(subgraphs.graph_ID(h_l) == subgraphs.graph_ID(h_r))
        return false;

    return true;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL_L2, cuttlefish::Graph_Partitioner)
