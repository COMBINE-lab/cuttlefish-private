
#include "Unitig_Collator.hpp"
#include "Unitig_File.hpp"
#include "Kmer.hpp"
#include "dBG_Utilities.hpp"
#include "Character_Buffer.hpp"
#include "FASTA_Record.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"
#include "xxHash/xxh3.h"

#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cassert>


// TODO: use more efficient data structures throughout the collator.


namespace cuttlefish
{

template <uint16_t k>
Unitig_Collator<k>::Unitig_Collator(P_e_t& P_e, const Data_Logistics& logistics, op_buf_list_t& op_buf, const std::size_t gmtig_bucket_count):
      P_e(P_e)
    , lmtig_buckets_path(logistics.lmtig_buckets_path())
    , unitig_coord_buckets_path(logistics.unitig_coord_buckets_path())
    , max_unitig_bucket_count(gmtig_bucket_count)
    , op_buf(op_buf)
{
     // TODO: fix better policy?
    if((max_unitig_bucket_count & (max_unitig_bucket_count - 1)) != 0)
    {
        std::cerr << "Maximal unitig bucket count needs to be a power of 2. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    max_bucket_sz = 0;
    std::size_t sum_bucket_sz = 0;
    std::for_each(P_e.cbegin(), P_e.cend(),
        [&](const auto& bucket){ max_bucket_sz = std::max(max_bucket_sz, bucket.data().size()); sum_bucket_sz += bucket.data().size(); });

    std::cerr << "Sum edge-bucket size: " << sum_bucket_sz << "\n";
    std::cerr << "Maximum edge-bucket size: " << max_bucket_sz << "\n";
}


template <uint16_t k>
void Unitig_Collator<k>::collate()
{
    const auto t_0 = timer::now();
    map();

    const auto t_1 = timer::now();
    std::cerr << "Time taken in mapping: " << timer::duration(t_1 - t_0) << "s.\n";

    reduce();
    const auto t_2 = timer::now();
    std::cerr << "Time taken in reduction: " << timer::duration(t_2 - t_1) << "s.\n";

    // TODO: print meta-information over the maximal unitigs'.
}


template <uint16_t k>
void Unitig_Collator<k>::map()
{
    typedef Buffer<Path_Info<k>> map_t;
    typedef Buffer<unitig_path_info_t> buf_t;

    std::vector<Padded_Data<map_t>> M_vec(parlay::num_workers());   // Worker-local DATs of edge-index to edge path-info.
    std::vector<Padded_Data<buf_t>> buf_vec(parlay::num_workers()); // Worker-local buffers to read in edge path-info.

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            M_vec[w_id].data().resize(max_bucket_sz);   // TODO: thread-local allocation suits best here.
        }, 1);


    max_unitig_bucket.reserve(max_unitig_bucket_count);
    for(std::size_t i = 0; i < max_unitig_bucket_count; ++i)
        max_unitig_bucket.emplace_back(unitig_coord_buckets_path + "_" + std::to_string(i));

    std::atomic_uint64_t edge_c = 0;    // Number of edges (i.e. unitigs) found.
#ifndef NDEBUG
    std::atomic_uint64_t h_p_e = 0; // Hash of the edges' path-information.
#endif

    // Maps the edges (i.e. unitigs) from bucket `b` to maximal-unitig buckets.
    const auto map_to_max_unitig_bucket =
    [&](const std::size_t b)
    {
        const auto w_id = parlay::worker_id();
        auto const M = M_vec[w_id].data().data();
        auto& buf = buf_vec[w_id].data();

        const auto b_sz = load_path_info(b, M, buf);
        edge_c += b_sz;
        P_e[b].data().remove();

#ifndef NDEBUG
        uint64_t h = 0;
        for(std::size_t i = 0; i < b_sz; ++i)
            h ^= M[i].hash();

        h_p_e ^= h;
#endif


        const auto bucket_path = lmtig_buckets_path + "_" + std::to_string(b);
        assert(file_exists(bucket_path));
        Unitig_File_Reader unitig_reader(bucket_path);
        Buffer<char> unitig;    // Read-off unitig.
        uni_idx_t idx = 0;  // The unitig's sequential ID in the bucket.
        std::size_t uni_len;    // The unitig's length in bases.
        for(; (uni_len = unitig_reader.read_next_unitig(unitig)) > 0; idx++)
        {
            assert(idx < b_sz);
            const auto p = M[idx].p();  // Path-ID of this unitig.
            const auto mapped_b_id = XXH3_64bits(&p, sizeof(p)) & (max_unitig_bucket_count - 1); // Hashed maximal unitig bucket.
            // TODO: use `wyhash`.

            max_unitig_bucket[mapped_b_id].data().add(M[idx], unitig.data(), uni_len);
        }

        assert(idx == b_sz);
        unitig_reader.remove_files();
    };

    parlay::parallel_for(1, P_e.size(), map_to_max_unitig_bucket, 1);

    std::cerr << "Found " << edge_c << " edges.\n";
#ifndef NDEBUG
    std::cerr << "Edges' path-information signature: " << h_p_e << "\n";
#endif
}


template <uint16_t k>
void Unitig_Collator<k>::reduce()
{
    std::size_t max_max_uni_b_sz = 0;   // Maximum unitig-count in some maximal unitig bucket.
    std::size_t max_max_uni_b_label_len = 0;    // Maximum dump-string length in some maximal unitig bucket.
    std::for_each(max_unitig_bucket.cbegin(), max_unitig_bucket.cend(),
        [&](const auto& b)
        {
            max_max_uni_b_sz = std::max(max_max_uni_b_sz, b.data().size()),
            max_max_uni_b_label_len = std::max(max_max_uni_b_label_len, b.data().label_len());
        });

    std::cerr << "Maximum maximal unitig bucket size:  " << max_max_uni_b_sz << "\n";
    std::cerr << "Maximum maximal unitig label length: " << max_max_uni_b_label_len << "\n";


    typedef Buffer<Unitig_Coord<k>> coord_buf_t;
    typedef Buffer<char> label_buf_t;
    std::vector<Padded_Data<coord_buf_t>> U_vec(parlay::num_workers()); // Worker-local buffers for unitig coordinate information.
    std::vector<Padded_Data<label_buf_t>> L_vec(parlay::num_workers()); // Worker-local buffers for dump-strings in buckets.

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            // TODO: thread-local allocations here suit best.
            U_vec[w_id].data().resize(max_max_uni_b_sz);
            L_vec[w_id].data().resize(max_max_uni_b_label_len);
        }, 1);

    // TODO: add per-worker progress tracker.

    const auto collate_max_unitig_bucket =
    [&](const std::size_t b)
    {
        const auto w_id = parlay::worker_id();
        auto const U = U_vec[w_id].data().data();   // Coordinate info of the unitigs.
        auto const L = L_vec[w_id].data().data();   // Dump-strings of the unitig labels.
        auto& output = op_buf[w_id].data(); // Output buffer for the maximal unitigs.

        const auto b_sz = max_unitig_bucket[b].data().load_coords(U);
        const auto len = max_unitig_bucket[b].data().load_labels(L);
        max_unitig_bucket[b].data().remove();
        (void)len;

        std::sort(U, U + b_sz);


        std::string max_unitig, max_u_rc;
        std::size_t i, j;
        for(i = 0; i < b_sz; i = j)
        {
            max_unitig.clear();

            const bool is_cycle = U[i].is_cycle();
            const std::size_t s = i;
            std::size_t e = s + 1;

            // Find the current maximal unitig's stretch in the bucket.
            while(e < b_sz && U[e].p() == U[s].p())
            {
                assert(U[e].o() != side_t::unspecified); assert(U[e].is_cycle() == is_cycle);
                e++;
            }

            if(e - s == 2 && !is_cycle) // Special case for handling orientation of path-traversals in the discontinuity graph.
            {
                assert(U[s].r() == 0); assert(U[s + 1].r() == 0);

                std::string u0(L + U[s].label_idx(), U[s].label_len());
                std::string u1(L + U[s + 1].label_idx(), U[s + 1].label_len());

                if(U[s].o() == side_t::front)
                    reverse_complement(u0);

                if(U[s + 1].o() == side_t::front)
                    reverse_complement(u1);

                reverse_complement(u1);
                max_unitig = u0 + std::string(u1.begin() + k, u1.end());
            }
            else
                for(j = s; j < e; ++j)
                {
                    std::string u_label(L + U[j].label_idx(), U[j].label_len());
                    if(U[j].o() == side_t::front)
                        reverse_complement(u_label);

                    max_unitig += (max_unitig.empty() ? u_label : std::string(u_label.begin() + k, u_label.end()));
                }


            if(is_cycle)
                max_unitig.pop_back();  // Cyclic maximal unitig traversals start and end at the same vertex, so one copy needs to be removed.


            max_u_rc = max_unitig;
            reverse_complement(max_u_rc);
            if(is_cycle)
            {
                Kmer<k> min_f, min_r;
                const auto min_idx_f = min_kmer(max_unitig, min_f);
                const auto min_idx_r = min_kmer(max_u_rc, min_r);
                if(min_f < min_r)
                    max_unitig = std::string(max_unitig.cbegin() + min_idx_f, max_unitig.cend()) + std::string(max_unitig.cbegin() + k - 1, max_unitig.cbegin() + min_idx_f + k - 1);
                else
                    max_unitig = std::string(max_u_rc.cbegin() + min_idx_r, max_u_rc.cend()) + std::string(max_u_rc.cbegin() + k - 1, max_u_rc.cbegin() + min_idx_r + k - 1);
            }
            else if(max_unitig > max_u_rc)
                max_unitig = max_u_rc;

            // TODO: decide record-ID choice.
            output += FASTA_Record(0, max_unitig);

            j = e;
        }
    };

    std::cerr << "Peak-RAM before collation: " << process_peak_memory() / (1024.0 * 1024.0 * 1024.0) << "\n";
    parlay::parallel_for(0, max_unitig_bucket_count, collate_max_unitig_bucket, 1);
    std::cerr << "Peak-RAM after collation:  " << process_peak_memory() / (1024.0 * 1024.0 * 1024.0) << "\n";
}


template <uint16_t k>
std::size_t Unitig_Collator<k>::load_path_info(const std::size_t b, Path_Info<k>* const M, Buffer<unitig_path_info_t>& buf)
{
    buf.reserve(P_e[b].data().size());  // TODO: perform one fixed resize beforehand, as the `P_e` buckets will not grow anymore.
    const std::size_t b_sz = P_e[b].data().load(buf.data());
    assert(b_sz <= max_bucket_sz);

    for(std::size_t idx = 0; idx < b_sz; ++idx)
    {
        const auto& p_e = buf[idx];
        assert(p_e.obj() < max_bucket_sz);

        M[p_e.obj()] = p_e.path_info();
    }

    return b_sz;
}


template <uint16_t k>
std::size_t Unitig_Collator<k>::min_kmer(const std::string& s, Kmer<k>& min)
{
    Kmer<k> kmer(s), dummy;
    std::size_t min_idx = 0;
    min = kmer;

    std::size_t curr_idx = 1;
    while(curr_idx + k <= s.length())
    {
        kmer.roll_to_next_kmer(s[curr_idx + k - 1], dummy);
        if(min > kmer)
            min = kmer, min_idx = curr_idx;

        curr_idx++;
    }

    return min_idx;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Unitig_Collator)
