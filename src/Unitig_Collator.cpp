
#include "Unitig_Collator.hpp"
#include "Unitig_File.hpp"
#include "Kmer.hpp"
#include "Spin_Lock.hpp"
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
#include <tuple>
#include <utility>
#include <algorithm>
#include <cassert>


// TODO: use more efficient data structures throughout the collator.


namespace cuttlefish
{

template <uint16_t k>
Unitig_Collator<k>::Unitig_Collator(const std::vector<Ext_Mem_Bucket_Concurrent<Obj_Path_Info_Pair<uni_idx_t, k>>>& P_e, const Data_Logistics& logistics, op_buf_list_t& op_buf):
      P_e(P_e)
    , lmtig_buckets_path(logistics.lmtig_buckets_path())
    , unitig_coord_buckets_path(logistics.unitig_coord_buckets_path())
    , M(nullptr)
    , p_e_buf(nullptr)
    , op_buf(op_buf)
{
    assert((max_unitig_bucket_count & (max_unitig_bucket_count - 1)) == 0);

    max_bucket_sz = 0;
    std::size_t sum_bucket_sz = 0;
    std::for_each(P_e.cbegin(), P_e.cend(), [&](const auto& bucket){ max_bucket_sz = std::max(max_bucket_sz, bucket.size()); sum_bucket_sz += bucket.size(); });

    std::cerr << "Sum edge-bucket size: " << sum_bucket_sz << "\n";
    std::cerr << "Maximum edge-bucket size: " << max_bucket_sz << "\n";
}


template <uint16_t k>
void Unitig_Collator<k>::par_collate()
{
    map();

    reduce();

/*
    std::string max_unitig, max_u_rc;
    std::size_t mu_tig = 0;
    auto& op = op_buf[parlay::worker_id()].data();
    Unitig_File_Reader mu_tig_reader(work_path + std::string("mutig"));
    while(mu_tig_reader.read_next_unitig(max_unitig))
    {
        max_u_rc = max_unitig;
        reverse_complement(max_u_rc);
        if(max_unitig > max_u_rc)
            max_unitig = max_u_rc;

        op += FASTA_Record(0, max_unitig);

        mu_tig++;
    }

*/
    // TODO: print meta-information over the maximal unitigs'.
}


template <uint16_t k>
void Unitig_Collator<k>::map()
{
    typedef Path_Info<k>* map_t;
    typedef unitig_path_info_t* buf_t;

    std::vector<Padded_Data<map_t>> M_vec(parlay::num_workers());   // Worker-local DATs of edge-index to edge path-info.
    std::vector<Padded_Data<buf_t>> buf_vec(parlay::num_workers()); // Worker-local buffers to read in edge path-info.

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            M_vec[w_id] = allocate<Path_Info<k>>(max_bucket_sz);
            buf_vec[w_id] = allocate<unitig_path_info_t>(max_bucket_sz);
        }, 1);


    std::vector<Spin_Lock> lock(max_unitig_bucket_count);   // Locks for maximal-unitig buckets.
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
        auto const M = M_vec[w_id].data();
        auto const buf = buf_vec[w_id].data();

        const auto b_sz = load_path_info(b, M, buf);
        edge_c += b_sz;

#ifndef NDEBUG
        uint64_t h = 0;
        for(std::size_t i = 0; i < b_sz; ++i)
            h ^= M[i].hash();

        h_p_e ^= h;
#endif


        assert(file_exists(lmtig_buckets_path + "_" + std::to_string(b)));
        Unitig_File_Reader unitig_reader(lmtig_buckets_path + "_" + std::to_string(b));
        std::string unitig; // Read-off unitig. TODO: use better-suited container.
        uni_idx_t idx = 0;  // The unitig's sequential ID in the bucket.
        std::size_t uni_len;    // The unitig's length in bases.
        for(; (uni_len = unitig_reader.read_next_unitig(unitig)) > 0; idx++)
        {
            assert(idx < b_sz);
            const auto p = M[idx].p();  // Path-ID of this unitig.
            const auto mapped_b_id = XXH3_64bits(&p, sizeof(p)) & (max_unitig_bucket_count - 1); // Hashed maximal unitig bucket.

            lock[mapped_b_id].lock();
            max_unitig_bucket[mapped_b_id].data().add(M[idx], unitig.data(), uni_len);
            lock[mapped_b_id].unlock();
        }

        assert(idx == b_sz);
    };

    parlay::parallel_for(1, P_e.size(), map_to_max_unitig_bucket, 1);

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            deallocate(M_vec[w_id].data());
            deallocate(buf_vec[w_id].data());
        }, 1);


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
    for(std::size_t i = 0; i < max_unitig_bucket_count; ++i)
        max_max_uni_b_sz = std::max(max_max_uni_b_sz, max_unitig_bucket[i].data().size()),
        max_max_uni_b_label_len = std::max(max_max_uni_b_label_len, max_unitig_bucket[i].data().label_len());

    std::cerr << "Maximum maximal unitig bucket size:  " << max_max_uni_b_sz << "\n";
    std::cerr << "Maximum maximal unitig label length: " << max_max_uni_b_label_len << "\n";


    typedef Unitig_Coord<k>* coord_buf_t;
    typedef char* label_buf_t;

    std::vector<Padded_Data<coord_buf_t>> U_vec(parlay::num_workers()); // Worker-local buffers for unitig coordinate information.
    std::vector<Padded_Data<label_buf_t>> L_vec(parlay::num_workers()); // Worker-local buffers for dump-strings in buckets.

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            U_vec[w_id] = allocate<Unitig_Coord<k>>(max_max_uni_b_sz);
            L_vec[w_id] = allocate<char>(max_max_uni_b_label_len);
        }, 1);

    // TODO: add per-worker progress tracker.

    const auto collate_max_unitig_bucket =
    [&](const std::size_t b)
    {
        const auto w_id = parlay::worker_id();
        auto const U = U_vec[w_id].data();  // Coordinate info of the unitigs.
        auto const L = L_vec[w_id].data();  // Dump-strings of the unitig labels.
        auto& output = op_buf[w_id].data(); // Output buffer for the maximal unitigs.

        const auto b_sz = max_unitig_bucket[b].data().load_coords(U);
        const auto len = max_unitig_bucket[b].data().load_labels(L);
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

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            deallocate(U_vec[w_id].data());
            deallocate(L_vec[w_id].data());
        }, 1);
}


/*
template <uint16_t k>
void Unitig_Collator<k>::collate()
{
    const std::size_t unitig_bucket_count = P_e.size() - 1;

    typedef max_unitig_id_t<k> path_id_t;
    typedef std::tuple<weight_t, std::string, side_t> lmtig_info_t;
    typedef std::pair<path_id_t, lmtig_info_t> kv_t;    // (p, <r, u ,o>)
    std::vector<kv_t> kv_store;

    M = allocate<Path_Info<k>>(max_bucket_sz);
    p_e_buf = allocate<unitig_path_info_t>(max_bucket_sz);

    std::ofstream output(output_path);
    std::string unitig;
    uint64_t e_c = 0;
#ifndef NDEBUG
    uint64_t h_p_e = 0;
#endif
    for(std::size_t b = 1; b <= unitig_bucket_count; ++b)
    {
        const auto b_sz = load_path_info(b);
        e_c += b_sz;

#ifndef NDEBUG
        for(std::size_t i = 0; i < b_sz; ++i)
            h_p_e ^= M[i].hash();
#endif

        Unitig_File_Reader unitig_reader(work_path + std::string("lmutig_") + std::to_string(b));
        uni_idx_t uni_idx = 0;
        std::size_t uni_len;
        while((uni_len = unitig_reader.read_next_unitig(unitig)))
        {
            assert(uni_idx < b_sz);
            const auto p_e = M[uni_idx];

            kv_store.emplace_back(p_e.p(), lmtig_info_t(p_e.r(), unitig, p_e.o()));
            uni_idx++;
        }

        assert(uni_idx == b_sz);
    }

    std::cerr << "Found " << e_c << " edges.\n";
#ifndef NDEBUG
    std::cerr << "Edges' path-information signature: " << h_p_e << "\n";
#endif

    std::cerr << "KV-store has " << kv_store.size() << " pairs.\n";
    std::sort(kv_store.begin(), kv_store.end());


    std::string max_unitig, max_u_rc;
    std::size_t i, j;
    std::size_t max_path_sz = 0;
    std::size_t min_sz = std::numeric_limits<std::size_t>::max();
    std::size_t max_sz = 0;
    std::size_t max_lmtig_sz = 0;
    for(i = 0; i < kv_store.size(); i = j)
    {
        max_unitig.clear();

        const std::size_t s = i;
        std::size_t e;
        for(e = s + 1; e < kv_store.size() && kv_store[e].first == kv_store[s].first; ++e);

        if(e - s == 2)
        {
            auto& [r0, u0, o0] = kv_store[s].second;
            if(o0 == side_t::front)
                reverse_complement(u0);

            auto& [r1, u1, o1] = kv_store[s + 1].second;
            if(o1 == side_t::front)
                reverse_complement(u1);

            reverse_complement(u1);

            max_unitig = u0 + std::string(u1.begin() + k, u1.end());

            j = e;
            max_lmtig_sz = std::max(max_lmtig_sz, std::max(u0.size(), u1.size()));
        }
        else
            for(j = i; j < kv_store.size() && kv_store[j].first == kv_store[i].first; ++j)
            {
                auto& [r, u, o] = kv_store[j].second;
                if(o == side_t::front)
                    reverse_complement(u);

                max_unitig += (max_unitig.empty() ? u : std::string(u.begin() + k, u.end()));

                max_lmtig_sz = std::max(max_lmtig_sz, u.size());
            }

        max_u_rc = max_unitig;
        reverse_complement(max_u_rc);
        if(max_unitig > max_u_rc)
            max_unitig = max_u_rc;

        output << max_unitig << "\n";

        max_path_sz = std::max(max_path_sz, j - i);
        min_sz = std::min(min_sz, max_unitig.size());
        max_sz = std::max(max_sz, max_unitig.size());
    }

    std::size_t mu_tig = 0;
    Unitig_File_Reader mu_tig_reader(work_path + std::string("mutig"));
    while(mu_tig_reader.read_next_unitig(max_unitig))
    {
        max_u_rc = max_unitig;
        reverse_complement(max_u_rc);
        if(max_unitig > max_u_rc)
            max_unitig = max_u_rc;

        mu_tig++;
        output << max_unitig << "\n";

        min_sz = std::min(min_sz, max_unitig.size());
        max_sz = std::max(max_sz, max_unitig.size());
    }

    output.close();

    deallocate(M);
    deallocate(p_e_buf);

    std::cerr << "Read " << mu_tig << " trivially maximal unitigs.\n";


    std::cerr << "Longest path length in discontinuity graph: " << max_path_sz << "\n";
    std::cerr << "Minimum maximal-unitig size: " << min_sz << "\n";
    std::cerr << "Maximum maximal-unitig size: " << max_sz << "\n";
    std::cerr << "Maximum locally-maximal unitig size: " << max_lmtig_sz << "\n";
}
*/


template <uint16_t k>
std::size_t Unitig_Collator<k>::load_path_info(const std::size_t b, Path_Info<k>* const M, unitig_path_info_t* const buf)
{
    const std::size_t b_sz = P_e[b].load(buf);
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
std::size_t Unitig_Collator<k>::load_path_info(const std::size_t b)
{
    return load_path_info(b, M, p_e_buf);
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
