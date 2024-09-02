
#ifndef UNITIG_COLLATOR_HPP
#define UNITIG_COLLATOR_HPP



#include "dBG_Contractor.hpp"
#include "Path_Info.hpp"
#include "Unitig_Coord_Bucket.hpp"
#include "DNA_Utility.hpp"
#include "Directed_Vertex.hpp"
#include "Discontinuity_Graph.hpp"
#include "globals.hpp"
#include "utility.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <vector>
#include <string>
#include <cassert>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Collates locally-maximal unitigs from different unitig-buckets as per their
// path-information in a discontinuity graph of `k`-mers.
template <uint16_t k, bool Colored_>
class Unitig_Collator
{
private:

    Discontinuity_Graph<k, Colored_>& G;    // The discontinuity-graph.

    typedef typename dBG_Contractor<k>::unitig_path_info_t unitig_path_info_t;

    typedef typename dBG_Contractor<k>::P_e_t P_e_t;
    P_e_t& P_e; // `P_e[b]` contains path-info for edges in bucket `b`.

    const std::string lmtig_buckets_path;   // Path-prefix to the lm-tig buckets.
    const std::string unitig_coord_buckets_path;    // Path-prefix to the unitig-coordinate buckets produced in map-reduce.

    std::size_t max_bucket_sz;  // Maximum size of the edge-buckets.

    const std::size_t max_unitig_bucket_count;  // Number of buckets storing literal globally-maximal unitigs.
    // std::vector<Padded<Unitig_Coord_Bucket<k>>> max_unitig_bucket; // Key-value collation buckets for lm-unitigs.
    std::vector<Padded<Unitig_Coord_Bucket_Concurrent<k, Colored_>>> max_unitig_bucket; // Key-value collation buckets for lm-unitigs.

    typedef typename dBG_Contractor<k>::op_buf_list_t op_buf_list_t;
    op_buf_list_t& op_buf;  // Worker-specific output buffers.

    class Maximal_Unitig_Label;


    // Maps each locally-maximal unitig to its maximal unitig's corresponding
    // bucket.
    void map();

    // Reduces each maximal unitig bucket to its contained maximal unitigs.
    void reduce();

    // Loads the path-info of edges from bucket `b` into the table `M`, and
    // returns the size of the bucket. Uses the buffer `buf` to transfer the
    // information from the bucket to the table.
    std::size_t load_path_info(std::size_t b, Path_Info<k>* M, Buffer<unitig_path_info_t>& buf);

    // Emits the trivially maximal unitigs to the output stream. Only
    // applicable in the colored case.
    void emit_trivial_mtigs();


public:

    // Constructs a unitig-collator for unitigs with their associated path-info
    // at `P_e`, i.e. `P_e[b]` contains path-information of the unitigs'
    // corresponding edges at bucket `b`. `logistics` is the data logistics
    // manager for the algorithm execution. Worker-specific maximal unitigs are
    // written to the buffers in `op_buf`. `gmtig_bucket_count` many buckets
    // are used to partition the lm-tigs to their maximal unitigs. `G` is the
    // associated discontinuity graph.
    Unitig_Collator(Discontinuity_Graph<k, Colored_>& G, P_e_t& P_e, const Data_Logistics& logistics, op_buf_list_t& op_buf, std::size_t gmtig_bucket_count);

    // Collates the locally-maximal unitigs into global ones.
    void collate();
};


// Label of a maximal unitig.
template  <uint16_t k, bool Colored_>
class Unitig_Collator<k, Colored_>::Maximal_Unitig_Label
{
private:

    Buffer<char> label_;    // Label-sequence.
    std::size_t sz; // Size of the label.

    Buffer<char> cycle_buf; // Working-space to process cyclic maximal unitigs.


    // Appends the `len_s`-length sequence `s` to the label. `RC_` specifies
    // whether `s` needs to be reverse-complemented or not.
    template <bool RC_> void append(const char* s, std::size_t len_s);


public:

    // Returns the label sequence.
    auto data() const { return label_.data(); }

    // Returns the size of the label.
    auto size() const { return  sz; }

    // Clears the label.
    auto clear() { sz = 0; }

    // Returns `true` iff the label is empty.
    auto empty() const { return sz == 0; }

    // Initializes the label with the sequence `unitig`. `rc` specifies whether
    // `unitig` needs to be put in its reverse-complemented form.
    void init(const std::string_view& unitig, bool rc);

    // Appends the `k`-overlapping sequence `unitig` to the label. `rc`
    // specifies whether `unitig` needs to be added in its reverse-
    // complemented form.
    void append(const std::string_view& unitig, bool rc);

    // Removes the last character of the label.
    void pop_back() { assert(sz > 0); sz--; }

    // Transforms the label to its canonical form.
    void canonicalize();

    // Transforms the label to its canonical form given that the maximal unitig
    // is a cycle.
    void canonicalize_cycle();
};


template <uint16_t k, bool Colored_>
template <bool RC_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig_Label::append(const char* const s, const std::size_t len_s)
{
    if constexpr(!RC_)
        std::memcpy(label_.data() + sz, s, len_s);
    else
        for(std::size_t i = 0; i < len_s; ++i)
            label_[sz + i] = DNA_Utility::complement(s[len_s - 1 - i]);

    sz += len_s;
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig_Label::init(const std::string_view& unitig, const bool rc)
{
    clear();

    label_.reserve_uninit(unitig.length());
    rc ?    append<true>(unitig.data(), unitig.length()) :
            append<false>(unitig.data(), unitig.length());

    sz = unitig.length();
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig_Label::append(const std::string_view& unitig, const bool rc)
{
    label_.reserve(sz + unitig.length() - k);
    if(!rc)
        append<false>(unitig.data() + k, unitig.length() - k);
    else
        append<true>(unitig.data(), unitig.length() - k);
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig_Label::canonicalize()
{
    for(std::size_t i = 0; i < k; ++i)
    {
        const auto b_fw = label_[i];
        const auto b_bw = DNA_Utility::complement(label_[sz - 1 - i]);

        if(b_fw < b_bw) // Already in canonical form.
            return;

        if(b_fw > b_bw) // Reverse-complement is the canonical form.
        {
            for(std::size_t j = 0; j < sz / 2; ++j)
            {
                const auto b_l = label_[j];
                const auto b_r = label_[sz - 1 - j];
                label_[j] = DNA_Utility::complement(b_r),
                label_[sz - 1 - j] = DNA_Utility::complement(b_l);
            }

            if(sz & 1)
                label_[sz / 2] = DNA_Utility::complement(label_[sz / 2]);

            return;
        }
    }
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig_Label::canonicalize_cycle()
{
    Directed_Vertex<k> v(label_.data());
    Kmer<k> min_fw(v.kmer());   // Minimum k-mer in the forward strand.
    Kmer<k> min_bw(v.kmer_bar());   // Minimum k-mer in the backward strand.
    std::size_t min_idx_fw = 0;
    std::size_t min_idx_bw = k - 1;

    for(std::size_t i = 1; i + k <= sz; ++i)
    {
        v.roll_forward(DNA_Utility::map_base(label_[i + k - 1]));

        if(min_fw > v.kmer())
            min_fw = v.kmer(), min_idx_fw = i;

        if(min_bw > v.kmer_bar())
            min_bw = v.kmer_bar(), min_idx_bw = i + k - 1;
    }


    cycle_buf.reserve_uninit(sz);
    if(min_fw < min_bw)
    {
        const auto len_r = sz - min_idx_fw;
        const auto len_l = sz - len_r;
        std::memcpy(cycle_buf.data(), label_.data() + len_l, len_r);
        std::memcpy(cycle_buf.data() + len_r, label_.data() + k - 1, len_l);
    }
    else
    {
        const int64_t len_l = min_idx_bw + 1;
        const int64_t len_r = sz - len_l;
        std::size_t p = 0;

        for(int64_t i = len_l - 1; i >= 0; --i)
            cycle_buf[p++] = DNA_Utility::complement(label_[i]);

        for(int64_t i = len_r - 1; i >= 0; --i)
            cycle_buf[p++] = DNA_Utility::complement(label_[len_l + i - (k - 1)]);
    }

    std::memcpy(label_.data(), cycle_buf.data(), sz);
}

}



#endif
