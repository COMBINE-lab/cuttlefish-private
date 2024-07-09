
#ifndef KMER_HASHTABLE_HPP
#define KMER_HASHTABLE_HPP



#include "Kmer.hpp"
#include "State_Config.hpp"
#include "utility.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <limits>
#include <cstring>
#include <vector>
#include <tuple>
#include <cmath>
#include <cassert>


namespace cuttlefish
{


// A fixed-sized hashtable for (k-mer, state) key-value pairs.
template <uint16_t k>
class Kmer_Hashtable
{
    class Iterator;
    friend class Iterator;

private:

    // Internal type of entries in the table.
    struct Key_Val_Entry
    {
        Kmer<k> key;
        State_Config val;
        uint8_t timestamp;  // Timestamp of the table version this entry belongs to.
    };

    // Type of each individual update operation into the table.
    struct Update_Entry
    {
        const Kmer<k> kmer;
        base_t front;
        base_t back;
        side_t disc_0;
        side_t disc_1;

        Update_Entry(const Kmer<k>& kmer, base_t front, base_t back, side_t disc_0, side_t disc_1);
    };


    static constexpr double lf_default = 0.75;  // Default maximum load-factor.

    const std::size_t capacity_;    // True capacity of the table; adjusted to be a power of 2.
    const std::size_t wrapper_mask; // Bitmask to wrap indexing into the table.

    Key_Val_Entry* const T; // The flat table of key-value collection.
    std::size_t sz; // Size of the hashtable.

    uint8_t cur_stamp;  // Current timestamp, or the version of the table.
    // TODO: it need not be just 1 byte, and can be more based on alignment spacing requirement of `Key_Val_Entry`.

    std::vector<Update_Entry> U;    // Buffer of update-entries.
    static constexpr std::size_t buf_cap = 1 * 1024 * 1024 / (sizeof(Update_Entry));    // 1MB worth of update entries.


    // Maps the hash value `h` to an index into the table.
    auto hash_to_idx(const std::size_t h) const { return h & wrapper_mask; }

    // Returns the next (wrapped) index to `i`.
    auto next_idx(const std::size_t idx) const { return hash_to_idx(idx + 1); }


public:

    // Constructs a hash table to support upto `max_n` k-mers, with a maximum
    // load-factor of `lf`.
    Kmer_Hashtable(std::size_t max_n, double lf = lf_default);

    Kmer_Hashtable(const Kmer_Hashtable&) = delete;
    Kmer_Hashtable(Kmer_Hashtable&&);
    Kmer_Hashtable& operator=(const Kmer_Hashtable&) = delete;
    Kmer_Hashtable& operator=(Kmer_Hashtable&&) = delete;

    ~Kmer_Hashtable() { deallocate(T); }

    // Returns the size of the hashtable.
    auto size() const { return sz; }

    // Returns the true capacity of the hashtable.
    auto capacity() const { return capacity_; }

    // Clears the hash table.
    void clear();

    // Signals an update to the hashtable for the k-mer `kmer`: for `front`-
    // and `back`-encoded edges at its front and back respectively, and
    // discontinuous sides `disc_0` and `disc_1`.
    void update(const Kmer<k>& kmer, base_t front, base_t back, side_t disc_0, side_t disc_1);

    // Flushes all buffered pending updates in the hashtable.
    void flush_updates();
};


template <uint16_t k>
inline Kmer_Hashtable<k>::Kmer_Hashtable(const std::size_t max_n, const double lf):
      capacity_(ceil_pow_2(static_cast<std::size_t>(std::ceil(max_n / lf))))
    , wrapper_mask(capacity_ - 1)
    , T(aligned_allocate<Key_Val_Entry>(capacity_))
    , sz(0)
    , cur_stamp(std::numeric_limits<uint8_t>::max())
{
    assert(capacity_ > 1);
    static_assert(buf_cap > 0);

    clear();
}


template <uint16_t k>
inline Kmer_Hashtable<k>::Kmer_Hashtable(Kmer_Hashtable&& rhs):
      capacity_(rhs.capacity_)
    , wrapper_mask(rhs.wrapper_mask)
    , T(rhs.T)
    , sz(rhs.sz)
    , cur_stamp(rhs.cur_stamp)
{
    // Moved objects are not really moved in C++.
    const_cast<Key_Val_Entry*&>(rhs.T) = nullptr;   // NOLINT(cppcoreguidelines-pro-type-const-cast)
}


template <uint16_t k>
inline Kmer_Hashtable<k>::Update_Entry::Update_Entry(const Kmer<k>& kmer, const base_t front, const base_t back, const side_t disc_0, const side_t disc_1):
      kmer(kmer)
    , front(front)
    , back(back)
    , disc_0(disc_0)
    , disc_1(disc_1)
{}


template <uint16_t k>
inline void Kmer_Hashtable<k>::clear()
{
    sz = 0;
    cur_stamp++;
    if(cur_stamp == 0)
    {
        std::memset(reinterpret_cast<char*>(T), 0, capacity_ * sizeof(Key_Val_Entry));
        cur_stamp = 1;
    }
}


template <uint16_t k>
inline void Kmer_Hashtable<k>::update(const Kmer<k>& kmer, const base_t front, const base_t back, const side_t disc_0, const side_t disc_1)
{
    U.emplace_back(kmer, front, back, disc_0, disc_1);

    if(U.size() == buf_cap)
        flush_updates();
}


template <uint16_t k>
inline void Kmer_Hashtable<k>::flush_updates()
{
    constexpr std::size_t batch_size = 64;

    uint64_t H[2][batch_size];

    std::size_t cur_batch_sz;
    std::size_t next_batch_sz = std::min(batch_size, U.size());
    for(std::size_t i = 0; i < next_batch_sz; ++i)
    {
        H[0][i] = hash_to_idx(U[i].kmer.to_u64());
        // __builtin_prefetch(T + H[0][i], 1, 0);
    }


    for(std::size_t g = 0; g * batch_size < U.size(); g++)
    {
        const std::size_t g_base = g * batch_size;
        cur_batch_sz = std::min(batch_size, U.size() - g_base);
        next_batch_sz = std::min(batch_size, U.size() - (g_base + cur_batch_sz));
        const auto* const h_cur = &H[g & 1][0];
        auto* const h_next = &H[(g + 1) & 1][0];

        // Prefetch next batch.
        const auto next_g_base = g_base + cur_batch_sz;
        for(std::size_t i = 0; i < next_batch_sz; ++i)
        {
            h_next[i] = hash_to_idx(U[next_g_base + i].kmer.to_u64());
            __builtin_prefetch(T + h_next[i], 1, 0);
        }

        // Process current batch.
        for(std::size_t i = 0; i < cur_batch_sz; ++i)
        {
            const auto& u = U[g_base + i];

            for(std::size_t j = h_cur[i]; ; j = next_idx(j))
                if(T[j].timestamp != cur_stamp) // Empty slot as it contains an entry from a previous version of the table.
                {
                    assert(size() < capacity());

                    T[j].timestamp = cur_stamp;

                    T[j].key = u.kmer;
                    T[j].val = State_Config();
                    T[j].val.update(u.front, u.back, u.disc_0, u.disc_1);
                    sz++;
                    break;
                }
                else if(T[j].key == u.kmer)
                {
                    T[j].val.update(u.front, u.back, u.disc_0, u.disc_1);
                    break;
                }
        }
    }


    U.clear();
}

}



#endif