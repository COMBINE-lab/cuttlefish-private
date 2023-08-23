
#ifndef CONCURRENT_HASH_TABLE_HPP
#define CONCURRENT_HASH_TABLE_HPP



#include "xxHash/xxhash.h"

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <cmath>


namespace cuttlefish
{

// =============================================================================
template <typename T_key_, typename T_val_, typename T_hasher_>
class Concurrent_Hash_Table
{
private:

    struct Key_Val_Pair
    {
        T_key_ key;
        T_val_ val;
    };

    static constexpr double lf_default = 0.75;   // Default maximum load-factor supported.

    T_key_ empty_key_;  // The empty key; currently it's set to all 1-bits.

    const T_hasher_ hash;   // The hasher object.

    const std::size_t capacity_;    // True capacity of the table; adjusted to be a power of 2.
    const std::size_t idx_wrapper_mask; // Bitmask to wrap indexing into the table.

    Key_Val_Pair* const T;  // The flat table of key-value collection.
    // std::size_t size_;  // Number of elements in the table. Probably not usable cheaply with concurrency.

    // Maps the hash value `h` to an index into the table.
    std::size_t hash_to_idx(std::size_t h) const { return h & idx_wrapper_mask; }

    // Returns the next (wrapped) index for `i`.
    std::size_t next_index(std::size_t i) const { return hash_to_idx(i + 1); }

    // Atomically compare-and-swaps the memory location `ptr`'s value to
    // `new_key` from `old_key`. Returns `true` iff this succeeds.
    static bool CAS(T_key_* ptr, T_key_ old_key, T_key_ new_key);

    // Returns the memory-equivalent value of `val` in type `T_to`.
    template <typename T_to_, typename T_from_>
    static T_to_ pun_type(const T_from_ val) { return *reinterpret_cast<const T_to_*>(&val); }


public:

    // Constructs a concurrent hash table to support upto `max_n` elements, with
    // a maximum load-factor of `load_factor`. The object `hasher` is used to
    // hash the keys in the table.
    Concurrent_Hash_Table(std::size_t max_n, double load_factor = lf_default, T_hasher_ hasher = T_hasher_());

    ~Concurrent_Hash_Table() { std::free(T); }

    // T_key_ empty_key() const { return empty_key_; }

    // Clears the hash table.
    // TODO: use `parlay::parallel_for` to efficiently fill the table with custom empty-key.
    void clear() { std::memset(reinterpret_cast<char*>(T), -1, capacity_ * sizeof(Key_Val_Pair)); }

    // Inserts the key `key` with value `val` into the table. Returns `false` if
    // the key already exists in the table. Otherwise returns `true` iff the
    // insertion succeeds, i.e. free space was found for the insertion. `mt_`
    // denotes whether multiple threads may access the hash table or not.
    template <bool mt_ = true> bool insert(T_key_ key, T_val_ val);

    // Inserts the key `key` with value `val` into the table. Returns `false` if
    // the key already exists in the table; and in that case, the address of the
    // existing value for the key is stored in `val_add`. Otherwise returns
    // `true` iff the insertion succeeds, i.e. free space was found for the
    // insertion. `mt_` denotes whether multiple threads may access the hash
    // table or not.
    template <bool mt_ = true> bool insert(T_key_ key, T_val_ val, T_val_*& val_add);

    // Inserts the key `key` with value `val` into the table. Returns `false` if
    // the key already exists in the table; and in that case, the existing value
    // associated to `key` is overwritten with `val`. Otherwise returns `true`
    // iff the insertion succeeds, i.e. free space was found for the insertion.
    // `mt_` denotes whether multiple threads may access the hash table or not.
    template <bool mt_ = true> bool insert_overwrite(T_key_ key, T_val_ val);

    // Searches for `key` in the table and returns the address of the value
    // associated to it iff it is found. Returns `nullptr` otherwise.
    const T_val_* find(const T_key_ key) const;

    // Searches for `key` in the table and returns `true` iff it is found. If
    // found, the associated value is stored in `val`. `val` remains unchanged
    // otherwise.
    bool find(const T_key_ key, T_val_& val) const;

    // Returns a 64-bit signature of the key-set of the hash table.
    uint64_t signature() const;

    // Returns a 64-bit signature of the values in the hash table.
    uint64_t signature_vals() const;
};


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Concurrent_Hash_Table(const std::size_t max_n, const double load_factor, const T_hasher_ hasher):
      empty_key_()
    , hash(hasher)
    , capacity_(static_cast<std::size_t>(1) << static_cast<std::size_t>(std::ceil(std::log2(max_n / load_factor))))
    , idx_wrapper_mask(capacity_ - 1)
    , T(static_cast<Key_Val_Pair*>(std::malloc(capacity_ * sizeof(Key_Val_Pair))))
{
    constexpr auto key_bytes = sizeof(T_key_);
    static_assert(  key_bytes == 1 || key_bytes == 2 || key_bytes == 4 || key_bytes == 8 || key_bytes == 16,
                    "Unsupported-sized keys for CAS in hash-table.");

    std::memset(reinterpret_cast<char*>(&empty_key_), -1, sizeof(empty_key_));

    clear();
}


template <typename T_key_, typename T_val_, typename T_hasher_>
template <bool mt_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert(const T_key_ key, const T_val_ val)
{
    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    {
        if(T[i].key == empty_key_)
        {
            if constexpr(!mt_)
            {
                T[i].key = key, T[i].val = val;
                return true;
            }
            else if(CAS(&T[i].key, empty_key_, key))
            {
                T[i].val = val;
                return true;
            }
        }

        if(T[i].key == key)
            return false;
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
template <bool mt_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert(const T_key_ key, const T_val_ val, T_val_*& val_add)
{
    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    {
        if(T[i].key == empty_key_)
        {
            if constexpr(!mt_)
            {
                T[i].key = key, T[i].val = val;
                return true;
            }
            else if(CAS(&T[i].key, empty_key_, key))
            {
                T[i].val = val;
                return true;
            }
        }

        if(T[i].key == key)
        {
            val_add = &T[i].val;
            return false;
        }
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
template <bool mt_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert_overwrite(const T_key_ key, const T_val_ val)
{
    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    {
        if(T[i].key == empty_key_)
        {
            if constexpr(!mt_)
            {
                T[i].key = key, T[i].val = val;
                return true;
            }
            else if(CAS(&T[i].key, empty_key_, key))
            {
                T[i].val = val;
                return true;
            }
        }

        if(T[i].key == key)
        {
            T[i].val = val;
            return false;
        }
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline const T_val_* Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find(const T_key_ key) const
{
    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
        if(T[i].key == key)
            return &T[i].val;
        else if(T[i].key == empty_key_)
            break;

    return nullptr;
}



template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find(const T_key_ key, T_val_& val) const
{
    const auto val_add = find(key);
    if(val_add != nullptr)
    {
        val = *val_add;
        return true;
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::CAS(T_key_* const ptr, const T_key_ old_key, const T_key_ new_key)
{
    constexpr auto key_bytes = sizeof(T_key_);
    static_assert(  key_bytes == 1 || key_bytes == 2 || key_bytes == 4 || key_bytes == 8 || key_bytes == 16,
                    "Unsupported-sized keys for CAS in hash-table.");

    if constexpr(sizeof(T_key_) == 1)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(ptr), pun_type<uint8_t>(old_key), pun_type<uint8_t>(new_key));

    if constexpr(sizeof(T_key_) == 2)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint16_t*>(ptr), pun_type<uint16_t>(old_key), pun_type<uint16_t>(new_key));

    if constexpr(sizeof(T_key_) == 4)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(ptr), pun_type<uint32_t>(old_key), pun_type<uint32_t>(new_key));

    if constexpr(sizeof(T_key_) == 8)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(ptr), pun_type<uint64_t>(old_key), pun_type<uint64_t>(new_key));

    if constexpr(sizeof(T_key_) == 16)
        return __sync_bool_compare_and_swap(reinterpret_cast<__uint128_t*>(ptr), pun_type<__uint128_t>(old_key), pun_type<__uint128_t>(new_key));

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline uint64_t Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::signature() const
{
    uint64_t sign = 0;
    for(std::size_t i = 0; i < capacity_; ++i)
        if(T[i].key != empty_key_)
            sign = sign ^ (XXH3_64bits(&T[i].key, sizeof(T[i].key)));

    return sign;
}


// Returns a 64-bit signature of the values in the hash table.
template <typename T_key_, typename T_val_, typename T_hasher_>
inline uint64_t Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>:: signature_vals() const
{
    uint64_t sign = 0;
    for(std::size_t i = 0; i < capacity_; ++i)
        if(T[i].key != empty_key_)
            sign = sign ^ XXH3_64bits(&T[i].val, sizeof(T[i].val));

    return sign;
}

}



#endif
