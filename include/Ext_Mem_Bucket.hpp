
#ifndef EXT_MEM_BUCKET_HPP
#define EXT_MEM_BUCKET_HPP



#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <cstddef>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// An external-memory-backed bucket for elements of type `T_`.
template <typename T_>
class Ext_Mem_Bucket
{
public:

    static constexpr std::size_t in_memory_bytes = 16lu * 1024; // 16KB.

    const std::string file_path;    // Path to the file storing the bucket.
    const std::size_t max_buf_bytes;    // Maximum size of the in-memory write-buffer in bytes.
    const std::size_t max_buf_elems;    // Maximum size of the in-memory write-buffer in elements.

    T_* buf;    // In-memory buffer of the bucket-elements.
    std::size_t size_;  // Number of elements added to the bucket.

    std::size_t in_mem_size;    // Number of elements in the in-memory buffer.

    std::ofstream file; // The bucket-file.


    // Flushes the in-memory buffer content to external memory.
    void flush();


public:

    // Constructs an external-memory bucket at path `file_path`. An optional in-
    // memory buffer size (in bytes) `buf_sz` for the bucket can be specified.
    Ext_Mem_Bucket(const std::string& file_path, const std::size_t buf_sz = in_memory_bytes);

    // Constructs a placeholder bucket.
    Ext_Mem_Bucket(): Ext_Mem_Bucket("", 0)
    {}

    Ext_Mem_Bucket(Ext_Mem_Bucket&& rhs);

    ~Ext_Mem_Bucket() { deallocate(buf); }

    Ext_Mem_Bucket(const Ext_Mem_Bucket&) = delete;
    Ext_Mem_Bucket& operator=(const Ext_Mem_Bucket&) = delete;
    Ext_Mem_Bucket& operator=(Ext_Mem_Bucket&&) = delete;

    // Returns the size of the bucket.
    auto size() const { return size_; }

    // Adds the element `elem` to the bucket.
    void add(const T_& elem);

    // Adds `sz` elements from `buf` into the bucket.
    void add(const T_* buf, std::size_t sz);

    // TODO
    // Adds `sz` elements from `buf` into the bucket. The order of the elements
    // per their addition to the bucket may not be preserved.
    void add_unordered(const T_* buf, std::size_t sz);

    // Emplaces an element, with its constructor-arguments being `args`, into
    // the bucket.
    template <typename... Args> void emplace(Args&&... args);

    // Serializes and closes the bucket. Elements should not be added anymore
    // once this has been invoked. This method is required only if the entirety
    // of the bucket needs to live in external-memory after the parent process
    // finishes.
    void serialize();

    // Loads the bucket into `b` and returns its size.
    std::size_t load(T_* b) const;

    // Clears the bucket.
    void clear();

    // Removes the bucket.
    void remove();

    // Returns the resident set size of the space-dominant components of this
    // bucket.
    std::size_t RSS() const;
};


template <typename T_>
inline Ext_Mem_Bucket<T_>::Ext_Mem_Bucket(const std::string& file_path, const std::size_t buf_sz):
      file_path(file_path)
    , max_buf_bytes(buf_sz)
    , max_buf_elems(buf_sz / sizeof(T_))
    , buf(allocate<T_>(max_buf_elems))
    , size_(0)
    , in_mem_size(0)
{
    assert(file_path.empty() || max_buf_elems > 0);

    if(!file_path.empty())
        file.open(file_path, std::ios::binary);

    if(!file)
    {
        std::cerr << "Error opening external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_>
inline Ext_Mem_Bucket<T_>::Ext_Mem_Bucket(Ext_Mem_Bucket&& rhs):
      file_path(std::move(rhs.file_path))
    , max_buf_bytes(std::move(rhs.max_buf_bytes))
    , max_buf_elems(std::move(rhs.max_buf_elems))
    , buf(std::move(rhs.buf))
    , size_(std::move(rhs.size_))
    , in_mem_size(std::move(rhs.in_mem_size))
    , file(std::move(rhs.file))
{
    // Moved objects are not really moved in C++.
    rhs.buf = nullptr;
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::add(const T_& elem)
{
    buf[in_mem_size++] = elem;
    size_++;

    assert(in_mem_size <= max_buf_elems);
    if(in_mem_size == max_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::add(const T_* const buf, const std::size_t sz)
{
    std::size_t rem_sz = sz;
    std::size_t added = 0;
    while(rem_sz > 0)
    {
        const auto to_add = std::min(rem_sz, max_buf_elems - in_mem_size);
        std::memcpy(reinterpret_cast<char*>(this->buf + in_mem_size), reinterpret_cast<const char*>(buf + added), to_add * sizeof(T_));
        in_mem_size += to_add, added += to_add, rem_sz -= to_add;

        assert(in_mem_size <= max_buf_elems);
        if(in_mem_size == max_buf_elems)
            flush();
    }

    size_ += sz;
}


template <typename T_>
template <typename... Args>
inline void Ext_Mem_Bucket<T_>::emplace(Args&&... args)
{
    new(buf + in_mem_size) T_(args...);
    in_mem_size++;
    size_++;

    assert(in_mem_size <= max_buf_elems);
    if(in_mem_size == max_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::flush()
{
    assert(in_mem_size <= max_buf_elems);

    file.write(reinterpret_cast<const char*>(buf), in_mem_size * sizeof(T_));
    if(!file)
    {
        std::cerr << "Error writing to external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    in_mem_size = 0;
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::serialize()
{
    if(in_mem_size != 0)
        flush();

    deallocate(buf);
    buf = nullptr;

    file.close();
    if(!file)
    {
        std::cerr << "Error closing external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_>
inline std::size_t Ext_Mem_Bucket<T_>::load(T_* b) const
{
    const auto file_sz = (size_ - in_mem_size) * sizeof(T_);
    assert(file_sz <= file_size(file_path));
    load_file(file_path, file_sz, reinterpret_cast<char*>(b));

    assert(in_mem_size < max_buf_elems);
    if(in_mem_size > 0)
        std::memcpy(reinterpret_cast<char*>(b) + file_sz, reinterpret_cast<const char*>(buf), in_mem_size * sizeof(T_));

    return size_;
}


template <typename  T_>
inline void Ext_Mem_Bucket<T_>::clear()
{
    size_ = 0;
    in_mem_size = 0;
    file.seekp(std::ios_base::beg);
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::remove()
{
    if(!file_path.empty())
    {
        file.close();
        if(!file || !remove_file(file_path))
        {
            std::cerr << "Error removing file at " << file_path << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    deallocate(buf);
    buf = nullptr;
}


template <typename T_>
inline std::size_t Ext_Mem_Bucket<T_>::RSS() const
{
    return max_buf_elems * sizeof(T_);
}



// A concurrent external-memory bucket for elements of type `T_`.
template <typename T_>
class Ext_Mem_Bucket_Concurrent
{
private:

    static constexpr std::size_t in_memory_bytes = 32 * 1024;   // 32KB.

    const std::string file_path;    // Path to the file storing the bucket.
    const std::size_t max_buf_bytes;    // Maximum size of the in-memory worker-local write-buffers in bytes.
    const std::size_t max_buf_elems;    // Maximum size of the in-memory worker-local write-buffers in elements.

    std::size_t flushed;    // Number of elements added to the bucket and flushed to external-memory.

    std::vector<Padded<std::vector<T_>>> buf_w_local;   // In-memory worker-local buffers of the bucket-elements.

    std::ofstream file; // The bucket-file.
    Spin_Lock lock_;    // Lock to the bucket-file.


    // Flushes the in-memory buffer content of the invoking worker to external-
    // memory.
    void flush();


public:

    // Constructs a concurrent external-memory bucket at path `file_path`. An
    // optional in-memory buffer size (in bytes) `buf_sz` for each worker can
    // be specified.
    Ext_Mem_Bucket_Concurrent(const std::string& file_path, const std::size_t buf_sz = in_memory_bytes);

    // Constructs a placeholder bucket.
    Ext_Mem_Bucket_Concurrent(): Ext_Mem_Bucket_Concurrent("", 0)
    {}

    Ext_Mem_Bucket_Concurrent(Ext_Mem_Bucket_Concurrent&& rhs);

    Ext_Mem_Bucket_Concurrent(const Ext_Mem_Bucket_Concurrent&) = delete;
    Ext_Mem_Bucket_Concurrent& operator=(const Ext_Mem_Bucket_Concurrent&) = delete;
    Ext_Mem_Bucket_Concurrent& operator=(Ext_Mem_Bucket_Concurrent&&) = delete;

    // Returns the size of the bucket. It is exact only when the bucket is not
    // being updated. Otherwise it is not necessarily exact and runs the risk
    // of data races.
    std::size_t size() const;

    // Adds the element `elem` to the bucket.
    void add(const T_& elem);

    // Emplaces an element, with its constructor-arguments being `args`, into
    // the bucket.
    template <typename... Args> void emplace(Args&&... args);

    // Serializes and closes the bucket. Elements should not be added anymore
    // once this has been invoked. This method is required only if the entirety
    // of the bucket needs to live in external-memory after the parent process
    // finishes.
    void close();

    // Loads the bucket into the vector `v`. It is safe only when the bucket is
    // not being updated, otherwise runs the risk of data races.
    // TODO: remove this method and support only raw containers, such that user modules are forced to avoid adversarial resizing.
    void load(std::vector<T_>& v) const;

    // Loads the bucket into `b` and returns its size. `b` must have enough
    // space allocated for all the bucket elements. It is safe only when the
    // bucket is not being updated, otherwise runs the risk of data races.
    std::size_t load(T_* b) const;

    // Removes the bucket.
    void remove();

    // Returns the resident set size of the space-dominant components of this
    // bucket.
    std::size_t RSS() const;
};


template <typename T_>
inline Ext_Mem_Bucket_Concurrent<T_>::Ext_Mem_Bucket_Concurrent(const std::string& file_path, const std::size_t buf_sz):
      file_path(file_path)
    , max_buf_bytes(buf_sz)
    , max_buf_elems(max_buf_bytes / sizeof(T_))
    , flushed(0)
    , buf_w_local(parlay::num_workers())
{
    assert(file_path.empty() || max_buf_elems > 0);

    if(!file_path.empty())
        file.open(file_path, std::ios::out | std::ios::binary);

    if(!file)
    {
        std::cerr << "Error opening external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::for_each(buf_w_local.begin(), buf_w_local.end(), [&](auto& v){ v.unwrap().reserve(max_buf_elems); });
}


template <typename T_>
inline Ext_Mem_Bucket_Concurrent<T_>::Ext_Mem_Bucket_Concurrent(Ext_Mem_Bucket_Concurrent&& rhs):
      file_path(std::move(rhs.file_path))
    , max_buf_bytes(std::move(rhs.max_buf_bytes))
    , max_buf_elems(std::move(rhs.max_buf_elems))
    , flushed(std::move(rhs.flushed))
    , buf_w_local(std::move(rhs.buf_w_local))
    , file(std::move(rhs.file))
{}


template <typename T_>
inline std::size_t Ext_Mem_Bucket_Concurrent<T_>::size() const
{
    std::size_t in_buf_sz = 0;
    std::for_each(buf_w_local.cbegin(), buf_w_local.cend(), [&](const auto& b){ in_buf_sz += b.unwrap().size(); });

    return flushed + in_buf_sz;
}


template <typename T_>
inline void Ext_Mem_Bucket_Concurrent<T_>::add(const T_& elem)
{
    auto& buf = buf_w_local[parlay::worker_id()].unwrap();
    buf.push_back(elem);

    assert(buf.size() <= max_buf_elems);
    if(buf.size() == max_buf_elems)
        flush();
}


template <typename T_>
template <typename... Args>
inline void Ext_Mem_Bucket_Concurrent<T_>::emplace(Args&&... args)
{
    auto& buf = buf_w_local[parlay::worker_id()].unwrap();
    buf.emplace_back(args...);

    assert(buf.size() <= max_buf_elems);
    if(buf.size() == max_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket_Concurrent<T_>::close()
{
    std::size_t in_mem = 0;
    for(std::size_t w = 0; w < buf_w_local.size(); ++w)
        in_mem += buf_w_local[w].unwrap().size();

    auto& buf = buf_w_local[0].unwrap();
    buf.reserve(in_mem);
    for(std::size_t w = 1; w < buf_w_local.size(); ++w)
    {
        auto& b = buf_w_local[w].unwrap();
        buf.insert(buf.end(), b.cbegin(), b.cend());

        b.clear();
        force_free(b);
    }

    file.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
    if(!file)
    {
        std::cerr << "Error writing to external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    flushed += buf.size();

    buf.clear();
    force_free(buf);
}


template <typename T_>
inline void Ext_Mem_Bucket_Concurrent<T_>::flush()
{
    auto& buf = buf_w_local[parlay::worker_id()].unwrap();
    assert(buf.size() <= max_buf_elems);

    if(buf.empty())
        return;

    lock_.lock();

    // TODO: use async-write.
    file.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
    if(!file)
    {
        std::cerr << "Error writing to external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    flushed += buf.size();

    lock_.unlock();

    buf.clear();
}


template <typename T_>
inline void Ext_Mem_Bucket_Concurrent<T_>::load(std::vector<T_>& v) const
{
    const auto sz = size();
    v.resize(sz);

    // Load from the bucket-file.

    const auto file_sz = load_file(file_path, reinterpret_cast<char*>(v.data()));
    assert(file_sz == flushed * sizeof(T_));
    (void)file_sz;

    // Load the elements pending in the worker-local buffers.

    auto curr_end = v.data() + flushed;
    for(const auto& b : buf_w_local)
    {
        const auto buf = b.unwrap();
        if(CF_LIKELY(!buf.empty())) // Conditional to avoid UB on `nullptr` being passed to `memcpy`.
            std::memcpy(reinterpret_cast<char*>(curr_end), reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
        curr_end += buf.size();
    }
}


template <typename T_>
inline std::size_t Ext_Mem_Bucket_Concurrent<T_>::load(T_* b) const
{
    std::size_t sz = flushed;

    // Load from the bucket-file.

    const auto file_sz = load_file(file_path, reinterpret_cast<char*>(b));
    assert(file_sz == flushed * sizeof(T_));
    (void)file_sz;

    // Load the elements pending in the worker-local buffers.

    auto curr_end = b + flushed;
    for(const auto& buf_w : buf_w_local)
    {
        const auto buf = buf_w.unwrap();
        if(CF_LIKELY(!buf.empty())) // Conditional to avoid UB on `nullptr` being passed to `memcpy`.
            std::memcpy(reinterpret_cast<char*>(curr_end), reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
        curr_end += buf.size();
        sz += buf.size();
    }

    return sz;
}


template <typename T_>
inline void Ext_Mem_Bucket_Concurrent<T_>::remove()
{
    if(!file_path.empty())
    {
        file.close();
        if(!file || !remove_file(file_path))
        {
            std::cerr << "Error removing file at " << file_path << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    }


    std::for_each(buf_w_local.begin(), buf_w_local.end(), [](auto& w_buf){ force_free(w_buf.unwrap()); });
}


template <typename T_>
inline std::size_t Ext_Mem_Bucket_Concurrent<T_>::RSS() const
{
    std::size_t buf_bytes = 0;
    std::for_each(buf_w_local.cbegin(), buf_w_local.cend(), [&](const auto& b){ buf_bytes += b.unwrap().capacity() * sizeof(T_); });

    return buf_bytes;
}

}



#endif
