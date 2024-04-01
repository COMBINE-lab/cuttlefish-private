
#ifndef EXT_MEM_BUCKET_HPP
#define EXT_MEM_BUCKET_HPP



#include <cstddef>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <filesystem>
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
private:

    static constexpr std::size_t in_memory_bytes = 16lu * 1024; // 16KB.

    const std::string file_path;    // Path to the file storing the bucket.
    const std::size_t max_write_buf_bytes;  // Maximum size of the in-memory write-buffer in bytes.
    const std::size_t max_write_buf_elems;  // Maximum size of the in-memory write-buffer in elements.

    // TODO: replace vector with custom container.
    std::vector<T_> buf;    // In-memory buffer of the bucket-elements.
    std::size_t size_;  // Number of elements added to the bucket.

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

    // Returns the size of the bucket.
    std::size_t size() const { return size_; }

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

    // Closes the bucket. Elements should not be added anymore once this has
    // been invoked. This method is required only if the entirety of the bucket
    // needs to live in external-memory after the parent process finishes.
    void close();

    // Loads the bucket into the vector `v`.
    void load(std::vector<T_>& v) const;

    // Loads the bucket into `b` and returns its size.
    std::size_t load(T_* b) const;
};


template <typename T_>
inline Ext_Mem_Bucket<T_>::Ext_Mem_Bucket(const std::string& file_path, const std::size_t buf_sz):
      file_path(file_path)
    , max_write_buf_bytes(buf_sz)
    , max_write_buf_elems(buf_sz / sizeof(T_))
    , size_(0)
{
    buf.reserve(max_write_buf_elems);

    if(!file_path.empty())
        file.open(file_path, std::ios::out | std::ios::binary);

    if(!file)
    {
        std::cerr << "Error opening external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::add(const T_& elem)
{
    buf.push_back(elem);
    size_++;
    if(buf.size() >= max_write_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::add(const T_* const buf, const std::size_t sz)
{
    std::for_each(buf, buf + sz, [&](const auto elem){ add(elem); });

    // size_ += sz;
    // if(this->buf.size() + sz >= max_write_buf_elems)
    //     file.write(reinterpret_cast<const char*>(buf), sz * sizeof(T_));    // TODO: thrashing possible with an almost full buffer.
    // else
    //     std::for_each(buf, buf + sz, [&](const auto elem){ this->buf.push_back(elem); });
}


template <typename T_>
template <typename... Args>
inline void Ext_Mem_Bucket<T_>::emplace(Args&&... args)
{
    buf.emplace_back(std::forward<Args>(args)...);  // TODO: learn details.
    size_++;
    if(buf.size() >= max_write_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::flush()
{
    file.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
    if(!file)
    {
        std::cerr << "Error writing to external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    buf.clear();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::close()
{
    if(!buf.empty())
        flush();

    file.close();
    if(!file)
    {
        std::cerr << "Error closing external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::load(std::vector<T_>& v) const
{
    std::error_code ec;
    const auto file_sz = std::filesystem::file_size(file_path, ec);

    assert(file_sz % sizeof(T_) == 0);
    assert(file_sz / sizeof(T_) + buf.size() == size_);

    v.resize(size_);

    std::ifstream input(file_path);
    input.read(reinterpret_cast<char*>(v.data()), file_sz);
    input.close();

    if(ec || !input)
    {
        std::cerr << "Error reading of external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    std::memcpy(reinterpret_cast<char*>(v.data()) + file_sz, reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
}


template <typename T_>
inline std::size_t Ext_Mem_Bucket<T_>::load(T_* b) const
{
    std::error_code ec;
    const auto file_sz = std::filesystem::file_size(file_path);

    assert(file_sz % sizeof(T_) == 0);
    assert(file_sz / sizeof(T_) + buf.size() == size_);

    std::ifstream input(file_path);
    input.read(reinterpret_cast<char*>(b), file_sz);
    input.close();

    if(ec || !input)
    {
        std::cerr << "Error reading of external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    if(!buf.empty())
        std::memcpy(reinterpret_cast<char*>(b) + file_sz, reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));

    return size_;
}

}



#endif
