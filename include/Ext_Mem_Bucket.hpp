
#ifndef EXT_MEM_BUCKET_HPP
#define EXT_MEM_BUCKET_HPP



#include <cstddef>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <cstdlib>


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
    const std::size_t max_buf_bytes;    // Maximum size of the in-memory buffer in bytes.
    const std::size_t max_buf_elems;    // Maximum size of the in-memory buffer in elements.

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

    // Adds the element `elem` to the bucket.
    void add(const T_& elem);

    // Emplaces an element, with its constructor-arguments being `args`, into
    // the bucket.
    template <typename... Args> void emplace(Args&&... args);

    // Closes the bucket.
    void close();

    // Loads the bucket into the vector `v`.
    void load(std::vector<T_>& v) const;
};


template <typename T_>
inline Ext_Mem_Bucket<T_>::Ext_Mem_Bucket(const std::string& file_path, const std::size_t buf_sz):
      file_path(file_path)
    , max_buf_bytes(buf_sz)
    , max_buf_elems(buf_sz / sizeof(T_))
    , size_(0)
{
    buf.reserve(max_buf_elems);

    if(!file_path.empty())
        file.open(file_path);
}


template <typename T_>
void Ext_Mem_Bucket<T_>::add(const T_& elem)
{
    buf.push_back(elem);
    size_++;
    if(buf.size() >= max_buf_elems)
        flush();
}


template <typename T_>
template <typename... Args>
void Ext_Mem_Bucket<T_>::emplace(Args&&... args)
{
    buf.emplace_back(std::forward<Args>(args)...);  // TODO: learn details.
    size_++;
    if(buf.size() >= max_buf_elems)
        flush();
}


template <typename T_>
inline void Ext_Mem_Bucket<T_>::flush()
{
    file.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(T_));
    buf.clear();
}


template <typename T_>
void Ext_Mem_Bucket<T_>::close()
{
    if(!buf.empty())
        flush();

    file.close();
}


template <typename T_>
void Ext_Mem_Bucket<T_>::load(std::vector<T_>& v) const
{
    std::error_code ec;
    const auto file_sz = std::filesystem::file_size(file_path, ec);

    v.resize(file_sz / sizeof(T_));

    std::ifstream input(file_path);
    input.read(reinterpret_cast<char*>(v.data()), file_sz);
    input.close();

    if(ec || !input)
    {
        std::cerr << "Error reading of external-memory bucket at " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}

}



#endif
