
#ifndef UNITIG_FILE_HPP
#define UNITIG_FILE_HPP



#include <cstddef>
#include <vector>
#include <string>
#include <fstream>


namespace cuttlefish
{

// =============================================================================
// Unitig-file writer manager.
class Unitig_File_Writer
{
    typedef uint16_t uni_len_t; // Type of the length of a unitig in a bucket.

private:

    static constexpr std::size_t in_memory_bytes = 16lu * 1024; // 16 KB.
    static constexpr auto in_memory_len = in_memory_bytes / sizeof(uni_len_t);

    const std::string file_path;    // Path to the file for the unitig content.

    std::vector<char> buf;  // In-memory buffer for the unitig content.

    std::size_t total_sz;   // Total size of the added unitig content.
    std::vector<uni_len_t> len; // Lengths of the unitigs in the file.
    std::size_t unitig_c;   // Number of unitigs added.
    std::ofstream output; // The unitig file.
    std::ofstream output_len;   // The lengths file.


    // Returns path to the file containing the lengths of the unitigs.
    const std::string length_file_path() const { return file_path + std::string(".len"); }

    // Flushes the in-memory unitig content to external memory.
    void flush_unitigs();

    // Flushes the in-memory unitig lengths to external memory.
    void flush_lengths();


public:

    // Constructs a unitig-writer to the file at path `file_path`.
    Unitig_File_Writer(const std::string& file_path);

    // Returns the total size of the added unitig content.
    auto size() const { return total_sz; }

    // Returns the number of unitigs added.
    auto unitig_count() const { return unitig_c; }

    // Adds the unitig content in the range `[beg, end)` to the writer.
    template <typename T_it_> void add(T_it_ beg, T_it_ end);

    // Closes the stream.
    void close();
};


template <typename T_it_>
inline void Unitig_File_Writer::add(const T_it_ beg, const T_it_ end)
{
    len.push_back(end - beg);
    buf.insert(buf.end(), beg, end);
    total_sz += (end - beg);
    unitig_c++;

    if(buf.size() >= in_memory_bytes)
        flush_unitigs();

    if(len.size() >= in_memory_len)
        flush_lengths();
}


inline void Unitig_File_Writer::flush_unitigs()
{
    output.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(decltype(buf)::value_type));
    buf.clear();
}


inline void Unitig_File_Writer::flush_lengths()
{
    output_len.write(reinterpret_cast<const char*>(len.data()), len.size() * sizeof(decltype(len)::value_type));
    len.clear();
}

}



#endif
