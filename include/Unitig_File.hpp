
#ifndef UNITIG_FILE_HPP
#define UNITIG_FILE_HPP



#include "Virtual_File.hpp"

#include <cstddef>
#include <limits>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// Unitig-file writer manager.
class Unitig_File_Writer
{
    // TODO: use `u16` after testing done with `u32`.
    // typedef uint16_t uni_len_t; // Type of the length of a unitig in a bucket.
    typedef uint32_t uni_len_t; // Type of the length of a unitig in a bucket.

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


// =============================================================================
// Unitig-file reader manager.
class Unitig_File_Reader
{
    typedef uint32_t uni_idx_t; // Type of the index of a unitig in a bucket.
    // TODO: use `u16` after testing done with `u32`.
    // typedef uint16_t uni_len_t; // Type of the length of a unitig in a bucket.
    typedef uint32_t uni_len_t; // Type of the length of a unitig in a bucket.

private:

    static constexpr std::streamsize in_memory_bytes = 16lu * 1024; // 16 KB.

    const std::string file_path;    // Path to the file with the unitig content.

    std::vector<char> buf;  // In-memory buffer for the unitig content.
    std::vector<uni_len_t> uni_len; // Sizes of the unitigs in the current buffer.

    std::ifstream input;    // The unitigs-file.
    Virtual_File<uni_len_t> len;    // The lengths-file.

    std::size_t buf_idx;    // Index into the unitig-buffer for the next unitig to read-in.
    uni_idx_t uni_idx_in_file;  // Index of the next unitig to read from file.
    uni_idx_t uni_idx_in_mem;   // Index of the next unitig to parse from buffer.

    const uni_idx_t unitig_count_;  // Number of unitigs in the file.
    uni_idx_t unitig_parsed_;   // Number of unitigs parsed.
    std::size_t total_sz;   // Total size of the read unitig content.


    // Returns path to the file containing the lengths of the unitigs.
    const std::string length_file_path() const { return file_path + std::string(".len"); }


public:

    // Constructs a unitig-reader for the file at path `file_path`.
    Unitig_File_Reader(const std::string& file_path);

    // Returns the number of unitigs in the file.
    auto unitig_count() const { return unitig_count_; }

    // Reads the next unitig into `unitig`; returns `true` iff there were
    // unitigs remaining to be read.
    template <typename T_> bool read_next_unitig(T_& unitig);
};


template <typename T_it_>
inline void Unitig_File_Writer::add(const T_it_ beg, const T_it_ end)
{
    len.push_back(end - beg);
    buf.insert(buf.end(), beg, end);
    total_sz += (end - beg);
    unitig_c++;
    assert((end - beg) <= std::numeric_limits<uni_len_t>::max());

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


template <typename T_>
inline bool Unitig_File_Reader::read_next_unitig(T_& unitig)
{
    if(buf_idx == buf.size())   // Buffer has been parsed completely; try a re-read.
    {
        assert(uni_idx_in_mem == uni_len.size());

        if(uni_idx_in_file == unitig_count_)    // All unitigs have been read-off.
            return false;

        std::streamsize bytes_to_read = 0;
        uni_len.clear();

        while(uni_idx_in_file < unitig_count_ && bytes_to_read < in_memory_bytes)
        {
            uni_len.push_back(len[uni_idx_in_file]);
            bytes_to_read += uni_len.back();

            uni_idx_in_file++;
        }

        assert(bytes_to_read > 0);
        buf.resize(bytes_to_read);
        input.read(buf.data(), bytes_to_read * sizeof(char));
        assert(input.gcount() == bytes_to_read);

        buf_idx = 0;
        uni_idx_in_mem = 0;
    }


    const auto len = uni_len[uni_idx_in_mem];
    unitig.resize(len);
    std::memcpy(unitig.data(), buf.data() + buf_idx, len * sizeof(decltype(buf)::value_type));

    buf_idx += len;
    uni_idx_in_mem++;
    unitig_parsed_++;
    total_sz += len;
    assert(buf_idx <= buf.size());

    return true;
}

}



#endif
