
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
private:

    static constexpr std::size_t in_memory_bytes = 16lu * 1024; // 16 KB.
    static constexpr auto in_memory_off = in_memory_bytes / sizeof(std::size_t);

    const std::string file_path;    // Path to the file for the unitig content.

    std::vector<char> buf;  // In-memory buffer for the unitig content.

    std::size_t total_sz;   // Total size of the added unitig content.
    std::vector<std::size_t> offset;    // Offsets of the unitigs into the file.
    std::size_t unitig_c;   // Number of unitigs.
    std::ofstream output; // The unitig file.
    std::ofstream output_off;   // The offsets file.


    // Returns path to the file containing the offsets of the unitigs in their
    // concatenated file.
    const std::string offset_file_path() const { return file_path + std::string(".off"); }

    // Flushes the in-memory unitig content to external memory.
    void flush_unitigs();

    // Flushes the in-memory unitig offsets to external memory.
    void flush_offsets();


public:

    Unitig_File_Writer(const std::string& file_path);

    auto unitig_count() const { return unitig_c; }

    template <typename T_it_> void add(T_it_ beg, T_it_ end);

    void close();
};


template <typename T_it_>
inline void Unitig_File_Writer::add(const T_it_ beg, const T_it_ end)
{
    offset.push_back(total_sz);
    buf.insert(buf.end(), beg, end);
    total_sz += (end - beg);
    unitig_c++;

    if(buf.size() >= in_memory_bytes)
        flush_unitigs();

    if(offset.size() >= in_memory_off)
        flush_offsets();
}


inline void Unitig_File_Writer::flush_unitigs()
{
    output.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(char));
    buf.clear();
}


inline void Unitig_File_Writer::flush_offsets()
{
    output_off.write(reinterpret_cast<const char*>(offset.data()), offset.size() * sizeof(decltype(offset)::value_type));
    offset.clear();
}

}



#endif
