
#ifndef KMER_INDEX_UTILITY_HPP
#define KMER_INDEX_UTILITY_HPP



#include <cstdint>
#include <cstddef>
#include <vector>
#include <type_traits>
#include <fstream>
#include <iostream>


// =============================================================================
// A class containing various utility fields and methods for the k-mer indexing
// scheme of de Bruijn graphs based on minimizers.
class Kmer_Index_Utility
{
protected:

    constexpr static std::size_t buf_sz_th = 5 * 1024 * 1024;   // Threshold for the total size (in bytes) of the buffers per path-sequence producer: 5 MB.
    constexpr static double gamma = 2.0;    // The gamma parameter for the minimizer-MPHF.
    constexpr static std::size_t idx_lock_count = 65536;    // Number of locks in the sparse-locks used in various steps.

    // TODO: move to build-params.
    constexpr static std::size_t overflow_threshold = (1 << 5); // Threshold size for minimizer instances to be put in the overflow index.

    static constexpr char PATH_FILE_EXT[] = ".paths";
    static constexpr char PATH_END_FILE_EXT[] = ".ends";
    static constexpr char MPHF_FILE_EXT[] = ".min.mphf";
    static constexpr char COUNT_FILE_EXT[] = ".min.count";
    static constexpr char OFFSET_FILE_EXT[] = ".min.offset";
    static constexpr char CONFIG_FILE_EXT[] = ".idx.conf";  // TODO: maybe also output configuration to the CF json file?
    static constexpr char MIN_INST_FILE_EXT[] = ".mins";
    static constexpr char OVERFLOW_KMER[] = ".overflow.kmers";
    static constexpr char OVERFLOW_MIN_INST_IDX[] = ".overflow.offset";
    static constexpr char OVERFLOW_MPHF_FILE_EXT[] = ".overflow.mphf";
    static constexpr char OVERFLOW_KMER_MAP_EXT[] = ".overflow.map";


    // Returns the configuration file path for the k-mer index present at path
    // prefix `idx_path`.
    static const std::string config_file_path(const std::string& idx_path);

    // Reads the k-mer length of some k-mer index from its configuration file
    // at path `config_path` and returns it.
    static uint16_t kmer_len(const std::string& config_path);

    // Reads the minimizer length of some k-mer index from its configuration
    // file at path `config_path` and returns it.
    static uint16_t minimizer_len(const std::string& config_path);

    // Dumps the data from `container` to the stream `output`, clearing
    // `container`.
    template <typename T_container_>
    static void dump(T_container_& container, std::ofstream& output);

    // Binary searches for the maximum rightmost value in the container
    // `container` within the index range `[left, right]` (both ends inclusive)
    // that is at most as `val`. If such a value exists, returns its index.
    // Otherwise, returns `left - 1`.
    template <typename T_container_>
    static int64_t lower_bound(const T_container_& container, int64_t left, int64_t right, typename T_container_::value_type val);

    // Binary searches for the minimum leftmost value in the container
    // `container` within the index range `[left, right]` (both ends inclusive)
    // that is larger than `val`. If such a value exists, returns its index.
    // Otherwise, returns `right + 1`.
    template <typename T_container_>
    static int64_t upper_bound(const T_container_& container, int64_t left, int64_t right, typename T_container_::value_type val);
};


inline const std::string Kmer_Index_Utility::config_file_path(const std::string& idx_path)
{
    return idx_path + CONFIG_FILE_EXT;
}


inline uint16_t Kmer_Index_Utility::kmer_len(const std::string& config_path)
{
    std::ifstream config_file(config_path.c_str(), std::ios::in | std::ios::binary);

    uint16_t k;
    config_file.read(reinterpret_cast<char*>(&k), sizeof(k));

    if(!config_file)
    {
        std::cerr << "Error reading from the configuration file at " << config_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    return k;
}


inline uint16_t Kmer_Index_Utility::minimizer_len(const std::string& config_path)
{
    std::ifstream config_file(config_path.c_str(), std::ios::in | std::ios::binary);

    config_file.seekg(sizeof(uint16_t), std::ios::beg); // Skip the `k`-value.

    uint16_t l;
    config_file.read(reinterpret_cast<char*>(&l), sizeof(l));

    if(!config_file)
    {
        std::cerr << "Error reading from the configuration file at " << config_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    return l;
}


template <typename T_container_>
inline void Kmer_Index_Utility::dump(T_container_& container, std::ofstream& output)
{
    if(!output.write(reinterpret_cast<const char*>(container.data()), container.size() * sizeof(typename T_container_::value_type)))
    {
        std::cerr << "Error writing to file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    container.clear();
}


template <typename T_container_>
inline int64_t Kmer_Index_Utility::lower_bound(const T_container_& container, int64_t left, int64_t right, const typename T_container_::value_type val)
{
    int64_t mid, result = left - 1;

    while(left <= right)
    {
        mid = (left + right) >> 1;
        if(container[mid] > val)
            right = mid - 1;
        else
        {
            result = mid;
            left = mid + 1;
        }
    }

    return result;
}


template <typename T_container_>
inline int64_t Kmer_Index_Utility::upper_bound(const T_container_& container, int64_t left, int64_t right, const typename T_container_::value_type val)
{
    int64_t mid, result = right + 1;

    while(left <= right)
    {
        mid = (left + right) >> 1;
        if(container[mid] <= val)
            left = mid + 1;
        else
        {
            result = mid;
            right = mid - 1;
        }
    }

    return result;
}



#endif
