
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

    static constexpr char PATH_FILE_EXT[] = ".paths";
    static constexpr char PATH_END_FILE_EXT[] = ".ends";
    static constexpr char MPHF_FILE_EXT[] = ".min.mphf";
    static constexpr char COUNT_FILE_EXT[] = ".min.count";
    static constexpr char OFFSET_FILE_EXT[] = ".min.offset";
    static constexpr char CONFIG_FILE_EXT[] = ".idx.conf";  // TODO: maybe also output configuration to the CF json file?
    static constexpr char MIN_INST_FILE_EXT[] = ".mins";


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
