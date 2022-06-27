
#include "Kmer_Index.hpp"
#include "Input_Defaults.hpp"
#include "utility.hpp"

#include <string>
#include <algorithm>
#include <thread>


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const uint16_t l, const uint16_t producer_count):
    // TODO: reserve space for `paths`, preferably from additional k-mer count field
    l(l),
    producer_count(producer_count),
    producer_path_buf(producer_count),
    producer_minimizer_buf(producer_count),
    producer_minimizer_file(producer_count),
    curr_token{0}
{
    for(uint16_t id = 0; id < producer_count; ++id)
        producer_minimizer_file[id].open(minimizer_file_path(id), std::ios::out | std::ios::binary);
}


template <uint16_t k>
const std::string Kmer_Index<k>::minimizer_file_path(const uint16_t producer_id)
{
    return std::to_string(producer_id) + cuttlefish::_default::MINIMIZER_FILE_EXT;
}


template <uint16_t k>
const typename Kmer_Index<k>::Producer_Token Kmer_Index<k>::get_token()
{
    lock.lock();
    const std::size_t token = curr_token++;
    lock.unlock();

    return Producer_Token(token);
}


template <uint16_t k>
void Kmer_Index<k>::finalize_production()
{
    for(std::size_t id = 0; id < producer_count; ++id)
    {
        if(!producer_path_buf[id].empty())
            flush(id);

        producer_minimizer_file[id].close();
    }

    std::ofstream path_file(cuttlefish::_default::PATH_FILE_EXT, std::ios::out | std::ios::binary); // TODO: add ext. to o/p file name (from `Build_Params`).
    paths.serialize(path_file);
    // path_file.write(reinterpret_cast<const char*>(paths.data()), paths.size());  // For testing.
    path_file.close();
}


template <uint16_t k>
void Kmer_Index<k>::consolidate_minimizers()
{
    std::vector<Minimizer_Offset_Pair*> min_buf(producer_count);    // Separate buffer for each minimizer file.
    std::vector<std::thread> worker;    // Worker threads.

    // TODO: bound memory usage within a provided argument.

    // Launch sorter threads.
    for(uint16_t id = 0; id < producer_count; ++id)
    {
        const auto min_buf_bytes = file_size(minimizer_file_path(id));   // What happens if it's zero?
        min_buf[id] = static_cast<Minimizer_Offset_Pair*>(std::malloc(min_buf_bytes));

        std::ifstream input(minimizer_file_path(id).c_str(), std::ios::in | std::ios::binary);
        if(!input.read(reinterpret_cast<char*>(min_buf[id]), min_buf_bytes))
        {
            std::cerr << "Error reading the minimizer files. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        const std::size_t min_buf_size = min_buf_bytes / sizeof(Minimizer_Offset_Pair);
        worker.emplace_back(std::sort<Minimizer_Offset_Pair*>, min_buf[id], min_buf[id] + min_buf_size);
    }

    // Wait for the sorting to complete.
    for(uint16_t id = 0; id < producer_count; ++id)
    {
        if(!worker[id].joinable())
        {
            std::cerr << "Early termination encountered for some worker thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        worker[id].join();
    }


    // TODO: multiway-merge.
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
