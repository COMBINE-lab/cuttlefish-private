
#include "Kmer_Index.hpp"
#include "Input_Defaults.hpp"
#include "utility.hpp"

#include <string>
#include <algorithm>


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const uint16_t l, const uint16_t producer_count):
    // TODO: reserve space for `paths`, preferably from additional k-mer count field
    l(l),
    producer_count(producer_count),
    producer_path_buf(producer_count),
    producer_minimizer_buf(producer_count),
    producer_minimizer_file(producer_count),
    min_group(producer_count, nullptr),
    min_group_size(producer_count),
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
    // Load and sort the minimizer files' content.
    read_and_sort_minimizers();
    std::cout << "Sorted the minimizers files.\n";

    // Multiway-merge the minimizer instances.
    merge_minimizers();
    std::cout << "Merged the minimizers files.\n";
}


template <uint16_t k>
void Kmer_Index<k>::read_and_sort_minimizers()
{
    // TODO: bound memory usage within a provided argument.

    // Launch sorter threads.
    for(uint16_t id = 0; id < producer_count; ++id)
    {
        const auto min_buf_bytes = file_size(minimizer_file_path(id));   // What happens if it's zero?
        min_group[id] = static_cast<Minimizer_Instance*>(std::malloc(min_buf_bytes));

        std::ifstream input(minimizer_file_path(id).c_str(), std::ios::in | std::ios::binary);
        if(!input.read(reinterpret_cast<char*>(min_group[id]), min_buf_bytes))
        {
            std::cerr << "Error reading the minimizer files. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        min_group_size[id] = min_buf_bytes / sizeof(Minimizer_Instance);
        worker.emplace_back(std::sort<Minimizer_Instance*>, min_group[id], min_group[id] + min_group_size[id]);
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
}


template <uint16_t k>
void Kmer_Index<k>::merge_minimizers()
{
    typedef std::pair<Minimizer_Instance*, std::size_t> min_container_t;

    // Gather the minimizer containers.
    std::vector<min_container_t> min_container;
    for(uint16_t id = 0; id < producer_count; ++id)
        min_container.emplace_back(min_group[id], min_group_size[id]);

    // Multiway merge the minimizer instance containers.

    std::vector<Minimizer_Instance> merged_min_buf;
    constexpr std::size_t max_buf_sz = buf_sz_th / sizeof(Minimizer_Instance);
    merged_min_buf.reserve(max_buf_sz);

    Minimizer_Instance_Multiway_Merger multiway_merger(min_container);
    Minimizer_Instance min;
    uint64_t min_inst_count = 0;

    std::ofstream min_file(cuttlefish::_default::MINIMIZER_FILE_EXT, std::ios::out | std::ios::binary);

    multiway_merger.peek(min);
    minimizer_t last_min = min.minimizer();
    uint64_t min_uniq_count = 1;
    while(multiway_merger.next(min))
    {
        merged_min_buf.push_back(min);
        if(merged_min_buf.size() >= max_buf_sz)
            dump(merged_min_buf, min_file);

        min_inst_count++;

        if(last_min != min.minimizer())
            min_uniq_count++,
            last_min = min.minimizer();
    }

    std::cout << "Minimizer instance count: " << min_inst_count << "\n";
    std::cout << "Unique minimizer count:   " << min_uniq_count << "\n";

    if(!merged_min_buf.empty())
        dump(merged_min_buf, min_file);

    min_file.close();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
