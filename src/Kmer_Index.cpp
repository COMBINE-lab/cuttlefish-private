
#include "Kmer_Index.hpp"
#include "Minimizer_Instance_Merger.hpp"
#include "Input_Defaults.hpp"
#include "utility.hpp"

#include <string>
#include <algorithm>
#include <cmath>
#include <cassert>


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const uint16_t l, const uint16_t producer_count):
    // TODO: reserve space for `paths`, preferably from additional k-mer count field
    l(l),
    producer_count(producer_count),
    num_instances(0),
    min_count(0),
    max_inst_count(0),
    producer_path_buf(producer_count),
    producer_path_end_buf(producer_count),
    producer_minimizer_buf(producer_count),
    producer_minimizer_file(producer_count),
    min_group(producer_count, nullptr),
    min_group_size(producer_count),
    min_mphf(nullptr),
    min_instance_count(nullptr),
    min_offset(nullptr),
    curr_token{0}
{
    assert(l <= 32);
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

    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(paths.size())));
    path_ends = new index_vector_t(bits_per_entry, path_ends_vec.size());
    for(std::size_t i = 0; i < path_ends_vec.size(); ++i)
        path_ends->at(i) = path_ends_vec[i];

    std::ofstream path_ends_file(cuttlefish::_default::PATH_ENDS_FILE_EXT, std::ios::out | std::ios::binary); // TODO: add ext. to o/p file name (from `Build_Params`).
    path_ends->serialize(path_ends_file);
    path_ends_file.close();
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

    // Construct an MPHF over the unique minimizers.
    construct_minimizer_mphf();
    const uint64_t total_bits = min_mphf->totalBitSize();
    std::cout <<    "Constructed the minimizer-MPHF.\n"
                    "Total size: " << total_bits / (8 * 1024) << " KB."
                    " Bits per k-mer: " << static_cast<double>(total_bits) / min_count << ".\n";

    // Count the instances per minimizer.
    count_minimizer_instances();
    std::cout << "Gathered the minimizers' instance-counts.\n";

    // Get the offsets of each minimizer.
    get_minimizer_offsets();
    std::cout << "Gathered the minimizers' offsets.\n";
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

    Minimizer_Instance_Merger multiway_merger(min_container);
    Minimizer_Instance min;

    std::ofstream min_file(cuttlefish::_default::MINIMIZER_FILE_EXT, std::ios::out | std::ios::binary); // TODO: fix placeholder file name.

    multiway_merger.peek(min);
    minimizer_t last_min = min.minimizer();
    min_count = 1;
    while(multiway_merger.next(min))
    {
        merged_min_buf.push_back(min);
        if(merged_min_buf.size() >= max_buf_sz)
            dump(merged_min_buf, min_file);

        if(last_min != min.minimizer())
            min_count++,
            last_min = min.minimizer();
    }

    std::cout << "Minimizer instance count: " << num_instances << "\n";
    std::cout << "Unique minimizer count:   " << min_count << "\n";

    if(!merged_min_buf.empty())
        dump(merged_min_buf, min_file);

    min_file.close();
}


template <uint16_t k>
void Kmer_Index<k>::construct_minimizer_mphf()
{
    typedef Minimizer_Instance_Iterator<std::FILE*> min_iter_t;
    std::FILE* const min_file = std::fopen(cuttlefish::_default::MINIMIZER_FILE_EXT, "rb");   // TODO: fix placeholder file name.
    const auto data_iterator = boomphf::range(min_iter_t(min_file), min_iter_t(nullptr));

    const char* const working_dir_path = ".";   // TODO: placeholder for now.
    min_mphf = new boomphf::mphf<minimizer_t, minimizer_hasher_t, false>(min_count, data_iterator, working_dir_path, producer_count, gamma);


    const std::string mphf_file_path = "min.mphf";  // TODO: placeholder for now.
    std::ofstream output(mphf_file_path.c_str(), std::ofstream::out);
    if(output.fail())
    {
        std::cerr << "Error writing to file " << mphf_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    min_mphf->save(output);
    output.close();

    if(std::fclose(min_file))
    {
        std::cout << "Error closing the minimizer file. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
void Kmer_Index<k>::count_minimizer_instances()
{
    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(num_instances)));
    assert(bits_per_entry > 0);
    min_instance_count = new index_vector_t(bits_per_entry, min_count + 1);

    std::FILE* const min_file = std::fopen(cuttlefish::_default::MINIMIZER_FILE_EXT, "rb");   // TODO: fix placeholder file name.
    Minimizer_Instance_Iterator<FILE*> min_inst_iter(min_file);
    minimizer_t min;
    std::size_t count;

    while(min_inst_iter.next(min, count))
    {
        min_instance_count->at(hash(min)) = count;
        if(max_inst_count < count)
            max_inst_count = count;
    }

    std::cout << "Maximum instance count of a minimizer: " << max_inst_count << ".\n";

    uint64_t cum_count = (min_instance_count->at(0) = 0);
    for(std::size_t i = 1; i <= min_count; ++i)
        min_instance_count->at(i) = (cum_count += min_instance_count->at(i));

    const std::string counts_file_path = "min.counts";  // TODO: placeholder for now.
    std::ofstream counts_file(counts_file_path.c_str(), std::ios::out | std::ios::binary);
    min_instance_count->serialize(counts_file);
    counts_file.close();
}


template <uint16_t k>
void Kmer_Index<k>::get_minimizer_offsets()
{
    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(paths.size())));
    assert(bits_per_entry > 0);
    min_offset = new index_vector_t(bits_per_entry, num_instances);

    std::FILE* const min_file = std::fopen(cuttlefish::_default::MINIMIZER_FILE_EXT, "rb");   // TODO: fix placeholder file name.
    Minimizer_Instance_Iterator<FILE*> min_inst_iter(min_file);
    minimizer_t min;
    std::vector<std::size_t> offsets;
    offsets.reserve(max_inst_count);

    while(min_inst_iter.next(min, offsets))
    {
        const std::size_t meta_idx = min_instance_count->at(hash(min) - 1);
        for(std::size_t i = 0; i < offsets.size(); ++i)
            min_offset->at(meta_idx + i) = offsets[i];
    }


    const std::string offsets_file_path = "min.offsets";    // TODO: placeholder for now.
    std::ofstream offsets_file(offsets_file_path.c_str(), std::ios::out | std::ios::binary);
    min_offset->serialize(offsets_file);
    offsets_file.close();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
