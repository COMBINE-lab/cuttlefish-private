
#include "Kmer_Index.hpp"
#include "Sparse_Lock.hpp"
#include "Minimizer_Instance_Merger.hpp"
#include "Input_Defaults.hpp"
#include "utility.hpp"

#include <string>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <chrono>


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
void Kmer_Index<k>::index()
{
    typedef std::chrono::high_resolution_clock::time_point time_point_t;
    constexpr auto now = std::chrono::high_resolution_clock::now;
    constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };

    // Dump the remaining in-memory content to disk.
    close_deposit_stream();
    std::cout << "Closed the sequence deposit stream.\n";
    const time_point_t t_deposit = now();

    // Load and sort the minimizer files' content.
    read_and_sort_minimizers();
    const time_point_t t_sort = now();
    std::cout << "Read and sorted the minimizers files. Time taken = " << duration(t_sort - t_deposit) << " seconds.\n";

    // Multiway-merge the minimizer instances.
    merge_minimizers();
    const time_point_t t_merge = now();
    std::cout << "Merged the minimizers files. Time taken = " << duration(t_merge - t_sort) << " seconds.\n";

    // Construct an MPHF over the unique minimizers.
    construct_minimizer_mphf();
    const uint64_t total_bits = min_mphf->totalBitSize();
    std::cout <<    "For the minimizer-MPHF:\n"
                    "\tTotal size: " << total_bits / (8U * 1024U * 1024U) << " MB."
                    " Bits per k-mer: " << static_cast<double>(total_bits) / min_count << ".\n";
    const time_point_t t_mphf = now();
    std::cout << "Constructed the minimizer-MPHF. Time taken = " << duration(t_mphf - t_merge) << " seconds.\n";

    // Count the instances per minimizer.
    count_minimizer_instances();
    const time_point_t t_count = now();
    std::cout << "Counted the minimizer instances. Time taken = " << duration(t_count - t_mphf) << " seconds.\n";

    // Get the offsets of each minimizer.
    get_minimizer_offsets();
    const time_point_t t_off = now();
    std::cout << "Gathered the minimizer instances' offsets. Time taken = " << duration(t_off - t_count) << " seconds.\n";
}


template <uint16_t k>
void Kmer_Index<k>::close_deposit_stream()
{
    // Flush the remaining buffer content, and release memory of all the path-sequence buffers of the producers.
    for(std::size_t id = 0; id < producer_count; ++id)
    {
        if(!producer_path_buf[id].empty())
            flush(id);

        force_free(producer_path_buf[id]);
        force_free(producer_path_end_buf[id]);

        producer_minimizer_file[id].close();
    }

    paths.serialize(cuttlefish::_default::PATH_FILE_EXT);   // TODO: add ext. to o/p file name (from `Build_Params`).

    // Release memory of the concatenated paths.
    sum_paths_len = paths.size();
    force_free(paths);


    // Compact the path endpoints to just as many bits as required.
    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(sum_paths_len)));
    assert(bits_per_entry > 0);
    path_count = path_ends_vec.size();
    path_ends = new path_end_vector_t(bits_per_entry, path_count);
    auto& p_end = *path_ends;
    for(std::size_t i = 0; i < path_count; ++i)
        p_end[i] = path_ends_vec[i];

    force_free(path_ends_vec);

    p_end.serialize(cuttlefish::_default::PATH_ENDS_FILE_EXT); // TODO: add ext. to o/p file name (from `Build_Params`).
    force_free(p_end);
}


template <uint16_t k>
void Kmer_Index<k>::read_and_sort_minimizers()
{
    // TODO: bound memory usage within a provided argument.

    // Launch sorter threads.
    for(uint16_t id = 0; id < producer_count; ++id)
        worker.emplace_back([this](const uint16_t id)
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

                std::sort(min_group[id], min_group[id] + min_group_size[id]);
            },
            id);

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


    // Release memory of the minimizer containers.
    for(uint16_t id = 0; id < producer_count; ++id)
        std::free(min_group[id]);

    force_free(min_group), force_free(min_group_size);
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
    min_instance_count = new min_vector_t(bits_per_entry, min_count + 1);

    std::FILE* const min_file = std::fopen(cuttlefish::_default::MINIMIZER_FILE_EXT, "rb");   // TODO: fix placeholder file name.
    Minimizer_Instance_Iterator<std::FILE*> min_inst_iter(min_file);
    constexpr std::size_t min_buf_sz = 5U * 1024U * 1024U / sizeof(Minimizer_Instance); // Minimizer instance chunk size, in elements: 5 MB total.

    std::vector<std::thread> worker;
    Sparse_Lock<Spin_Lock> idx_lock(min_count + 1, idx_lock_count);
    std::vector<std::size_t> max_count(producer_count);

    for(uint16_t t = 0; t < producer_count; ++t)
        worker.emplace_back(
            [this, &min_inst_iter, &idx_lock](std::size_t& max_c)
            {
                Minimizer_Instance* min_chunk = static_cast<Minimizer_Instance*>(std::malloc(min_buf_sz * sizeof(Minimizer_Instance)));
                std::size_t buf_elem_count;
                uint64_t h;
                std::size_t count;
                auto& mi_count = *min_instance_count;

                max_c = 0;
                while((buf_elem_count = min_inst_iter.next(min_chunk, min_buf_sz)) != 0)
                    for(std::size_t idx = 0; idx < buf_elem_count; ++idx)
                    {
                        h = hash(min_chunk[idx].minimizer());

                        idx_lock.lock(h);
                        count = mi_count[h] = mi_count[h] + 1;
                        idx_lock.unlock(h);

                        if(max_c < count)
                            max_c = count;
                    }

                std::free(min_chunk);
            },
            std::ref(max_count[t]));

    for(uint16_t t = 0; t < producer_count; ++t)
    {
        if(!worker[t].joinable())
        {
            std::cerr << "Early termination encountered for some worker thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        worker[t].join();
        if(max_inst_count < max_count[t])
            max_inst_count = max_count[t];
    }

    std::cout << "Maximum instance count of a minimizer: " << max_inst_count << ".\n";


    // Transform the instance counts to cumulative counts, ordered per the minimizers' hash.
    auto& mi_count = *min_instance_count;
    uint64_t cum_count = (mi_count[0] = 0);
    for(std::size_t i = 1; i <= min_count; ++i)
        mi_count[i] = (cum_count += mi_count[i]);

    const std::string counts_file_path = "min.counts";  // TODO: placeholder for now.
    mi_count.serialize(counts_file_path.c_str());
}


template <uint16_t k>
void Kmer_Index<k>::get_minimizer_offsets()
{
    const uint32_t bits_per_entry = static_cast<uint32_t>(std::ceil(std::log2(sum_paths_len)));
    assert(bits_per_entry > 0);
    min_offset = new min_vector_t(bits_per_entry, num_instances);

    std::FILE* const min_file = std::fopen(cuttlefish::_default::MINIMIZER_FILE_EXT, "rb");   // TODO: fix placeholder file name.
    Minimizer_Instance_Iterator<std::FILE*> min_inst_iter(min_file);
    constexpr std::size_t min_buf_sz = 5U * 1024U * 1024U / sizeof(Minimizer_Instance); // Minimizer instance chunk size, in elements: 5 MB total.

    std::vector<std::thread> worker;
    Sparse_Lock<Spin_Lock> idx_lock(min_count, idx_lock_count);

    for(uint16_t t = 0; t < producer_count; ++t)
        worker.emplace_back(
            [this, &min_inst_iter, &idx_lock]()
            {
                Minimizer_Instance* min_chunk = static_cast<Minimizer_Instance*>(std::malloc(min_buf_sz * sizeof(Minimizer_Instance)));
                std::size_t buf_elem_count;
                uint64_t h;
                std::size_t offset;
                auto& mi_count = *min_instance_count;
                auto& m_offset = *min_offset;

                while((buf_elem_count = min_inst_iter.next(min_chunk, min_buf_sz)) != 0)
                    for(std::size_t idx = 0; idx < buf_elem_count; ++idx)
                    {
                        h = hash(min_chunk[idx].minimizer()) - 1;

                        // Tracking the next offsets to put the minimizers' offsets, within `min_instance_count` in-place,
                        // transforms it by sliding it one index to the left. That is, now, `min_instance_count[MPH(min)]`
                        // has the prefix-sum, opposed to the index `MPH(min) + 1` as earlier.
                        idx_lock.lock(h);
                        offset = mi_count[h];
                        mi_count[h] = offset + 1;
                        idx_lock.unlock(h);

                        m_offset[offset] = min_chunk[idx].offset();
                    }
            }
        );

    for(uint16_t t = 0; t < producer_count; ++t)
    {
        if(!worker[t].joinable())
        {
            std::cerr << "Early termination encountered for some worker thread. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        worker[t].join();
    }


    const std::string offsets_file_path = "min.offsets";    // TODO: placeholder for now.
    min_offset->serialize(offsets_file_path.c_str());


    // Release the portion of the index still in memory.
    delete min_mphf;
    force_free(*min_instance_count);
    force_free((*min_offset));
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
