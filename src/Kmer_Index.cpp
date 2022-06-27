
#include "Kmer_Index.hpp"
#include "Input_Defaults.hpp"

#include <string>


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const uint16_t l, const uint16_t worker_count):
    // TODO: reserve space for `paths`, preferably from additional k-mer count field
    l(l),
    worker_count(worker_count),
    worker_path_buf(worker_count),
    worker_minimizer_buf(worker_count),
    worker_minimizer_file(worker_count),
    curr_token{0}
{
    for(uint16_t id = 0; id < worker_count; ++id)
        worker_minimizer_file[id].open(std::to_string(id) + cuttlefish::_default::MINIMIZER_FILE_EXT, std::ios::out | std::ios::binary);
}


template <uint16_t k>
const typename Kmer_Index<k>::Worker_Token Kmer_Index<k>::get_token()
{
    lock.lock();
    const std::size_t token = curr_token++;
    lock.unlock();

    return Worker_Token(token);
}


template <uint16_t k>
void Kmer_Index<k>::finalize_workers()
{
    for(std::size_t id = 0; id < worker_count; ++id)
    {
        if(!worker_path_buf[id].empty())
            flush(id);

        worker_minimizer_file[id].close();
    }

    std::ofstream path_file(cuttlefish::_default::PATH_FILE_EXT, std::ios::out | std::ios::binary); // TODO: add ext. to o/p file name (from `Build_Params`).
    paths.serialize(path_file);
    // path_file.write(reinterpret_cast<const char*>(paths.data()), paths.size());  // For testing.
    path_file.close();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
