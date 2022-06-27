
#include "Kmer_Index.hpp"
#include "Input_Defaults.hpp"

#include <string>


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
        producer_minimizer_file[id].open(std::to_string(id) + cuttlefish::_default::MINIMIZER_FILE_EXT, std::ios::out | std::ios::binary);
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



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
