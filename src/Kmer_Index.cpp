
#include "Kmer_Index.hpp"


template <uint16_t k>
Kmer_Index<k>::Kmer_Index(const uint16_t worker_count):
    worker_count(worker_count),
    worker_path_buf(worker_count),
    worker_minimizer_buf(worker_count),
    curr_token{0}
{}


template <uint16_t k>
const typename Kmer_Index<k>::Worker_Token Kmer_Index<k>::get_token()
{
    lock.lock();
    const std::size_t token = curr_token++;
    lock.unlock();

    return Worker_Token(token);
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Index)
