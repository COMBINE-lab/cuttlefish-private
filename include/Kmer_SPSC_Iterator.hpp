
#ifndef KMER_SPSC_ITERATOR_HPP
#define KMER_SPSC_ITERATOR_HPP



#include "Kmer.hpp"
#include "Kmer_Container.hpp"
#include "Spin_Lock.hpp"
#include "kmc_api/kmc_file.h"

#include <cstdint>
#include <vector>
#include <cstddef>


// =============================================================================


// An "iterator" class to iterate over a k-mer database on disk, where a single
// producer thread (sequentially) reads the raw binary representations of the
// k-mers from disk, and the same thread also acts as the consumer, i.e. it
// fetches and parses the raw k-mers.
// Note: in a technical-C++ sense, it is not an "iterator".
template <uint16_t k>
class Kmer_SPSC_Iterator
{
private:

    const Kmer_Container<k>* const kmer_container;  // The associated k-mer container over which to iterate.
    CKMC_DB kmer_database;  // The k-mer database object.

    const uint64_t kmer_count;  // Number of k-mers present in the underlying database.

    uint64_t kmers_read;    // Number of raw k-mers read off disk by the iterator.

    uint8_t* suff_buf;  // Buffer for the raw binary suffixes of the k-mers.
    std::size_t suf_buf_size = (1lu << 24); // Size of the buffer (in bytes) for raw k-mer suffixes: 16 MB initially, can expand if need be.

    std::vector<std::pair<uint64_t, uint64_t>> pref_buf;    // Buffer for the raw binary prefixes of the k-mers, in the form: `<prefix, #corresponding_suffix>`.
    std::vector<std::pair<uint64_t, uint64_t>>::iterator pref_it;   // Pointer to the prefix to start parsing the remaining unparsed k-mers.

    uint64_t kmers_parsed_off_buf;  // Number of k-mers parsed from the current buffer content.
    uint64_t kmers_available_in_buf;    // Number of k-mers available in the current buffer content.

    // Statuses of the buffer.
    enum class Buffer_Status: uint8_t
    {
        pending,    // currently empty, k-mers yet to be read in from disk;
        available,  // k-mers are available and waiting to be parsed and processed;
        no_more,    // no k-mers will be read in anymore.
    };

    volatile Buffer_Status buf_stat;    // Status of the buffer.
    Spin_Lock buf_lock; // Exclusive-access lock to the buffer.


    // Opens the k-mer database at the path prefix `db_path`.
    void open_kmer_database(const std::string& db_path);

    // Closes the underlying k-mer database.
    void close_kmer_database();

    // Reads raw binary k-mer representations from the underlying database.
    // Reads as much raw data as it fits within the suffix buffer `suff_buf`,
    // and also puts the corresponding prefixes (and their count) in `pref_buf`.
    // Returns `true` iff there were unread data left in the database. If
    // `PREF_ATOMIC` is `true`, then the k-mer content of each prefix (from the
    // prefile-file of the k-mer database) is read into the buffer in its
    // entirety, i.e. if some prefix `pref` has `f` corresponding k-mers, then
    // all the suffixes of those k-mers are read into `suff_buf` at the same
    // time: the suffixes are not read in partial chunks. The buffer is resized
    // if required in such cases.
    template <bool PREF_ATOMIC = false>
    bool read_kmer_data();

    // Sets the buffer status `buf_stat` to `new_stat`, in a thread-safe manner.
    void set_buffer_status(Buffer_Status new_stat);


public:

    // Constructs an iterator for the provided container `kmer_container`.
    Kmer_SPSC_Iterator(const Kmer_Container<k>* kmer_container);

    // Copy-construction is prohibited. This is a complex object with a KMC
    // database object in some *arbitrary* state, which should not be allowed
    // to be copied. Besides, it contains certain constant fields.
    Kmer_SPSC_Iterator(const Kmer_SPMC_Iterator<k>& other) = delete;

    // Destructs the iterator.
    ~Kmer_SPSC_Iterator();

    // Copy-assignment is prohibited. This is a complex object with a KMC
    // database object in some *arbitrary* state, which should not be allowed
    // to be copied. Besides, it contains certain constant fields.
    Kmer_SPSC_Iterator<k>& operator=(const Kmer_SPSC_Iterator<k>& rhs) = delete;

    // Equality operators are prohibited.
    bool operator==(const Kmer_SPSC_Iterator<k>& rhs) = delete;

    // Equality operators are prohibited.
    bool operator!=(const Kmer_SPSC_Iterator<k>& rhs) = delete;

    // Tries to parse the next k-mer from the buffer. Returns `true` iff k-mers
    // were remaining, i.e. the underlying database has not been depleted.
    bool parse_kmer(Kmer<k>& kmer);

    // Launches the iterator.
    void launch();

    // Closes the iterator.
    void close();
};


template <uint16_t k>
inline Kmer_SPSC_Iterator<k>::Kmer_SPSC_Iterator(const Kmer_Container<k>* const kmer_container):
    kmer_container{kmer_container},
    kmer_count{kmer_container->size()},
    kmers_read{0},
    suff_buf{nullptr},
    kmers_parsed_off_buf{0},
    kmers_available_in_buf{0},
    buf_stat{Buffer_Status::pending}
{}


template <uint16_t k>
inline Kmer_SPSC_Iterator<k>::~Kmer_SPSC_Iterator()
{
    if(suff_buf != nullptr)
        delete[] suff_buf;

    pref_buf.clear();
    pref_buf.shrink_to_fit();
}


template <uint16_t k>
inline void Kmer_SPSC_Iterator<k>::open_kmer_database(const std::string& db_path)
{
    if(!kmer_database.open_for_cuttlefish_listing(db_path))
    {
        std::cerr << "\nError opening k-mer database with path prefix " << db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline void Kmer_SPSC_Iterator<k>::close_kmer_database()
{
    if(!kmer_database.Close())
    {
        std::cerr << "\nError closing k-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline bool Kmer_SPSC_Iterator<k>::parse_kmer(Kmer<k>& kmer)
{
    if(kmers_parsed_off_buf == kmers_available_in_buf)
    {
        buf_lock.lock();

        bool data_read = false;

        if(buf_stat == Buffer_Status::available)
            data_read = true;
        else if(buf_stat == Buffer_Status::pending)
            data_read = read_kmer_data();

        buf_lock.unlock();

        if(!data_read)
            return false;
    }


    kmer_database.parse_kmer_buf<k>(pref_it, suff_buf, kmers_parsed_off_buf * kmer_database.suff_record_size(), kmer);
    kmers_parsed_off_buf++;

    if(kmers_parsed_off_buf == kmers_available_in_buf)
        set_buffer_status(Buffer_Status::pending);


    return true;
}


template <uint16_t k>
template <bool PREF_ATOMIC>
inline bool Kmer_SPSC_Iterator<k>::read_kmer_data()
{
    if(kmer_database.Eof())
    {
        buf_stat = Buffer_Status::no_more;  // No locks being acquired for the buffer—since the buffer is already in use,
                                            // it is assumed that the lock has been acquired by the invoker.
        return false;
    }


    kmers_available_in_buf =    ([&]() constexpr -> uint64_t
                                { 
                                    if constexpr (PREF_ATOMIC)
                                        return kmer_database.read_raw_suffixes_atomic(suff_buf, pref_buf, suf_buf_size);
                                    
                                    return kmer_database.read_raw_suffixes(suff_buf, pref_buf, suf_buf_size);
                                }());
    pref_it = pref_buf.begin();

    if(!kmers_available_in_buf)
    {
        std::cerr << "\nError reading the k-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    kmers_read += kmers_available_in_buf;
    kmers_parsed_off_buf = 0;
    buf_stat = Buffer_Status::available;    // See note above regarding the buffer lock.

    return true;
}


template <uint16_t k>
inline void Kmer_SPSC_Iterator<k>::set_buffer_status(const Buffer_Status new_stat)
{
    buf_lock.lock();
    
    buf_stat = new_stat;

    buf_lock.unlock();
}


template <uint16_t k>
inline void Kmer_SPSC_Iterator<k>::launch()
{
    open_kmer_database(kmer_container->container_location());

    suff_buf = new uint8_t[suf_buf_size];
    pref_buf.clear();
    set_buffer_status(Buffer_Status::pending);
}


template <uint16_t k>
inline void Kmer_SPSC_Iterator<k>::close()
{
    close_kmer_database();
}



#endif
