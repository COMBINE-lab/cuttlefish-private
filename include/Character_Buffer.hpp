
#ifndef CHARACTER_BUFFER_HPP
#define CHARACTER_BUFFER_HPP



#include "Spin_Lock.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "FASTA_Record.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>


// A buffer class to contain contiguous characters. The buffer is to have a maximum
// capacity of `CAPACITY` (although it is non-binding when a string with length 
// larger than that is added), and it flushes to a sink of type `T_sink_` when it
// overflows or is destructed. Writing to the provided sink (in the constructor)
// is thread-safe. By default, `CAPACITY = 100 KB` (soft limit) worth of characters
// can be retained in memory, at most, before flushing.
template <typename T_sink_, std::size_t CAPACITY = 100 * 1024ULL>
class Character_Buffer
{
private:

    // TODO: replace with a custom-container.
    std::vector<char> buffer;   // The character buffer.
    T_sink_& sink;  // Reference to the sink to flush the buffer content to.


    // Ensures that `buffer` has enough space for additional `append_size`
    // number of bytes, using flush and allocation as necessary.
    void ensure_space(std::size_t append_size);

    // Flushes the buffer content to the sink, and clears the buffer.
    void flush();


public:

    // Constructs a character buffer object that would flush its content to `sink`.
    Character_Buffer(T_sink_& sink);

    Character_Buffer(const Character_Buffer& rhs) = default;

    Character_Buffer(Character_Buffer&& rhs) = default;

    // Appends the content of `str` to the buffer. Flushes are possible.
    template <typename T_container_>
    void operator+=(const T_container_& str);

    // Appends the content of the FASTA record `fasta_rec` to the buffer. Flushes
    // are possible.
    template <typename T_container_>
    void operator+=(const FASTA_Record<T_container_>& fasta_rec);

    // Appends the content of the FASTA record `fasta_cycle` to the buffer, that is
    // supposed to be a cycle in a de Bruijn graph `G(·, k)`. The cyclic FASTA
    // sequence is rotated around its index `pivot` — the entire sequence is
    // right-rotated so that the `pivot`-index character is at index 0 finally. A
    // line-break is added at the end of the sequence, since the user might not be
    // able to provide it with the "to be rotated" sequence.
    template <uint16_t k, typename T_container_>
    void rotate_append_cycle(const FASTA_Record<T_container_>& fasta_cycle, std::size_t pivot);

    // Returns the `len`-length suffix of the buffer.
    const char* suffix(std::size_t len) const;

    // Flushes the buffer if not empty.
    void close();

    // Destructs the buffer object, flushing it if content are present.
    ~Character_Buffer();
};


// Helper class to actually flush the content of the `Character_Buffer` class to its
// sink of type `T_sink`.
// It's used to circumvent the C++ constraint that partial specialization of a
// a member function is not possible without partially specializing the entire
// class. We need to specialize the actual flushing mechanism to support various
// types of sinks, e.g. `std::ofstream`, `spdlog::logger` etc.
// Since the sole purpose of the class is to support the `Character_Buffer` class
// circumvent some contraint, everything is encapsulated in its specializations
// as private, with `Character_Buffer` as friend.
template <typename T_sink_>
class Character_Buffer_Flusher
{};


template <>
class Character_Buffer_Flusher<std::ofstream>
{
    template <typename, std::size_t> friend class Character_Buffer;

private:

    // Mutual-exclusion lock to control multi-threaded access to otherwise not thread-
    // safe sinks (e.g. `std::ofstream`). Note that, the lock is per sink-type, not per
    // actual sink — which is a limitation.
    static Spin_Lock lock;


    // Writes the content of the vector `buf` to the sink `sink`.
    static void write(std::vector<char>& buf, std::ofstream& sink);
};


template <>
class Character_Buffer_Flusher<Async_Logger_Wrapper>
{
    template <typename, std::size_t> friend class Character_Buffer;

private:

    // Writes the content of the vector `buf` to the sink `sink`. Note that the vector
    // `buf` is modified in the process — a null-terminator (`\0`) is appended at the
    // end — which is expected to be not problematic under the assumption that the
    // buffer is cleared after the write (i.e. flush).
    static void write(std::vector<char>& buf, const Async_Logger_Wrapper& sink);
};


template <typename T_sink_, std::size_t CAPACITY>
inline Character_Buffer<T_sink_, CAPACITY>::Character_Buffer(T_sink_& sink):
    sink(sink)
{
    buffer.reserve(CAPACITY);
}


template <typename T_sink_, std::size_t CAPACITY>
template <typename T_container_>
inline void Character_Buffer<T_sink_, CAPACITY>::operator+=(const T_container_& str)
{
    ensure_space(str.size());

    // `std::memcpy` at the end of `buffer` does not update the size of the vector `buffer`.
    buffer.insert(buffer.end(), str.begin(), str.end());
}


template <typename T_sink_, std::size_t CAPACITY>
template <typename T_container_>
inline void Character_Buffer<T_sink_, CAPACITY>::operator+=(const FASTA_Record<T_container_>& fasta_rec)
{
    ensure_space(fasta_rec.header_size() + 1 + fasta_rec.seq_size() + 1);   // Two extra bytes for the line-breaks.

    fasta_rec.append_header(buffer);    // Append the header.
    buffer.emplace_back('\n');  // Break line.
    fasta_rec.append_seq(buffer);   // Append the sequence.
    buffer.emplace_back('\n');  // Break line.
}


template <typename T_sink_, std::size_t CAPACITY>
template <uint16_t k, typename T_container_>
inline void Character_Buffer<T_sink_, CAPACITY>::rotate_append_cycle(const FASTA_Record<T_container_>& fasta_rec, const std::size_t pivot)
{
    ensure_space(fasta_rec.header_size() + 1 + fasta_rec.seq_size() + 1);   // Two extra bytes for two line-breaks.

    fasta_rec.append_header(buffer);    // Append the header.
    buffer.emplace_back('\n');  // Break line.
    fasta_rec.template append_rotated_cycle<k>(buffer, pivot);  // Append the sequence right-rotated around index `pivot`.
    buffer.emplace_back('\n');  // End the sequence.
}


template <typename T_sink_, std::size_t CAPACITY>
inline const char* Character_Buffer<T_sink_, CAPACITY>::suffix(const std::size_t len) const
{
    return buffer.data() + (buffer.size() - len);
}


template <typename T_sink_, std::size_t CAPACITY>
inline void Character_Buffer<T_sink_, CAPACITY>::ensure_space(const std::size_t append_size)
{
    if(buffer.size() + append_size >= CAPACITY) // Using `>=` since for async logging, a `\0` is inserted at the end of `buffer`.
    {
        flush();
        
        if(append_size >= CAPACITY)
        {
            // std::cerr <<    "A single output string overflows the string-buffer capacity.\n"
            //                 "Output string length: " << str.size() << ", string-buffer capacity: " << CAPACITY << ".\n"
            //                 "Please consider increasing the buffer capacity parameter in build for future use.\n";
            
            buffer.reserve(append_size);    // TODO: fix it—thrashing possible with adversarial `append_size` sequence.
        }
    }
}


template <typename T_sink_, std::size_t CAPACITY>
inline void Character_Buffer<T_sink_, CAPACITY>::flush()
{
    Character_Buffer_Flusher<T_sink_>::write(buffer, sink);

    buffer.clear();
}


template <typename T_sink_, std::size_t CAPACITY>
inline void Character_Buffer<T_sink_, CAPACITY>::close()
{
    if(!buffer.empty())
        flush();
}


template <typename T_sink_, std::size_t CAPACITY>
inline Character_Buffer<T_sink_, CAPACITY>::~Character_Buffer()
{
    close();
}


inline void Character_Buffer_Flusher<std::ofstream>::write(std::vector<char>& buf, std::ofstream& output)
{
    lock.lock();

    output.write(buf.data(), buf.size());
    if(!output)
    {
        std::cerr << "Error writing the output. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    lock.unlock();
}


inline void Character_Buffer_Flusher<Async_Logger_Wrapper>::write(std::vector<char>& buf, const Async_Logger_Wrapper& sink)
{
    buf.emplace_back('\0');

    sink.write(buf.data());
}


#endif
