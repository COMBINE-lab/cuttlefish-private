
#ifndef FASTA_RECORD_HPP
#define FASTA_RECORD_HPP



#include "fmt/format.h"

#include <cstdint>
#include <cstddef>
#include <string>
#include <string_view>


// =============================================================================
// A class wrapping a basic FASTA record. The class is specifically designed for
// writing purposes of output maximal unitigs in the FASTA format.
class FASTA_Record
{
private:

    const fmt::format_int id_;  // Identifier for the FASTA sequence.
    const std::string_view seq_;    // The FASTA sequence.
    const std::string_view seq_add_;    // Additional FASTA sequence (in case the original sequence `*seq` is broken into two parts).
    const std::size_t offset_;  // Offset position into the sequence `seq_`—data before this index will be skipped in the record.
    const std::size_t offset_add_;  // Offset position into the additional sequence `seq_add`—data before this index will be skipped in the record.


    // Constructs a FASTA header with identifier `id`, along with the sequences
    // `seq` and `seq_add` (onward their indices `offset` and `offset_add`,
    // respectively). Only constant references to the sequences are captured,
    // so the record's correctness holds as long as the referred sequences themselves
    // remains unaltered.
    FASTA_Record(uint64_t id, const std::string_view& seq, const std::string_view& seq_add, std::size_t offset = 0, std::size_t offset_add = 0);


public:

    // Constructs a FASTA header with identifier `id` and the sequence `seq`
    // (onward its index `offset`). Only a constant reference to the sequence
    // is captured, so the record's correctness holds as long as the referred
    // sequence itself remains unaltered.
    FASTA_Record(uint64_t id, const std::string& str, std::size_t offset = 0);

    // Constructs a FASTA header with identifier `id` and the sequence `seq`
    // (onward its index `offset`). Only a constant reference to the sequence
    // is captured, so the record's correctness holds as long as the referred
    // sequence itself remains unaltered.
    FASTA_Record(uint64_t id, const std::string_view& str, std::size_t offset = 0);

    // Constructs a FASTA header with identifier `id`, along with the sequences
    // `seq` and `seq_add` (onward their indices `offset` and `offset_add`,
    // respectively). Only constant references to the sequences are captured,
    // so the record's correctness holds as long as the referred sequences themselves
    // remains unaltered.
    FASTA_Record(uint64_t id, const std::string& seq, const std::string& seq_add, std::size_t offset = 0, std::size_t offset_add = 0);

    // Returns the length of the header line of the record.
    std::size_t header_size() const;

    // Returns the length of the sequence of the record.
    std::size_t seq_size() const;

    // Appends the header line to the `buffer`.
    void append_header(std::string& buffer) const;

    // Appends the FASTA sequence to the `buffer`.
    void append_seq(std::string& buffer) const;

    // Appends the FASTA sequence to `buffer` in a rotated form — the
    // added sequence is supposed to be a cycle in a de Bruijn graph `G(·, k)`,
    // and it is right rotated so that the character at index `pivot` is at
    // index 0 finally.
    template <uint16_t k>
    void append_rotated_cycle(std::string& buffer, std::size_t pivot) const;
};


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string& seq, const std::size_t offset):
    FASTA_Record(id, seq, std::string_view(), offset)
{}


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string_view& seq, const std::size_t offset):
    FASTA_Record(id, seq, std::string_view(), offset)
{}


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string& seq, const std::string& seq_add, const std::size_t offset, const std::size_t offset_add):
    FASTA_Record(id, std::string_view(seq), std::string_view(seq_add), offset, offset_add)
{}


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string_view& seq, const std::string_view& seq_add, const std::size_t offset, const std::size_t offset_add):
    id_(id),
    seq_(seq),
    seq_add_(seq_add),
    offset_(offset),
    offset_add_(offset_add)
{}


inline std::size_t FASTA_Record::header_size() const
{
    return  id_.size() + 1; // One additional byte for `>`.
}


inline std::size_t FASTA_Record::seq_size() const
{
    return (seq_.size() - offset_) + (!seq_add_.empty() ? (seq_add_.size() - offset_add_) : 0);
}


inline void FASTA_Record::append_header(std::string& buffer) const
{
    buffer.push_back('>');

    buffer.insert(buffer.end(), id_.data(), id_.data() + id_.size());
}


inline void FASTA_Record::append_seq(std::string& buffer) const
{
    buffer.append(seq_.cbegin() + offset_, seq_.cend());
    if(!seq_add_.empty())
        buffer.append(seq_add_.cbegin() + offset_add_, seq_add_.cend());
}


template <uint16_t k>
inline void FASTA_Record::append_rotated_cycle(std::string& buffer, const std::size_t pivot) const
{
    buffer.append(seq_.cbegin() + pivot, seq_.cend());
    buffer.append(seq_.cbegin() + k - 1, seq_.cbegin() + k - 1 + pivot);
}



#endif
