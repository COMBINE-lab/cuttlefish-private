
#ifndef MINIMIZER_INSTANCE_HPP
#define MINIMIZER_INSTANCE_HPP



#include "Min_Heap.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <utility>


// =============================================================================


// A `(minimizer, offset)` pair, where `minimizer` is the l-minimizer (for a
// pre-defined `l`) of some k-mer `x` in some sequence `seq`, where `minimizer`
// is at the index `offset` in `seq`.
class Minimizer_Instance
{
private:

    typedef cuttlefish::minimizer_t minimizer_t;

    minimizer_t minimizer_; // The minimizer.
    std::size_t offset_;    // Offset of the minimizer in the underlying sequence.

public:

    // Constructs an empty minimizer instance.
    Minimizer_Instance()
    {}

    // Constructs an instance of the minimizer `minimizer`, situated at the
    // offset `offset` some sequence.
    Minimizer_Instance(const minimizer_t minimizer, const std::size_t offset):
        minimizer_(minimizer),
        offset_(offset)
    {}

    // Returns the minimizer.
    cuttlefish::minimizer_t minimizer() const { return minimizer_; }

    // Returns the offset.
    std::size_t offset() const { return offset_; }

    // Shifts (to the right) the offset of the minimizer by `offset_shit`.
    void shift(const std::size_t offset_shift)  { offset_ += offset_shift; }

    // Returns `true` iff this instance has a lesser minimizer than `rhs`; and
    // in the case of same minimizers, returns `true` iff it has a lower offset.
    bool operator<(const Minimizer_Instance& rhs) const
    {
        return minimizer_ != rhs.minimizer_ ? (minimizer_ < rhs.minimizer_) : (offset_ < rhs.offset_);
    }
};


// A class to iterate over a collection of type `T_container_` of minimizer
// instances.
template <typename T_container_>
class Minimizer_Instance_Collection_Iterator
{};


// The minimizer instances are stored in a `Minimizer_Instance*`-type container,
// having its size paired.
template <>
class Minimizer_Instance_Collection_Iterator<std::pair<Minimizer_Instance*, std::size_t>>
{
    typedef std::pair<Minimizer_Instance*, std::size_t> T_container_;

private:

    const Minimizer_Instance* const container;  // Pointer to the container.
    const Minimizer_Instance* const container_end;  // Pointer to the container-end.
    const Minimizer_Instance* ptr;  // Pointer to the current minimizer instance in a linear scan over the container.


public:

    Minimizer_Instance_Collection_Iterator(const T_container_& container):
        container(container.first),
        container_end(container.first + container.second),
        ptr(this->container)
    {}

    void operator++()
    {
        ptr++;
    }

    const Minimizer_Instance& operator*() const
    {
        return *ptr;
    }

    explicit operator bool() const
    {
        return ptr != container_end;
    }
};


// A class to multiway-merge a number of sorted minimizer-instance containers
// of type `T_container_` and collect the union instances in sorted order.
template <typename T_container_>
class Minimizer_Instance_Multiway_Merger
{
private:

    const uint32_t source_count;    // Number of source containers.
    std::vector<T_container_>& source;  // Collection of the input containers.
    std::vector<Minimizer_Instance_Collection_Iterator<T_container_>> iterator; // Collection of iterators over the sources.

    typedef std::pair<Minimizer_Instance, uint16_t> Min_Source_Pair;    // Type of elements to be sorted (merged) through the min heap.
    Min_Heap<Min_Source_Pair> min_heap; // Min heap of minimizer instances and their source container IDs.

public:

    // Constructs a multiway merger for the minimizer instance containers
    // collection `source`.
    Minimizer_Instance_Multiway_Merger(std::vector<T_container_>& source);

    // Peeks into the next "minimum" minimizer instance from the multiway
    // merger into `min`. Returns `true` iff an instance could be extracted,
    // i.e. the min heap were not empty.
    bool peek(Minimizer_Instance& min) const;

    // Fetches the next "minimum" minimizer instance from the multiway merger
    // into `min`. Returns `true` iff an instance could be extracted,
    // i.e. the min heap were not empty.
    bool next(Minimizer_Instance& min);
};


template <typename T_container_>
inline Minimizer_Instance_Multiway_Merger<T_container_>::Minimizer_Instance_Multiway_Merger(std::vector<T_container_>& source):
    source_count(source.size()),
    source(source)
{
    for(uint32_t i = 0; i < source_count; ++i)
    {
        iterator.emplace_back(source[i]);
        if(iterator[i])
        {
            // min_heap.emplace(*iterator[i], i);
            min_heap.push(Min_Source_Pair(*iterator[i], i));
            ++iterator[i];
        }
    }
}


template <typename T_container_>
inline bool Minimizer_Instance_Multiway_Merger<T_container_>::next(Minimizer_Instance& min)
{
    if(!peek(min))
        return false;

    const auto source_id = min_heap.top().second;
    min_heap.pop();

    if(iterator[source_id])
    {
        // min_heap.emplace(*iterator[source_id], source_id);
        min_heap.push(Min_Source_Pair(*iterator[source_id], source_id));
        ++iterator[source_id];
    }

    return true;
}


template <typename T_container_>
inline bool Minimizer_Instance_Multiway_Merger<T_container_>::peek(Minimizer_Instance& min) const
{
    if(min_heap.empty())
        return false;

    min = min_heap.top().first;
    return true;
}



#endif
