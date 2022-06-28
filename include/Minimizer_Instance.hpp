
#ifndef MINIMIZER_INSTANCE_HPP
#define MINIMIZER_INSTANCE_HPP



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

    minimizer_t minimizer;  // The minimizer.
    std::size_t offset;     // Offset of the minimizer in the underlying sequence.

public:

    // Constructs an instance of the minimizer `minimizer`, situated at the
    // offset `offset` some sequence.
    Minimizer_Instance(const minimizer_t minimizer, const std::size_t offset):
        minimizer(minimizer),
        offset(offset)
    {}

    // Shifts (to the right) the offset of the minimizer by `offset_shit`.
    void shift(const std::size_t offset_shift)  { offset += offset_shift; }

    // Returns `true` iff this instance has a lesser minimizer than `rhs`; and
    // in the case of same minimizers, returns `true` iff it has a lower offset.
    bool operator<(const Minimizer_Instance& rhs) const
    {
        return minimizer != rhs.minimizer ? (minimizer < rhs.minimizer) : (offset < rhs.offset);
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



#endif
