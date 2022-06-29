
#ifndef MINIMIZER_INSTANCE_ITERATOR_HPP
#define MINIMIZER_INSTANCE_ITERATOR_HPP



#include "Minimizer_Instance.hpp"

#include <cstddef>
#include <utility>

// =============================================================================


// A class to iterate over a collection of type `T_container_` of minimizer
// instances.
template <typename T_container_>
class Minimizer_Instance_Iterator
{};


// The minimizer instances are stored in a `Minimizer_Instance*`-type container,
// having its size paired.
template <>
class Minimizer_Instance_Iterator<std::pair<Minimizer_Instance*, std::size_t>>
{
    typedef std::pair<Minimizer_Instance*, std::size_t> T_container_;

private:

    const Minimizer_Instance* const container;  // Pointer to the container.
    const Minimizer_Instance* const container_end;  // Pointer to the container-end.
    const Minimizer_Instance* ptr;  // Pointer to the current minimizer instance in a linear scan over the container.


public:

    Minimizer_Instance_Iterator(const T_container_& container):
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
