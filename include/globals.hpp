
#ifndef GLOBALS_HPP
#define GLOBALS_HPP



#include "boost/preprocessor/repetition/repeat.hpp"

#include <cstdint>


// The macro `INSTANCE_COUNT` must be set exactly to `(MAX_K + 1) / 2` for a required maximum k-value.
// Also, the `MAX_K` value must be odd (as the k-values used in the algorithm) for correct results.
// TODO: use more user-friendly `MAX_K` definition as compile flag, and set `INSTANCE_COUNT` internally.
#ifndef INSTANCE_COUNT
    #ifndef FIXED_K
        #define INSTANCE_COUNT 32
    #else
        #define INSTANCE_COUNT ((FIXED_K + 1) / 2)
    #endif
#endif


// Forward declarations of the DNA code types.
namespace DNA
{
    enum Base: uint8_t;

    enum class Extended_Base: uint8_t;
}

// Forward declaration of the k-mer type.
template <uint16_t k> class Kmer;


namespace cuttlefish
{
    constexpr uint16_t MAX_K = (2 * INSTANCE_COUNT - 1);
    constexpr uint16_t MAX_L = 32;


    typedef bool dir_t;
    typedef DNA::Base base_t;
    typedef DNA::Extended_Base edge_encoding_t;
    typedef uint8_t state_code_t;


    constexpr dir_t FWD = true;
    constexpr dir_t BWD = false;


    enum class State_Class: uint8_t
    {
        single_in_single_out = 0,
        multi_in_single_out = 1,
        single_in_multi_out = 2,
        multi_in_multi_out = 3,
    };


    typedef enum class Side: uint8_t
    {
        front = 0,
        back = 1,
        unspecified = 2
    } side_t;

    constexpr auto inv_side = [](const side_t s) { return s == side_t::back ?   side_t::front :
                                                                                (s == side_t::front ? side_t::back : side_t::unspecified); };



    constexpr uint8_t BITS_PER_REF_KMER = 5;
    constexpr uint8_t BITS_PER_READ_KMER = 6;

    // Minimizers can be represented using 64-bit integers.
    typedef uint64_t minimizer_t;


    // YACC-specifics:

    // Type of weights of edges in the discontinuity-graph.
    typedef uint16_t weight_t;

    // Type of the ID of a maximal unitig.
    template <uint16_t k> using max_unitig_id_t = Kmer<k>;

    // Type of the index of a unitig in a bucket.
    typedef uint32_t uni_idx_t;

    // TODO: use `u16` after testing done with `u32`.
    // typedef uint16_t uni_len_t; // Type of the length of a unitig in a bucket.
    typedef uint32_t uni_len_t; // Type of the length of a unitig in a bucket.
}


// Metaprogramming macro-loops for instantiating required template instances.

// Given some `x`, explicitly instantiates the class `class_name` for the template parameter `k` with `2x + 1`;
// i.e. it is an instantiator for odd k-values.
#define INSTANTIATE(z, x, class_name) template class class_name<2 * x + 1>;

// Enumerates all the explicit instantiations of the template class `class_name` using `instantiator`, for all
// `x` in `[0, count)`. The `x`-value is used as appropriate by `instantiator`.
#ifndef FIXED_K
    #define ENUMERATE(count, instantiator, class_name) BOOST_PP_REPEAT(count, instantiator, class_name)
#else
    #define ENUMERATE(count, instantiator, class_name)  instantiator(0, (count - 1), class_name)
#endif

// Given some `x`, explicitly instantiates two instances of the class `class_name`, with the template parameters
// `k` = `2x + 1`, and `BITS_PER_KEY` with `BITS_PER_REF_KMER` and `BITS_PER_READ_KMER` for alternate instances;
// i.e. it is an instantiator for odd k-values and all the different bits-per-key requirements.
#define INSTANTIATE_PER_BIT(z, k, class_name)   template class class_name<2 * k + 1, cuttlefish::BITS_PER_REF_KMER>;\
                                                template class class_name<2 * k + 1, cuttlefish::BITS_PER_READ_KMER>;

// Given some `x`, explicitly instantiates two instances of the class `class_name`, using the template parameter `k`
// with `2x + 1` and `2x + 2`, i.e. it is an instantiator for both odd and even k-values.
#define INSTANTIATE_ALL(z, x, class_name)   template class class_name<2 * x + 1>;\
                                            template class class_name<2 * x + 2>;


// BOOST_PP_REPEAT reference: https://www.boost.org/doc/libs/1_55_0/libs/preprocessor/doc/ref/repeat.html



#endif
