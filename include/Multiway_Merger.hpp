
#ifndef MULTIWAY_MERGER_HPP
#define MULTIWAY_MERGER_HPP



#include "Kmer_SPSC_Iterator.hpp"
#include "Kmer.hpp"
#include "Min_Heap.hpp"

#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <algorithm>


// =============================================================================


// A class to multiway-merge a number of k-mer databases and collect the list of
// the union k-mers in sorted order.
template <uint16_t k>
class Multiway_Merger
{
private:

    const uint32_t db_count;    // Number of input k-mer databases.
    const std::vector<std::string> db_list; // List of the paths to the input k-mer databases.
    std::vector<Kmer_SPSC_Iterator<k>> iterator;    // Collection of iterators over the input k-mer databases.

    std::unordered_map<std::string, uint32_t> db_id;    // Numerical ID of the input k-mer databases.

    
    typedef uint16_t source_id_t;   // Type of the source-identifier with each k-mer coming from some database.
    typedef std::pair<Kmer<k>, source_id_t> heap_elem_t;    // Type of elements to be sorted (merged) through min heap.
    Min_Heap<heap_elem_t> min_heap; // Min heap of k-mers and their source database IDs.

    std::vector<uint32_t> kmer_count;   // Number of k-mers from each input database currently present in the heap.



public:

    // Constructs a multiway merger for the k-mer databases at the paths-
    // collection `db_list`.
    Multiway_Merger(const std::vector<std::string>& db_list);

    // Destructs the multiway merger.
    ~Multiway_Merger();

    // Copy-construction is prohibited. This is a complex object containing
    // k-mer database iterator collections, each in some *arbitrary* state,
    // which should not be allowed to be copied.
    Multiway_Merger(const Multiway_Merger&) = delete;

    // Copy-assignment is prohibited. This is a complex object containing k-mer
    // database iterator collections, each in some *arbitrary* state, which
    // should not be allowed to be copied.
    const Multiway_Merger& operator=(const Multiway_Merger&) = delete;

    // Launches the multiway merge over the provided database list.
    void launch();
};



#endif
