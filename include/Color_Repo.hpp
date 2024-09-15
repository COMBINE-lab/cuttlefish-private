
#ifndef COLOR_REPO_HPP
#define COLOR_REPO_HPP



#include "Ext_Mem_Bucket.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <string>
#include <vector>
#include <cassert>


namespace cuttlefish
{

// External-memory repository for color-sets.
class Color_Repo
{
    typedef Ext_Mem_Bucket<source_id_t> bucket_t;

private:

    std::vector<Padded<bucket_t>> B;    // Worker-specific color-buckets.


public:

    // Initializes the color-repository at path-prefix `path`.
    void init(const std::string& path);

    // Returns the appropriate color-bucket for the worker.
    bucket_t& bucket();
};


inline typename Color_Repo::bucket_t& Color_Repo::bucket()
{
    assert(B.size() == parlay::num_workers());
    return B[parlay::worker_id()].unwrap();
}

}



#endif