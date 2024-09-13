
#include "Color_Repo.hpp"
#include "parlay/parallel.h"

#include <cstdint>


namespace cuttlefish
{

void Color_Repo::init(const std::string& path)
{
    B.reserve(parlay::num_workers());
    for(uint32_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
        B.emplace_back(path + "." + std::to_string(w_id));
}

}
