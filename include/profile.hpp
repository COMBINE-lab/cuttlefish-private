
#ifndef PROFILE_HPP
#define PROFILE_HPP



#include <string>
#include <functional>


// Executes `f()` with optional profiling at record-file `tag`.
#ifdef PART_PROFILE
    #define EXECUTE(f, tag) (cuttlefish::profile([&](){ f(); }, (tag)));
#else
    #define EXECUTE(f, tag) (f());
#endif


namespace cuttlefish
{
    // `perf`-profiles the execution of the function `f` at file `record`.
    void profile(const std::function<void()>& f, const std::string& record);
};



#endif
