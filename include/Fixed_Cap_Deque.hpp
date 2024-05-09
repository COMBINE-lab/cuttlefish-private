
#ifndef FIXED_CAP_DEQUE_HPP
#define FIXED_CAP_DEQUE_HPP



#include "utility.hpp"

#include <cstddef>


namespace cuttlefish
{

// Fixed-capacity deque of capacity `cap_`, which needs to be a power of 2. UB
// when size exceeds capacity.
template <typename T_>
class Fixed_Cap_Deque
{
private:

    const std::size_t cap_; // Maximum capacity.
    std::size_t front_; // Current index to the front.
    std::size_t back_;    // Current index to the back.
    T_* const arr;  // Underlying memory pool.

    void grow_back() { back_ = (back_ < cap_ - 1 ? back_ + 1 : 0); }

    void grow_front() { front_ = (front_ > 0 ? front_ - 1 : cap_ - 1); }

    void shrink_back() { back_ = (back_ > 0 ? back_ - 1 : cap_ - 1); }

    void shrink_front() { front_ = (front_ < cap_ - 1 ? front_ + 1 : 0); }

public:

    Fixed_Cap_Deque(const std::size_t capacity):
          cap_(capacity)
        , front_(0)
        , back_(0)
        , arr(allocate<T_>(cap_ + 1))
    {}

    ~Fixed_Cap_Deque() { deallocate(arr); }

    bool empty() const { return front_ == back_; }

    const T_& front() const { return arr[front_]; }

    const T_& back() const { return arr[back_ > 0 ? back_ - 1 : cap_ - 1]; }

    void clear() { front_ = back_; }

    void push_back(const T_& val) { arr[back_] = val; grow_back(); }

    void push_front(const T_& val) { grow_front(); arr[front_] = val; }

    template <typename... Args>
    void emplace_back(Args&&... args) { new (arr + back_) T_(args...); grow_back(); }

    template <typename... Args>
    void emplace_front(Args&&... args) { grow_front(); new (arr + front_) T_(args...); }

    void pop_back() { shrink_back(); }

    void pop_front() { shrink_front(); }
};

}



#endif
