
#ifndef UTILITY_HPP
#define UTILITY_HPP



#include <cstddef>
#include <string>
#include <vector>
#include <type_traits>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <chrono>

// TODO: wrap everything here in some namespaces.
// =============================================================================

// Returns a random string of length `len`, using characters from `alphabet`.
std::string get_random_string(size_t len, const char* alphabet =    "0123456789"
                                                                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                    "abcdefghijklmnopqrstuvwxyz");

// Returns `true` iff `pref` is a prefix of `s`.
bool is_prefix(const std::string& s, const std::string& pref);

// Returns `true` iff there exists a file in the file system with the path
// `file_path`.
bool file_exists(const std::string& file_path);

// Returns `true` iff these exists a directory in the file system with the
// path `dir_path`.
bool dir_exists(const std::string& dir_path);

// Returns the file size is bytes of the file at path `file_path`. Returns
// `0` in case the file does not exist.
std::size_t file_size(const std::string& file_path);

// Returns `true` iff there exists some file in the file system path
// `path` with its name being prefixed by `prefix`.
bool file_prefix_exists(const std::string& path, const std::string& prefix);

// Returns a string that is a copy of `s` but has all the whitespaces removed.
std::string remove_whitespaces(const char* s);

// Given the collection of strings `s`, returns the concatenated string
// `s0 : s1 : ... : s_m`, where successive strings are separated by `delimiter`.
const std::string concat_strings(const std::vector<std::string>& s, const std::string& delimiter = ", ");

// Removes the file at path `file_path` from disk. Returns `true` iff the
// removal is successful.
bool remove_file(const std::string& file_path);

// Clears the content of the file at path `file_path`.
void clear_file(const std::string& file_path);

// Returns the name of the file present at the path `file_path`.
const std::string filename(const std::string& file_path);

// Returns the directory of the file present at the path `file_path`.
const std::string dirname(const std::string& file_path);

// Moves the file present at path `from_path` to the path `to_path`.
void move_file(const std::string& from_path, const std::string& to_path);

// Returns the maximum memory ("high-water-mark") used by the running
// process in bytes. Returns `0` in case of errors encountered.
std::size_t process_peak_memory();

// Force-frees the memory allocated to the container `container`.
template <typename T_container_>
void force_free(T_container_& container)
{
    T_container_().swap(container);
}

// Returns pointer to a memory-allocation for `size` elements of type `T_`.
template <typename T_>
T_* allocate(const std::size_t size)
{
    return static_cast<T_*>(std::malloc(size * sizeof(T_)));
}

// Deallocates the pointer `ptr`, allocated with `allocate`.
template <typename T_>
void deallocate(T_* const ptr)
{
    std::free(ptr);
}

// Resizes the type-`T_` container `container` geometrically with the growth
// factor `gf` such that it has a size of at least `sz`.
template <typename T_>
void resize_geometric(T_& container, const std::size_t sz, const double gf = 2.0)
{
    assert(gf > 1.0);

    std::size_t curr_sz = std::max(container.size(), 1lu);
    while(curr_sz < sz)
        curr_sz *= gf;

    if(container.size() < curr_sz)
        container.resize(curr_sz);
}

// Returns the corresponding integer value of type `T_to_` for a type-`T_`
// enum-value `enum_val`.
template <typename T_, typename T_to_ = std::size_t>
constexpr T_to_ as_int(const T_ enum_val)
{
    static_assert(std::is_enum_v<T_>);
    return static_cast<T_to_>(enum_val);
}

// TODO: Add thread-joiner wrapper.


// Wrapper class for a data-element of type `T_` to ensure that in a linear
// collection of `T_`'s, each element is aligned to a cache-line boundary.
template <typename T_>
class alignas(2 * L1_CACHE_LINE_SIZE)
    Padded_Data
{
private:

    T_ data_;


public:

    Padded_Data()
    {}

    Padded_Data(const T_& data):
      data_(data)
    {}

    Padded_Data(T_&& data):
        data_(std::move(data))
    {}

    T_& data() { return data_; }

    const T_& data() const { return data_; }
};


// Fixed-sized deque; UB when size exceeds capacity.
template <typename T_>
class deque
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

    deque(const std::size_t capacity):
          cap_(capacity)
        , front_(0)
        , back_(0)
        , arr(static_cast<T_*>(std::malloc((cap_ + 1) * sizeof(T_))))
    {}

    ~deque() { std::free(arr); }

    bool empty() const { return front_ == back_; }

    const T_& front() const { return arr[front_]; }

    const T_& back() const { return arr[back_ > 0 ? back_ - 1 : cap_ - 1]; }

    void clear() { front_ = back_ = 0; }

    void push_back(const T_& val) { arr[back_] = val; grow_back(); }

    void push_front(const T_& val) { grow_front(); arr[front_] = val; }

    template <typename... Args>
    void emplace_back(Args&&... args) { new (arr + back_) T_(std::forward<Args>(args)...); grow_back(); }

    template <typename... Args>
    void emplace_front(Args&&... args) { grow_front(); new (arr + front_) T_(std::forward<Args>(args)...); }

    void pop_back() { shrink_back(); }

    void pop_front() { shrink_front(); }
};


namespace timer
{
    inline auto now() { return std::chrono::high_resolution_clock::now(); }

    inline auto duration(const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };
}



#endif
