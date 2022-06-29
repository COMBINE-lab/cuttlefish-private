
#include "Minimizer_Instance_Iterator.hpp"

#include <iostream>


Minimizer_Instance_Iterator<std::FILE*>::Minimizer_Instance_Iterator():
    file_ptr(nullptr),
    pos(0),
    buffer(nullptr),
    buf_elem_count(0),
    buf_idx(0)
{}


Minimizer_Instance_Iterator<std::FILE*>::Minimizer_Instance_Iterator(std::FILE* const file_ptr): Minimizer_Instance_Iterator()
{
    this->file_ptr = file_ptr;

    if(std::fseek(file_ptr, 0, SEEK_SET))
    {
        std::cerr << "Error reading the minimizer files. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    peek(); // To handle the edge case of empty files.
}


Minimizer_Instance_Iterator<std::FILE*>::Minimizer_Instance_Iterator(const Minimizer_Instance_Iterator& other):
    file_ptr(other.file_ptr),
    pos(other.pos),
    buffer(nullptr),
    buf_elem_count(other.buf_elem_count),
    buf_idx(other.buf_idx),
    elem(other.elem)
{
    if(other.buffer != nullptr)
    {
        buffer = static_cast<Minimizer_Instance*>(std::malloc(buf_sz * sizeof(Minimizer_Instance)));
        std::memcpy(buffer, other.buffer, buf_sz * sizeof(Minimizer_Instance));
    }
}


Minimizer_Instance_Iterator<std::FILE*>::~Minimizer_Instance_Iterator()
{
    if(buffer != nullptr)
        std::free(buffer);
}


void Minimizer_Instance_Iterator<std::FILE*>::peek()
{
    const auto next_ch = std::fgetc(file_ptr);
    std::ungetc(next_ch, file_ptr);

    if(next_ch == EOF)
        file_ptr = nullptr;
}
