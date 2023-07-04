
#include "Unitig_File.hpp"


namespace cuttlefish
{

Unitig_File_Writer::Unitig_File_Writer(const std::string& file_path):
      file_path(file_path)
    , total_sz(0)
    , unitig_c(0)
    , output(file_path)
    , output_len(length_file_path())
{}


void Unitig_File_Writer::close()
{
    if(!buf.empty())
        flush_unitigs();

    if(!len.empty())
        flush_lengths();

    output.close();
    output_len.close();
}


Unitig_File_Reader::Unitig_File_Reader(const std::string& file_path):
      file_path(file_path)
    , input(file_path)
    , len(length_file_path().c_str())
    , buf_idx(0)
    , uni_idx_in_file(0)
    , uni_idx_in_mem(0)
    , unitig_count_(len.size())
    , unitig_parsed_(0)
    , total_sz(0)
{}

}
