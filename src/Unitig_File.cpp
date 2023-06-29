
#include "Unitig_File.hpp"


namespace cuttlefish
{

Unitig_File_Writer::Unitig_File_Writer(const std::string& file_path):
      file_path(file_path)
    , total_sz(0)
    , unitig_c(0)
    , output(file_path)
    , output_off(offset_file_path())
{}


void Unitig_File_Writer::close()
{
    if(!buf.empty())
        flush_unitigs();

    if(!offset.empty())
        flush_offsets();

    output.close();
    output_off.close();
}

}
