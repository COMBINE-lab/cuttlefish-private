
#include "Kmer_Index_Utility.hpp"


const std::string Kmer_Index_Utility::config_file_path(const std::string& idx_path)
{
    return idx_path + CONFIG_FILE_EXT;
}


uint16_t Kmer_Index_Utility::kmer_len(const std::string& config_path)
{
    std::ifstream config_file(config_path.c_str(), std::ios::in | std::ios::binary);

    uint16_t k;
    config_file.read(reinterpret_cast<char*>(&k), sizeof(k));

    if(!config_file)
    {
        std::cerr << "Error reading from the configuration file at " << config_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    return k;
}


uint16_t Kmer_Index_Utility::minimizer_len(const std::string& config_path)
{
    std::ifstream config_file(config_path.c_str(), std::ios::in | std::ios::binary);

    config_file.seekg(sizeof(uint16_t), std::ios::beg); // Skip the `k`-value.

    uint16_t l;
    config_file.read(reinterpret_cast<char*>(&l), sizeof(l));

    if(!config_file)
    {
        std::cerr << "Error reading from the configuration file at " << config_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    return l;
}
