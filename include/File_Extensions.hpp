
#ifndef FILE_EXTENSIONS_HPP
#define FILE_EXTENSIONS_HPP



namespace cuttlefish
{
    // File extensions for the data structures and files output by the algorithm.
    namespace file_ext
    {
        constexpr char edges_ext[] = ".cf_E";
        constexpr char vertices_ext[] = ".cf_V";
        constexpr char hash_ext[] = ".cf_hf";
        constexpr char buckets_ext[] = ".cf_hb";
        constexpr char unipaths_ext[] = ".fa";
        constexpr char json_ext[] = ".json";
        constexpr char temp[] = ".cf_op";
        
        // For reference dBGs only:

        // TODO: use these to replace the corresponding constants from `CdBG_Writer`.
        constexpr char seg_ext[] = ".cf_seg";
        constexpr char seq_ext[] = ".cf_seq";
    }
}



#endif