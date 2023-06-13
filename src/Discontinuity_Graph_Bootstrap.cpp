
#include "Discontinuity_Graph_Bootstrap.hpp"
#include "Kmer.hpp"
#include "Minimizer_Iterator.hpp"
#include "Discontinuity_Edge.hpp"
#include "Ref_Parser.hpp"
#include "globals.hpp"

#include <cstddef>


namespace cuttlefish
{

template <uint16_t k>
Discontinuity_Graph_Bootstrap<k>::Discontinuity_Graph_Bootstrap(const std::string& cdbg_path, const uint16_t l, const std::size_t part_count, const std::string& graph_path, const uint64_t seed):
      cdbg_path(cdbg_path)
    , l(l)
    , minimizer_seed(seed)
    , part_count(part_count)
    , graph_path(graph_path)
    , phi(phi_label)
    , E(part_count, graph_path)
{}


template <uint16_t k>
void Discontinuity_Graph_Bootstrap<k>::generate()
{
    uint64_t vertex_count = 0;  // Number of discontinuity k-mers.
    uint64_t edge_count = 0;    // Number of edges in the discontinuity graph.

    uint64_t trivial_unitig_count = 0;  // Number of maximal unitigs in the dBG without a discontinuity k-mer.
    std::size_t max_trivial_len = 0;    // Length of the longest maximal unitig without a discontinuity k-mer.

    Ref_Parser parser(cdbg_path);
    while(parser.read_next_seq())
    {
        const auto seq = parser.seq();
        Minimizer_Iterator<decltype(seq), true> min_it(seq, parser.seq_len(), k - 1, l, minimizer_seed);
        minimizer_t last_min, last_min_idx, min, min_idx;

        Kmer<k> first_vertex;   // First discontinuity k-mer in the sequence.
        Kmer<k> last_vertex;    // Last discontinuity k-mer preceding the current k-mer.
        Kmer<k> p = phi;

        std::size_t kmer_idx = 0;   // Index of the current k-mer being checked.
        Kmer<k> kmer(seq, 0), kmer_bar = kmer.reverse_complement(); // Current k-mer.
        side_t s_x, s_y;    // Sides of vertices to which a discontinuity-edge is incident to.
        bool on_chain = false;  // Whether at least one discontinuity vertex on the chain of the current sequence has been visited.

        // Iterate over minimizers of `(k - 1)`-mers.
        min_it.value_at(last_min, last_min_idx);
        while(++min_it)
        {
            if(kmer_idx > 0)
                kmer.roll_to_next_kmer(seq[kmer_idx + k - 1], kmer_bar);

            min_it.value_at(min, min_idx);
            if(min_idx != last_min_idx)
            {
                vertex_count++;

                if(on_chain)
                {
                    const Kmer<k>& x = last_vertex;
                    const Kmer<k>& y = kmer;
                    Kmer<k> x_hat = x.canonical();
                    Kmer<k> y_hat = y.canonical();
                    s_x = (x == x_hat ? side_t::back : side_t::front);
                    s_y = (y == y_hat ? side_t::front : side_t::back);

                    edge_count++;
                    // TODO: place locally-maximal unitig to some bucket `U_j`.
                    E.add(x_hat, s_x, y_hat, s_y, 1, 0);
                }
                else
                    first_vertex = kmer,
                    on_chain = true;

                last_vertex = kmer;
            }

            last_min = min, last_min_idx = min_idx;
            kmer_idx++;
        }


        if(!on_chain)
            trivial_unitig_count++,
            max_trivial_len = std::max(max_trivial_len, parser.seq_len());
        else
        {
            edge_count += 2;
            // TODO: place locally-maximal unitig to some bucket `U_j`.
            E.add(  p, side_t::back,
                    first_vertex.canonical(), first_vertex == first_vertex.canonical() ? side_t::front : side_t::back,
                    1, 0, true);
            E.add(  p, side_t::back,
                    last_vertex.canonical(), last_vertex == last_vertex.canonical() ? side_t::back : side_t::front,
                    1, 0, true);
        }
    }

    // E.serialize();


    std::cerr << "#vertices: " << vertex_count << "\n";
    std::cerr << "#edges: " << edge_count << "\n";
    std::cerr << "#trivial-unipaths: " << trivial_unitig_count << "\n";
    std::cerr << "max-trivial-len: " << max_trivial_len << "\n";
}

}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::Discontinuity_Graph_Bootstrap)
