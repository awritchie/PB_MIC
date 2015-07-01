#include "includes.hpp"
#include "read_pqr.hpp"
#include "make_surface_grid.hpp"
#include "grvy_input_parser.hpp"
#include "electrostatics.hpp"



int main(int argc, const char * argv[])
{
    // Parse through the input file
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input.in>.\n";
        std::exit(1);
    }
    t_options options;
    grvy_input_parser(argc, argv, options);
    readback_options(options);

    // Read the structure file
    std::vector<t_atoms> atoms;
    read_pqr(options, atoms);
    
    // Discretize to a grid
    make_grid(options, atoms);
     
    return 0;
}

