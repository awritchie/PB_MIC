#include "read_pqr.hpp"

void read_pqr(t_options &options, std::vector<t_atoms> &atoms)
{
    std::string line;
    std::ifstream file(options.pqr.c_str());
    // Format is atomline[0] = x, atomline[1] = y, atomline[2] = z,
    //           atomline[3] = q, atomline[4] = r,
    // with the possibility of adding dipole and quadrupole moments
    t_atoms atom;
    std::string junk;
    if (file.is_open())
    {
        while (file.good())
        {
            getline(file,line);
            if ((not line.empty()) && (line.substr(0,1) != ";") && (line.substr(0,1) != "#"))
            {
                std::stringstream linestream(line);
                linestream >> junk >> junk >> junk >> junk >> junk >> atom.xyz[0] >> atom.xyz[1] >> atom.xyz[2] >> atom.q[0] >> atom.r;
                atoms.push_back(atom);
            }
        }
    }
    else if (!file)
    {
        std::cerr << "\nError reading " << options.pqr << "." << std::endl;
        std::exit(1);
    }
    file.close();

    return;
}
