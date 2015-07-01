#include "grvy_input_parser.hpp"

void grvy_input_parser(int argc, const char *argv[], t_options &options)
{
    int igot;
    char *pqr,*method;
    // Initialize/read the file
    igot = grvy_input_fopen(argv[1]);
    if (igot == 0)
    {
        std::exit(1);
    }

    // Read specific input
    if (grvy_input_fread_char("pqr",&pqr))
    {
        std::string pqrname(pqr);
        options.pqr = pqrname;
        free(pqr);
    }
    if (! grvy_input_fread_int("nlebedev",&options.nlebedev))
        options.nlebedev = 50;
    
    // Electrostatic settings
    if (grvy_input_fread_char("elec/method",&method))
    {
        std::string methodname(method);
        options.methodname = methodname;
        if (std::strncmp(method, "DESMO", sizeof(method) - 1) == 0)
        {
            options.method = DESMO;
        }
        else if (std::strncmp(method, "GCOSMO", sizeof(method) - 1) == 0)
        {
            options.method = GCOSMO;
        }
        else if (std::strncmp(method, "IEF-PCM", sizeof(method) - 1) == 0)
        {
            options.method = IEF;
        }
        else
        {
            fprintf(stderr,"\nError: PCM %s not recognized\n",method);
            std::exit(1);
        }
        free(method);
    }
    else
    {
        options.methodname = "DESMO";
        options.method = DESMO;
    }
    if (! grvy_input_fread_float("elec/sdie",&options.sdie))
        options.sdie = 80;
    if (! grvy_input_fread_float("elec/temp",&options.temp))
        options.temp = 298.15;

    // Ions
    if (! grvy_input_fread_float("ions/pion_conc",&options.pion_conc))
        options.pion_conc = 0.0f;
    if (! grvy_input_fread_float("ions/pion_charge",&options.pion_q))
        options.pion_q = 1.0;
    if (! grvy_input_fread_float("ions/pion_radius",&options.pion_r))
        options.pion_r = 2.0;
    if (! grvy_input_fread_float("ions/nion_conc",&options.nion_conc))
        options.nion_conc = 0.0f;
    if (! grvy_input_fread_float("ions/nion_charge",&options.nion_q))
        options.nion_q = -1.0;
    if (! grvy_input_fread_float("ions/nion_radius",&options.nion_r))
        options.nion_r = 2.0;
    if (options.nion_conc != options.pion_conc)
    {
        fprintf(stderr,"\nError, you have a counterion charge imbalance!  Net charge conc. = %.3f\n",options.nion_conc + options.pion_conc);
        std::exit(1);
    }
    // Debug
    if (! grvy_input_fread_int("debug/print_surface_grid",&options.print_surface_grid))
        options.print_surface_grid = 0;
    if (! grvy_input_fread_int("debug/print_atom_grid",&options.print_atom_grid))
        options.print_atom_grid = 0;
    if (! grvy_input_fread_int("debug/print_method_name",&options.print_method_name))
        options.print_method_name = 0;
    
    // system
    if (! grvy_input_fread_int("system/nthreads",&options.nthreads))
        options.nthreads = 0;

     /* If I change this value, I need to hand-edit make_surface_grid.{hpp,cpp} */
     options.atomsize=5; // May change this to handle AMOEBA multipoles eventually

    // Calculate ionic strength = 1/2 sum([mol/L]_i * q_i**2)
     options.I = 0.5*(options.pion_conc*options.pion_q*options.pion_q + options.nion_conc*options.nion_q*options.nion_q);
     if (options.nion_r > options.pion_r)
     {
         options.R = options.nion_r;
     }
     else
     {
         options.R = options.pion_r;
     }

     options.kappa = getKappa(options);
     options.deblen = 1. / options.kappa;

    return;
}

void readback_options(t_options options)
{
    fprintf(stderr,"%s\n",std::string(80,'=').c_str());
    fprintf(stderr,"%-22s: %56s\n","PQR file",options.pqr.c_str());
    fprintf(stderr,"%-22s: %56d\n","Number of threads",options.nthreads);
    fprintf(stderr,"\n");
//    fprintf(stderr,"%s\n",std::string(80,'-').c_str());
    fprintf(stderr,"%-22s: %56s\n","PCM Used",options.methodname.c_str());
    fprintf(stderr,"%-22s: %56d\n","N Lebedev points",options.nlebedev);
    fprintf(stderr,"%-22s: %56.3f\n","Solvent Dielectric",options.sdie);
    fprintf(stderr,"%-22s: %54.3f K\n","Temperature",options.temp);
    fprintf(stderr,"%-22s: %54.6f A\n","Debye Length",options.deblen);
    fprintf(stderr,"\n");
//    fprintf(stderr,"%s\n",std::string(80,'-').c_str());
    fprintf(stderr,"%-22s: %39.3f M ionic strength\n","Ionic Species",options.I);
    fprintf(stderr,"%-22s: %47.3f A-radius\n","Max Ion Radius",options.R);
    fprintf(stderr,"%-22s: %6.3f A-radius, %6.3f e-charge, %6.3f M concentration\n"," ",options.pion_r, options.pion_q, options.pion_conc);
    fprintf(stderr,"%-22s: %6.3f A-radius, %6.3f e-charge, %6.3f M concentration\n"," ",options.nion_r, options.nion_q, options.nion_conc);
    fprintf(stderr,"\n");
//    fprintf(stderr,"%s\n",std::string(80,'-').c_str());
    fprintf(stderr,"%-22s: %56d\n","Print atom grid",options.print_atom_grid);
    fprintf(stderr,"%-22s: %56d\n","Print surface grid",options.print_surface_grid);
    fprintf(stderr,"%-22s: %56d\n","Print method name",options.print_method_name);
    
    fprintf(stderr,"%s\n",std::string(80,'=').c_str());
    return;
}
