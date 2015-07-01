#include "make_surface_grid.hpp"

int n_lebedev_points[31] = {    6,   14,   26,   38,   50,   74,   86,  110,
                              146,  170,  194,  230,  266,  302,  350,  434,
                              590,  770,  974, 1202, 1454, 1730, 2030, 2354,
                             2702, 3074, 3890, 4334, 4802, 5294, 5810 };

int make_grid( t_options &options, std::vector<t_atoms> &atoms )
{
    // Make sure the number of lebedev grid points chosen is calculated 
    // in the algorithm
    bool inLebedev = false;
    for (int i=0; i<31; i++)
    {
        if (options.nlebedev == n_lebedev_points[i])
        {
            inLebedev = true;
            break;
        }
    }
    if (! inLebedev )
    {
        std::cerr << "\nError: " << options.nlebedev << " selected for lebedev, but only [ ";
#ifdef __ICC
        std::cerr << n_lebedev_points[:] << " ";
#else
        for (int i=0;i<32;i++)
        {
            std::cout << n_lebedev_points[i] << " ";
        }
#endif
        std::cerr << "] points per atom allowable.\n";
        std::exit(1);
    }

    // Make Lebedev sphere
    double lx[options.nlebedev], ly[options.nlebedev], lz[options.nlebedev], lw[options.nlebedev];
    ld_by_order(options.nlebedev,lx,ly,lz,lw);
    
    // Apply grid points to structure
    int natoms = atoms.size();
    int n = natoms*options.nlebedev;
    
    // Make surface array
    std::vector<t_grid> surface(n);
    
    // Loop over every atom
    int NTHREADS = options.nthreads;
    #pragma omp parallel for num_threads(NTHREADS) reduction(-:n)
    for (int i=0; i<natoms; i++)
    {
        //int tid=omp_get_thread_num();
        //printf("Atom %i on thread-%d\n",i,tid);
        for (int p=0; p<options.nlebedev; p++)
        {
            int x = p+i*options.nlebedev;

            // R_p = lebedev_p * Radius_atom_j + R_atom_j
            surface[x].xyz[0] = (lx[p]*atoms[i].r+atoms[i].xyz[0]);
            surface[x].xyz[1] = (ly[p]*atoms[i].r+atoms[i].xyz[1]);
            surface[x].xyz[2] = (lz[p]*atoms[i].r+atoms[i].xyz[2]);
            // Need to know the radius of the grid point's originating atom
            surface[x].r_K = atoms[i].r;
            // a_i = w_i * R_I * R_I hwere w_i is the Lebedev quadrature weight and R_I is the atom radius
            surface[x].a = lw[p] * atoms[i].r * atoms[i].r;
            surface[x].keep = true;
            for (int j=0; j<natoms; j++)
            {
                float r = calc_r(3, &surface[x].xyz[0], &atoms[j].xyz[0]);
                // Discard points which are closer to an atom than that atom's
                // radius, so long as the point doesn't originate from that 
                // atom.
                if ( i != j && r < atoms[j].r )
                {
                    surface[x].keep = false;
                    n--;
                    break;
                }    
                else
                {
                    surface[x].Phi0 += getPhi(r, &atoms[j].q[0]);
                    surface[x].PhiK += getPhi_Kappa(r, &atoms[j].q[0], options.kappa, options.R, atoms[j].r, options.sdie);
                }
            }
            if (surface[x].keep)
            {
                // Calculate the surface normal vector
                surface[x].N[0] = lx[p];
                surface[x].N[1] = ly[p];
                surface[x].N[2] = lz[p];
            }
        }
    }
    
    // Print out the surface grid for every atom
    if (options.print_atom_grid == 1)
    {
        for (int i=0; i<options.nlebedev*natoms;i++)
        {
            for (int j=0; j<3; j++)
            {
                std::cout << surface[i].xyz[j] << " ";
            }
            std::cout << surface[i].Phi0 << " " << surface[i].PhiK << "\n";
        }
    }

    // Keep only surface points
    std::vector<t_grid> abs_surface(n);
    float *v = new float [n];
    int k = 0;
    for (int i=0; i<natoms*options.nlebedev; i++)
    {
        if (surface[i].keep)
        {
            abs_surface[k] = surface[i];
            v[k] = surface[i].Phi0;
            if (options.print_surface_grid)
            {
                for (int j=0;j<3;j++)
                {
                    std::cout << abs_surface[k].xyz[j] << " ";
                }
            }
            if (options.print_surface_grid == 1)
            {
                std::cout << abs_surface[k].Phi0 << " " << abs_surface[k].PhiK << "\n";
            }
            k++;
        }
    }

    fprintf(stderr,"Keeping %d of %d grid points.\n",n,natoms*options.nlebedev);
    
    std::size_t N = pow(n,2);
    float *K = new float [N];
    float *R = new float [N];
    
    switch(options.method)
    {
        case DESMO: desmo(options, abs_surface, K, R);
            break;
        case GCOSMO: gcosmo(options, abs_surface, K, R);
            break;
        case IEF: ief(options, abs_surface, K, R);
            break;
        default:
            fprintf(stderr,"\nError: PCM %s not recognized\n",options.methodname.c_str());
            std::exit(1);
    }
    
/*    To-Do list:
 * 1) Add sgesv to cpp_blas.{hpp,cpp}
 * 2) Include this routine in desmo, gcosmo, and ief
 *      2b.  Possibly rename those functions to specify using LU factorization to solve
 * 3) Calculate solvation energy
 * 4) Calculate field!
 */
    float *Rv = new float [n];
    omp_set_num_threads(NTHREADS);
#ifdef __ICC
    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,1,n,1,&R[0],n,&v[0],n,0,&Rv[0],n);
#else
    cblas_sgemm('N','N',n,1,n,1,&R[0],n,&v[0],n,0,&Rv[0],n);
#endif
    float *q = new float [n];
    char trans = 'N';
    int dim = n;
    int nrhs = 1;
    int LDA = n;
    int LDB = n;
    int info ;
    int *ipiv = new int [n];
    //sgesv(&dim, &nrhs, &K[0], &LDA, &ipiv[0], &Rv[0], &LDB, &info); 
    std::cout << " = array(( ";
    for (int i=0;i<10;i++)
    {
        std::cout << Rv[i] << ", ";
    }
    std::cout << "))\n";

    delete[] Rv;
    delete[] v;
    delete[] q;
    delete[] ipiv;
    delete[] K;
    delete[] R;
    
    return 0;
}
