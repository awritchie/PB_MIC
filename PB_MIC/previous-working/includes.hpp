#ifndef PB_MIC_includes_hpp
#define PB_MIC_includes_hpp

#include <iostream>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <cstdlib>
#include <stdio.h>

#ifdef __ICC
#include "mkl.h"
#else
#include "cpp_blas.hpp"
#include "cpp_lapack.hpp"
#endif

enum Model
{
    DESMO,
    GCOSMO,
    IEF
};

struct t_atoms
{
    float xyz[3];
    float q[10];
    float r;
    std::vector<int> pg;
};

struct t_grid
{
    float xyz[3];
    float N[3];
    float Phi0, PhiK, r_K, a;
    bool keep;
};

struct t_options
{
    int atomsize;

    std::string pqr, methodname;
    Model method;
    int nlebedev;
    float sdie, temp;
    float pion_conc, pion_q, pion_r;
    float nion_conc, nion_q, nion_r;
    float I, R;
    float deblen, kappa;

    //debug
    int print_surface_grid;
    int print_atom_grid;
    int print_method_name;
    //system
    int nthreads;
};


#endif


