#ifndef PB_MIC_electrostatics_hpp
#define PB_MIC_electrostatics_hpp

#include <iostream>
#include <omp.h>
#include "includes.hpp"

#ifndef M_PI
#define M_PI = atan(1.)*4.
#endif

#define e_c  4.803242384e-10
#define N_A  6.022045000e+23
#define k_B  1.380662000e-16

float getKappa(t_options options);

float getPhi(float &r, float *q);

float getPhi_Kappa(float &r, float *q, float &kappa, float &ionR, float &atomR, float &sdie);

float calc_r( int size, float *A, float *B );

void desmo(t_options &options, std::vector<t_grid> &surface, float *K, float *R);

void gcosmo(t_options &options, std::vector<t_grid> &surface, float *K, float *R);

void ief(t_options &options, std::vector<t_grid> &surface, float *K, float *R);





#endif


