#ifndef PB_MIC_make_surface_grid_hpp
#define PB_MIC_make_surface_grid_hpp

#include <iostream>
#include "sphere_lebedev_rule.hpp"
#include "includes.hpp"
#include <omp.h>
#include "electrostatics.hpp"

int make_grid( t_options &options, std::vector<t_atoms> &atoms );

#endif


