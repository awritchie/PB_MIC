#ifndef PB_MIC_grvy_input_parser_hpp
#define PB_MIC_grvy_input_parser_hpp

#include <iostream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <grvy.h>
#include <cstdlib>
#include "includes.hpp"
#include "electrostatics.hpp"

void grvy_input_parser(int argc, const char *argv[], t_options &options);
void readback_options(t_options options);

#endif


