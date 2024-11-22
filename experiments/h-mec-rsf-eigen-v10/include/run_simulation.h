#ifndef RUN_SIMULATION_H
#define RUN_SIMULATION_H

#include "constants.hpp"
#include "global_variables.hpp"
#include "functions.hpp"
#include "hdf5.hpp"
#include <iostream>
#include <H5Cpp.h>
#include "functions.hpp"
#include <fstream>
#include <iomanip>
#include "write_output.h"
#include "LR_setup.h"
#include "move_markers.h"
#include "process_p_cells.h"
#include "power_law.h"


#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/PardisoSupport>

using namespace std;

void run_simulation(int &timestep);

#endif