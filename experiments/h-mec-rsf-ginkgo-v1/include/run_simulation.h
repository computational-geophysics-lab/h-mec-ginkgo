#ifndef RUN_SIMULATION_H
#define RUN_SIMULATION_H
#include <vector>
#include "constants.hpp"
#include "global_variables.hpp"
#include "functions.hpp"
#include "hdf5.hpp"
#include <iostream>
#include <H5Cpp.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "write_output.h"
#include "LR_setup.h"
#include "move_markers.h"
#include "process_p_cells.h"
#include "power_law.h"
#include "ginkgo_exec.hpp"


#include <ginkgo/ginkgo.hpp>

using namespace std;

void run_simulation(int &timestep);

#endif