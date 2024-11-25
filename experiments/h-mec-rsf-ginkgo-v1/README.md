# H-MECs-Ginkgo-V1
## Organization
The simulation is organized in the following way:
```bash
h-mecs-rsf-ginkgo-v1
├── CMakeLists.txt
├── h_mec_rsf_ginkgo_v1.cpp     # main function
└── include
    ├── run_simulation_gko.cpp  # computation is in here
    ├── constants
    ├── global_variables.hpp    # Declare all global variables
    ├── ginkgo_exec.hpp         # Set up ginkgo exec. To switch from CPU to GPU make changes in here
    ├── hdf5                    # Handle hdf5 reading/writing
    ├── init_geometry           # Set up marker geometry
    ├── power_law               # Update parameters after timestep
    ├── process_p_cells         # Process pressure cells after timestep
    ├── read_in_matrices        # Read in previous simulation data to continue simulation
    └── write_output
```
## Switching from CPU to GPU
To build the code utilizing a GPU device, make sure you have made these two changes:
```cpp
// In the CMakeLists.txt file, switch the compile option for the respective GPU to ON:
target_compile_options(${PROGRAM_NAME}
PUBLIC  -DGINKGO_BUILD_CUDA=ON 
)
```
and in the file include/ginkgo_exec.cpp
```cpp
// Comment the Reference Executor for debugging and enable:

// OmpExecutor for CPUs
inline const auto exec = gko::OmpExecutor::create();

// CudaExecutor for Nvidia GPUs
inline const auto gpu_exec = gko::CudaExecutor::create(0, gko::OmpExecutor::create()); // If using AMD instead enable the HipExecutor
```

## The simulation timestep loop
Most of the actual computation is done inside the run_simulation_gko.cpp file. Inside the timestep loop we repeat the following steps:
 - Updating variables on the markers
 - Assembling the LSE.
 - Building Solver & solving the LSE.
 - Update all variables
 - Write to output files

## About ginkgo
In ginkgo the matrices (and vectors) are objects of type gko::matrix::Dense<double>. The entries of the matrices are accessed by pointers into the gko matrices. When constructing a matrix, a executor has to be chosen, i.e. the matrix lives either on the CPU or the GPU. Assembly of the LSE is done on the CPU, whereas the solving of the LSE is done on the GPU. This means after assembling the LSE we have to copy it to the GPU and afterwards we have to copy the solution from the GPU to the CPU to process it. 

The [Ginkgo](https://ginkgo-project.github.io)-website contains the [Documentation](https://ginkgo-project.github.io/ginkgo-generated-documentation/doc/develop/) and a [Tutorial](https://github.com/ginkgo-project/ginkgo/wiki/Tutorial-1:-Getting-Started) for getting started.

A list of solvers can be found under [Ginkgo Solvers](https://ginkgo-project.github.io/ginkgo-generated-documentation/doc/develop/namespacegko_1_1solver.html). The current direct LU solver is experimental and can be found [here](https://ginkgo-project.github.io/ginkgo-generated-documentation/doc/develop/classgko_1_1experimental_1_1solver_1_1Direct_1_1Factory.html)
