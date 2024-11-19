# H-MECs-Ginkgo

2D Hydro-Mechanical Earthquake Cycle

Computational Earthquake Physics ETH Zurich, 2022 Dal Zilio, L., Hegyi, B., Behr, W. M., Gerya, T. (2022) Hydro-mechanical earthquake cycles in a poro-visco-elasto-plastic fluid-bearing fault structure DOI: https://doi.org/10.1016/j.tecto.2022.229516

## Structure

This repository aims to create a first prototype for GPU-targeted H-MEC code, utilizing the Ginkgo library.

```bash
H-MECs-Ginkgo
├── cmake                 # .cmake files for build system
├── CMakeLists.txt
├── README.md
└── experiments           # various numerical experiments
    ├── demo-XX            
    └── h-mec-rsf-vX      # h-mec code with various versions `X`
```

## Getting Started

To get started, we assume you have already built and installed the Ginkgo library as a prerequisite. If you have not yet done it, please refer to the tutorial for [installing Ginkgo](https://github.com/ginkgo-project/ginkgo/wiki/Tutorial-1:-Getting-Started). In case you encounter any problems with building Ginkgo, you can refer to the [discussions section](https://github.com/ginkgo-project/ginkgo/discussions) for help.

Once you have Ginkgo installed in your system, now you can open a UNIX terminal and go to a dedicated directory where you want this tutorial series to locate at. Then issue the following command

```bash
git clone git@github.com:youwuyou/H-MECs-Ginkgo.git
```

This will create a subdirectory H-MECs-Ginkgo containing codes and other data needed for the tutorial examples.

Now we create a build directory and initialize the build system. By default, we compile the repository with tests, you could also switch it off by using the flag `-DBUILD_TESTS=OFF` if you just want to build the examples themselves without testing.

```bash
cd H-MECs-Ginkgo
mkdir build
cd build
cmake .. && make
```

Now let us run a demonstration code. For running this code, we assume you already have the [OpenCV library](https://opencv.org/) installed. And make sure you locate at `H-MECs-Ginkgo/build` directory.

```bash
cd experiments/demo-heat-equation
./demo-heat-equation
```

If the program successfully runs, it will creates a video file using OpenCV and a custom color mapping as described [here](https://ginkgo-project.github.io/ginkgo-generated-documentation/doc/develop/heat_equation.html).

## Other Dependencies

On Ubuntu OS, for examples involving Eigen and HDF5 libraries, you can install them using:

```bash
sudo apt update # update the package list
sudo apt install libeigen3-dev # install Eigen library
sudo apt-get install libhdf5-serial-dev # install HDF5
```

## Testing

If you have compiled the repository with `-DBUILD_TESTS=ON`, you can run the tests using `make test`.


# H-MECs
This repository contains the code used for Hydro-Mechanical Earthquake Cycles (H-MECs) simulations. It relies on the Ginkgo library, which makes it possible to run the code on CPU and GPU with minimal changes, Nvidia, AMD and Intel are supported GPUs.
The simulation is based on the follwing paper:
Computational Earthquake Physics ETH Zurich, 2022 Dal Zilio, L., Hegyi, B., Behr, W. M., Gerya, T. (2022) Hydro-mechanical earthquake cycles in a poro-visco-elasto-plastic fluid-bearing fault structure DOI: https://doi.org/10.1016/j.tecto.2022.229516

## Structure
The repository is structured as follows:
```bash
H-MECs
├── cmake                   # cmake helper script to help installing dependencies
├── CMakeLists.txt
├── README.md
└── experiments             # Future numerical experiments can be added here
    └── h-mec-rsf-ginkgo-v1
```

## Getting Started
The code depends on two libraries, Ginkgo for linear algebra and HDF5 as efficient file system. Before building the code repository make sure to have the dependencies installed:
### Installing Ginkgo
The code is built on the numerical algebra library Ginkgo. To install the Ginkgo library, refer to the install instructions from [Ginkgo](https://github.com/ginkgo-project/ginkgo/wiki/Tutorial-1:-Getting-Started).
### Installing HDF5
Instructions to install HDF5 can be found on the [HDF5](https://www.hdfgroup.org/download-hdf5/) website.

### Install H-MECs
Once you have the Ginkgo library installed you can proceed. In a terminal go to your preferred target directory and enter the following command to clone the git repository:
```bash
git clone git@github.com:NikMeier/H-MECs.git
```
The simulation is then built with standard cmake procedure:

```bash
cd H-MECs
mkdir build
cd build
cmake .. && cmake --build .
```

In the build directory you will find the experiment folder, containing all the experiments, for now only 'h-mecs-rsf-ginkgo-v1'.

## Running a experiment
To run the code go to the directory and run the executable:
```bash
cd experiments/h-mec-rsf-ginkgo-v1
./h_mec_rsf_ginkgo_v1
```
The output should look like this:
```bash
A 401 x 51 grid is simulated, starting at timestep = 0!
>> VW width = 33 (km)
>> Critical nucleation size = 2.47255 (km)
>> Cohesive zone = 500.691 [m]
Iteration #1!
dt = 5e+08s
Iteration #2!
dt = 5.87652e+06s
...
```

## Running the experiment on GPU
To run the code on a GPU make sure to have the following option turned on, for AMD/Intel GPUs turn the corresponding options to ON:
``` py title="H-MECs/experiments/h-mec-rsf-ginkgo-v1/CMakeLists.txt"
-DGINKGO_BUILD_CUDA=ON          # Build using CUDA for Nvidia devices
```

Finally
## Next steps
