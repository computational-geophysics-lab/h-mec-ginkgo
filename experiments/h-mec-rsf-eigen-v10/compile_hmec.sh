#!/bin/sh

MAIN_FILE="h_mec_rsf_v10.cpp"
LINK_FILE1="include/hdf5.cpp"
LINK_FILE2="include/functions.cpp"
LINK_FILE3="include/init_geometry.cpp"
LINK_FILE4="include/read_in_matrices.cpp"
LINK_FILE5="include/run_simulation.cpp"
LINK_FILE6="include/write_output.cpp"
LINK_FILE7="include/LR_setup.cpp"
LINK_FILE8="include/move_markers.cpp"
LINK_FILE9="include/process_p_cells.cpp"
LINK_FILE10="include/power_law.cpp"
declare -i NUM_THREADS=4
declare -i START_TIMESTEP=0

echo
echo "=============================================================================================="
echo "Submitting ${MAIN_FILE} using $NUM_THREADS cores and starting from timestep $START_TIMESTEP"
echo
echo "Terminate running job..."
echo
scancel --name=compile_hmec

echo
echo "Removing old files..."
rm -rf EVO*
rm -rf lsf*
rm -rf *txt
if [ $START_TIMESTEP -eq 0 ]
then
	echo "Starting new simulation -> removing old .h5 files..."
	rm -rf *.h5
fi

printf "$START_TIMESTEP" >> file.txt

echo
echo "Set max threads to 1 for compilation"
export OMP_NUM_THREADS=1


#SBATCH -n 2
#SBATCH --time=0:30:00
#SBATCH --mem-per-cpu=50000
#SBATCH --open-mode=truncate
#SBATCH --job-name=compile_hmec
#SBATCH --output=cpp_pardiso.out
#SBATCH --error=cpp_pardiso.err

echo
echo "Loading GCC Compiler, openmpi, Eigen Library, HDF5 library, libszip and petcs library  module..."
module load gcc/9.3.0 openmpi/4.1.4 eigen/3.3.9 hdf5/1.10.1 libszip/2.1.1 petsc/3.15.5

export PATH=/cluster/apps/nss/intel/oneapi/2022.1.2/mkl/2022.0.2/bin/intel64:$PATH
export LIBRARY_PATH=/cluster/apps/nss/intel/oneapi/2022.1.2/mkl/2022.0.2/lib/intel64:$LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/apps/nss/intel/oneapi/2022.1.2/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
export CPATH=/cluster/apps/nss/intel/oneapi/2022.1.2/mkl/2022.0.2/include:$CPATH
export INCLUDE=/cluster/apps/nss/intel/oneapi/2022.1.2/mkl/2022.0.2/include:$INCLUDE
export CPLUS_INCLUDE_PATH=/cluster/apps/nss/intel/oneapi/2022.1.2/mkl/2022.0.2/include:$CPLUS_INCLUDE_PATH
export MKLROOT=/cluster/apps/nss/intel/oneapi/2022.1.2/mkl/2022.0.2

echo "Compiling Program..."
echo

mpic++ -std=c++17 -o model_cpp -O3 -DEIGEN_USE_MKL_ALL -L${EIGEN_ROOT}/lib64 -L${HDF5_ROOT}/lib -lhdf5 -lhdf5_cpp -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -ftree-vectorize -march=native -DMKL_LP64 $MAIN_FILE $LINK_FILE1 $LINK_FILE2 $LINK_FILE3 $LINK_FILE4 $LINK_FILE5 $LINK_FILE6 $LINK_FILE7 $LINK_FILE8 $LINK_FILE9 $LINK_FILE10 -L${INTEL_ROOT}/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -L/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_hl_cpp.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_cpp.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_hl.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5.a -lhdf5 -lsz -lz -lrt -Wl,-rpath -Wl,/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib

echo "Run file"
sbatch --time=24:00:00 --ntasks=4 --wrap="mpirun model_cpp"
echo "=============================================================================================="
