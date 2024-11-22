/* ========================================
H - MEC: Hydro - Mechanical Earthquake Cycle
Computational Earthquake Physics
ETH Zurich, 2022

Dal Zilio, L., Hegyi, B., Behr, W. M., & Gerya, T. (2022)
Hydro - mechanical earthquake cycles in a
poro - visco - elasto - plastic fluid - bearing fault.
DOI: http://arxiv.org / abs / 2201.11786

========================================
Solving of compressible Stokes + continuity equations
with variable viscosity and Darsi + continuity equations
Total X - Stokes: dSIGMAxxt' / dx + dSIGMAxyt' / dy - dPt / dx = -RHOt * gx
Total Y - Stokes: dSIGMAyxt' / dx + dSIGMAyyt' / dy - dPt / dy = -RHOt * gy
Solid Continuity: dVxs / dx + dVys / dy + (Pt - Pf) / ETAbulk = 0
Fluid X - Darsi: - ETAfluid / K * VxD - dPf / dx = -RHOf * gx
Fluid Y - Darsi: - ETAfluid / K * VyD - dPf / dy = -RHOf * gy
Fluid Continuity: dVxD / dx + dVyD / dy - (Pt - Pf) / ETAbulk = 0
+ staggered grid
+ P - v formulation
+ advecting material fileds by markers
========================================
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/PardisoSupport>
#include <random>
#include <ctime>
#include <chrono>
#include <math.h>
#include <iomanip>
#include <H5Cpp.h>

#include "include/hdf5.hpp"
#include "include/constants.hpp"
#include "include/global_variables.hpp"
#include "include/functions.hpp"
#include "include/init_geometry.h"
#include "include/read_in_matrices.h"
#include "include/run_simulation.h"
#include "include/write_output.h"


using namespace std;
using namespace H5;



// Basic nodes
Eigen::MatrixXd OM0(Ny, Nx);  // Old state parameter
Eigen::MatrixXd OM(Ny, Nx);   // State parameter
Eigen::MatrixXd OM5(Ny, Nx);
Eigen::MatrixXd ARSF(Ny, Nx); // a - parameter of RSF
Eigen::MatrixXd BRSF(Ny, Nx); // b - parameter of RSF
Eigen::MatrixXd LRSF(Ny, Nx); // L - parameter of RSF

// Unknown parameters
Eigen::MatrixXd pt(Ny1, Nx1);  // Total pressure
Eigen::MatrixXd vxs(Ny1, Nx1); // Solid vx - velocity
Eigen::MatrixXd vys(Ny1, Nx1); // Solid vy - velocity
Eigen::MatrixXd vzs(Ny1, Nx1); // Solid vz - velocity
Eigen::MatrixXd pf(Ny1, Nx1);  // Fluid pressure
Eigen::MatrixXd vxD(Ny1, Nx1); // Darsi vx - velocity
Eigen::MatrixXd vyD(Ny1, Nx1); // Darsi vy - velocity

// Nodal matrices
// Basic nodes
Eigen::MatrixXd RHO(Ny, Nx), ETA(Ny, Nx), ETA0(Ny, Nx), ETA1(Ny, Nx), ETA5(Ny, Nx), ETA00(Ny, Nx), IETAPLB(Ny, Nx), SXY(Ny, Nx), SXY0(Ny, Nx), SZX0(Ny, Nx), SZY0(Ny, Nx), YNY0(Ny, Nx), KKK(Ny, Nx),
      GGG(Ny, Nx), COHC(Ny, Nx), COHT(Ny, Nx), FRIC(Ny, Nx), FRIT(Ny, Nx), DILC(Ny, Nx), TTT(Ny, Nx), EIIB(Ny, Nx), VSLIPB(Ny, Nx), SZX(Ny, Nx), SZY(Ny, Nx);

// Pressure nodes
Eigen::MatrixXd ETAB(Ny1, Nx1), ETAB0(Ny1, Nx1), ETAP(Ny1, Nx1), ETAP0(Ny1, Nx1), POR(Ny1, Nx1), GGGP(Ny1, Nx1), GGGB(Ny1, Nx1), PTF0(Ny1, Nx1), PT0(Ny1, Nx1), PF0(Ny1, Nx1), pt_ave(Ny1, Nx1),
      pf_ave(Ny1, Nx1), SXX(Ny1, Nx1), SXX0(Ny1, Nx1), SYY(Ny1, Nx1), SYY0(Ny1, Nx1), DILP(Ny1, Nx1);
// Vx nodes
Eigen::MatrixXd RHOX(Ny1, Nx1), RHOFX(Ny1, Nx1), ETADX(Ny1, Nx1), PORX(Ny1, Nx1), VX0(Ny1, Nx1), VXF0(Ny1, Nx1);
// Vy nodes
Eigen::MatrixXd RHOY(Ny1, Nx1), RHOFY(Ny1, Nx1), ETADY(Ny1, Nx1), PORY(Ny1, Nx1), VY0(Ny1, Nx1), VYF0(Ny1, Nx1);
// Vz nodes
Eigen::MatrixXd VZ0(Ny1, Nx1);

Eigen::MatrixXd ESP(Ny, Nx), EXY(Ny, Nx), EZX(Ny, Nx), EZY(Ny, Nx), EXX(Ny1, Nx1), EYY(Ny1, Nx1), EII(Ny1, Nx1), EIIVP(Ny1, Nx1), SII(Ny1, Nx1), DSII(Ny1, Nx1), DIS(Ny1, Nx1);

Eigen::MatrixXd EL_DECOM(Ny1, Nx1);   // Elastic (de)compaction
 Eigen::MatrixXd VIS_COMP(Ny1, Nx1);

// Lagrangian solid markers
Eigen::VectorXd t_marker = Eigen::VectorXd::Constant(marknum, 1);         // Marker rock type
Eigen::VectorXd rhom = Eigen::VectorXd::Constant(marknum, 2800);          // Density of solid
Eigen::VectorXd etasm = Eigen::VectorXd::Constant(marknum, 1e21);         // Standard shear viscosity of bulk
Eigen::VectorXd etam(marknum);                                  // Shear viscosity of bulk
Eigen::VectorXd cohescm = Eigen::VectorXd::Constant(marknum, cohes);      // Cohesion for confined fracture of solid
Eigen::VectorXd cohestm = Eigen::VectorXd::Constant(marknum, cohes);      // Cohesion for tensile fracture of solid
Eigen::VectorXd frictcm = Eigen::VectorXd::Constant(marknum, .5);         // friction for confined fracture of solid
Eigen::VectorXd dilatcm = Eigen::VectorXd::Constant(marknum, dilatation); // dilatation for confined fracture of solid
Eigen::VectorXd fricttm = Eigen::VectorXd::Constant(marknum, tensile);    // friction for tensile fracture of solid
Eigen::VectorXd porm(marknum);                                  // Porosity of solid
Eigen::VectorXd kkkm = Eigen::VectorXd::Constant(marknum, 2e-16);         // Standard permeability of solid
Eigen::VectorXd xm(marknum);                                    // Horizontal coordinates of solid markers
Eigen::VectorXd ym(marknum);                                    // Vertical coordinates of solid markers
Eigen::VectorXd sxxm(marknum);                                  // Marker SIGMAxx', Pa
Eigen::VectorXd syym(marknum);                                  // Marker SIGMAyy', Pa
Eigen::VectorXd sxym(marknum);                                  // Marker SIGMAxy', Pa

// (3) Defining global matrixes
// according to the global number of unknowns
// Sparse Matrix L is not yet defined as it will be built from a set of Triplets each step
Eigen::VectorXd R(N); // Vector of the right parts of equations
Eigen::VectorXd S(N);

// variable type declaration
int ynlast, iterstep;

double timesum = 0;
double dt = dtelastic0;
bool yndtdecrease = true;

double dt00, dtx, dty, dtz, dtlapusta, Vmax, maxvxy;

Eigen::VectorXd DSYLSQ(niterglobal);

Eigen::MatrixXd DVX0(Ny1, Nx1), DVY0(Ny1, Nx1), DVZ0(Ny1, Nx1), DSY(Ny, Nx), YNY(Ny, Nx), SigmaY(Ny, Nx), SII_fault(Ny, Nx), SIIB(Ny, Nx); // DVX0, DVY0, DVZ0 unused

Eigen::VectorXd timesumcur(num_timesteps), dtcur(num_timesteps);
Eigen::VectorXd maxvxsmod(num_timesteps), minvxsmod(num_timesteps), maxvysmod(num_timesteps), minvysmod(num_timesteps), maxvzsmod(num_timesteps), minvzsmod(num_timesteps);

// variables for timer 
std::chrono::time_point<std::chrono::system_clock> time_start, time_end; // total runtime

int main() {
    time_start = std::chrono::system_clock::now();
    Eigen::initParallel();
    srand(time(nullptr));

    // read input
    // Load file
    string timestep_str;
    ifstream input_timestep("./input_data/StartingTimestep.txt");
    getline(input_timestep, timestep_str);
    int timestep = stoi(timestep_str);
    input_timestep.close();

    init_geometry();

    read_in_matrices(timestep);
    
    // /////////////////////////////////////////////////////////////////////////////////////// 
    // actual computations start here
    // /////////////////////////////////////////////////////////////////////////////////////// 

    run_simulation(timestep);

    time_end = std::chrono::system_clock::now();

    // total run time
    std::chrono::duration<double> time_tot = time_end - time_start;
    cout << "====================================" << endl;
    cout << "total time: " << time_tot.count() << " sec" << endl;
    cout << "====================================" << endl;
    
    return 0;
}
