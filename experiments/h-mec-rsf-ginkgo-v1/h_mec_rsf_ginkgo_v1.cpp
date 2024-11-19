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


#include <random>
#include <ctime>
#include <chrono>
#include <math.h>
#include <iomanip>
#include <H5Cpp.h>
#include <ginkgo/ginkgo.hpp>

#include "include/constants.hpp"
#include "include/global_variables.hpp"
#include "include/functions.hpp"
#include "include/init_geometry.h"
#include "include/read_in_matrices.h"
#include "include/run_simulation.h"
#include "include/hdf5.hpp"
#include "include/write_output.h"


/* ====================== Ginkgo =============================
Changes:
 - Replaced all Eigen vectors/matrices with gko matrices*/

using namespace H5;

// Basic nodes
std::unique_ptr<gko::matrix::Dense<>> OM0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny,Nx)); // Old state parameter
double* OM0 = OM0_gko->get_values(); // Ny x Nx, Old state parameter
std::unique_ptr<gko::matrix::Dense<>> OM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny,Nx)); // Ny x Nx, State Parameter
double* OM = OM_gko->get_values(); // Ny x Nx, state parameter
std::unique_ptr<gko::matrix::Dense<>> OM5_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny,Nx)); // Ny x Nx
double* OM5 = OM5_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> ARSF_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny,Nx)); // Ny x Nx, a - parameter of RSF
double* ARSF = ARSF_gko->get_values(); // Ny x Nx, a - parameter of RSF
std::unique_ptr<gko::matrix::Dense<>> BRSF_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny,Nx)); // Ny x Nx, b - parameter of RSF
double* BRSF = BRSF_gko->get_values(); // Ny x Nx, b - parameter of RSF
std::unique_ptr<gko::matrix::Dense<>> LRSF_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny,Nx)); // Ny x Nx, L - parameter of RSF
double* LRSF = LRSF_gko->get_values(); // Ny x Nx, L - parameter of RSF

/*
// Basic nodes
Eigen::MatrixXd OM0(Ny, Nx);  // Old state parameter
Eigen::MatrixXd OM(Ny, Nx);   // State parameter
Eigen::MatrixXd OM5(Ny, Nx);
Eigen::MatrixXd ARSF(Ny, Nx); // a - parameter of RSF
Eigen::MatrixXd BRSF(Ny, Nx); // b - parameter of RSF
Eigen::MatrixXd LRSF(Ny, Nx); // L - parameter of RSF
*/

std::unique_ptr<gko::matrix::Dense<>> pt_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Total pressure
double* pt = pt_gko->get_values(); // Ny1 x Nx1, Total pressure
std::unique_ptr<gko::matrix::Dense<>> vxs_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Solid vx - velocity
double* vxs = vxs_gko->get_values(); // Ny1 x Nx1, Solid vx - velocity
std::unique_ptr<gko::matrix::Dense<>> vys_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Solid vy - velocity
double* vys = vys_gko->get_values(); // Ny1 x Nx1, Solid vy - velocity
std::unique_ptr<gko::matrix::Dense<>> vzs_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Solid vz - velocity
double* vzs = vzs_gko->get_values(); // Ny1 x Nx1, Solid vz - velocity
std::unique_ptr<gko::matrix::Dense<>> pf_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Fluid pressure
double* pf = pf_gko->get_values(); // Ny1 x Nx1, Fluid Pressure
std::unique_ptr<gko::matrix::Dense<>> vxD_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Darsi vx - velocity
double* vxD = vxD_gko->get_values(); // Ny1 x Nx1, Darsi vx - velocity
std::unique_ptr<gko::matrix::Dense<>> vyD_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Darsi vy - velocity
double* vyD = vyD_gko->get_values(); // Ny1 x Nx1, Darsi vy - velocity

/*
// Unknown parameters
Eigen::MatrixXd pt(Ny1, Nx1);  // Total pressure
Eigen::MatrixXd vxs(Ny1, Nx1); // Solid vx - velocity
Eigen::MatrixXd vys(Ny1, Nx1); // Solid vy - velocity
Eigen::MatrixXd vzs(Ny1, Nx1); // Solid vz - velocity
Eigen::MatrixXd pf(Ny1, Nx1);  // Fluid pressure
Eigen::MatrixXd vxD(Ny1, Nx1); // Darsi vx - velocity
Eigen::MatrixXd vyD(Ny1, Nx1); // Darsi vy - velocity
*/

std::unique_ptr<gko::matrix::Dense<>> RHO_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx // Ny x Nx
double* RHO = RHO_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> ETA_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* ETA = ETA_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> ETA0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* ETA0 = ETA0_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> ETA1_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* ETA1 = ETA1_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> ETA5_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* ETA5 = ETA5_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> ETA00_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* ETA00 = ETA00_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> IETAPLB_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* IETAPLB = IETAPLB_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SXY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SXY = SXY_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SXY0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SXY0 = SXY0_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SZX0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SZX0 = SZX0_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SZY0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SZY0 = SZY0_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> YNY0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* YNY0 = YNY0_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> KKK_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* KKK = KKK_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> GGG_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* GGG = GGG_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> COHC_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* COHC = COHC_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> COHT_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* COHT = COHT_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> FRIC_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* FRIC = FRIC_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> FRIT_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* FRIT = FRIT_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> DILC_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* DILC = DILC_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> TTT_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* TTT = TTT_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> EIIB_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* EIIB = EIIB_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> VSLIPB_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* VSLIPB = VSLIPB_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SZX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SZX = SZX_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SZY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SZY = SZY_gko->get_values(); // Ny x Nx

/*
// Nodal matrices
// Basic nodes
Eigen::MatrixXd RHO(Ny, Nx), ETA(Ny, Nx), ETA0(Ny, Nx), ETA1(Ny, Nx), ETA5(Ny, Nx), ETA00(Ny, Nx), IETAPLB(Ny, Nx), SXY(Ny, Nx), SXY0(Ny, Nx), SZX0(Ny, Nx), SZY0(Ny, Nx), YNY0(Ny, Nx), KKK(Ny, Nx),
      GGG(Ny, Nx), COHC(Ny, Nx), COHT(Ny, Nx), FRIC(Ny, Nx), FRIT(Ny, Nx), DILC(Ny, Nx), TTT(Ny, Nx), EIIB(Ny, Nx), VSLIPB(Ny, Nx), SZX(Ny, Nx), SZY(Ny, Nx);
*/

std::unique_ptr<gko::matrix::Dense<>> ETAB_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* ETAB = ETAB_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> ETAB0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* ETAB0 = ETAB0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> ETAP_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* ETAP = ETAP_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> ETAP0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* ETAP0 = ETAP0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> POR_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* POR = POR_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> GGGP_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* GGGP = GGGP_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> GGGB_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* GGGB = GGGB_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> PTF0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* PTF0 = PTF0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> PT0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* PT0 = PT0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> PF0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* PF0 = PF0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> pt_ave_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* pt_ave = pt_ave_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> pf_ave_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* pf_ave = pf_ave_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> SXX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* SXX = SXX_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> SXX0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* SXX0 = SXX0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> SYY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* SYY = SYY_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> SYY0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* SYY0 = SYY0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> DILP_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* DILP = DILP_gko->get_values(); // Ny1 x Nx1

/*
// Pressure nodes
Eigen::MatrixXd ETAB(Ny1, Nx1), ETAB0(Ny1, Nx1), ETAP(Ny1, Nx1), ETAP0(Ny1, Nx1), POR(Ny1, Nx1), GGGP(Ny1, Nx1), GGGB(Ny1, Nx1), PTF0(Ny1, Nx1), PT0(Ny1, Nx1), PF0(Ny1, Nx1), pt_ave(Ny1, Nx1),
      pf_ave(Ny1, Nx1), SXX(Ny1, Nx1), SXX0(Ny1, Nx1), SYY(Ny1, Nx1), SYY0(Ny1, Nx1), DILP(Ny1, Nx1);
*/

std::unique_ptr<gko::matrix::Dense<>> RHOX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* RHOX = RHOX_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> RHOFX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* RHOFX = RHOFX_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> ETADX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* ETADX = ETADX_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> PORX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* PORX = PORX_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> VX0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* VX0 = VX0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> VXF0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* VXF0 = VXF0_gko->get_values(); // Ny1 x Nx1
/*
// Vx nodes
Eigen::MatrixXd RHOX(Ny1, Nx1), RHOFX(Ny1, Nx1), ETADX(Ny1, Nx1), PORX(Ny1, Nx1), VX0(Ny1, Nx1), VXF0(Ny1, Nx1);
*/

std::unique_ptr<gko::matrix::Dense<>> RHOY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* RHOY = RHOY_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> RHOFY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* RHOFY = RHOFY_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> ETADY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* ETADY = ETADY_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> PORY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* PORY = PORY_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> VY0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* VY0 = VY0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> VYF0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* VYF0 = VYF0_gko->get_values(); // Ny1 x Nx1
/*
// Vy nodes
Eigen::MatrixXd RHOY(Ny1, Nx1), RHOFY(Ny1, Nx1), ETADY(Ny1, Nx1), PORY(Ny1, Nx1), VY0(Ny1, Nx1), VYF0(Ny1, Nx1);
*/

std::unique_ptr<gko::matrix::Dense<>> VZ0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* VZ0 = VZ0_gko->get_values(); // Ny1 x Nx1

// Vz nodes
//Eigen::MatrixXd VZ0(Ny1, Nx1);


std::unique_ptr<gko::matrix::Dense<>> ESP_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* ESP = ESP_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> EXY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* EXY = EXY_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> EZX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* EZX = EZX_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> EZY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* EZY = EZY_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> EXX_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* EXX = EXX_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> EYY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* EYY = EYY_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> EII_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* EII = EII_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> EIIVP_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* EIIVP = EIIVP_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> SII_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* SII = SII_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> DSII_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* DSII = DSII_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> DIS_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* DIS = DIS_gko->get_values(); // Ny1 x Nx1

std::unique_ptr<gko::matrix::Dense<>> EL_DECOM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1, Elastic (De)compaction
double* EL_DECOM = EL_DECOM_gko->get_values(); // Ny1 x Nx1, Elastic (De)compaction
std::unique_ptr<gko::matrix::Dense<>> VIS_COMP_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* VIS_COMP = VIS_COMP_gko->get_values(); // Ny1 x Nx1


std::unique_ptr<gko::matrix::Dense<>> t_marker_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* t_marker = t_marker_gko->get_values(); // marknum x 1 vector, Marker rock type
std::unique_ptr<gko::matrix::Dense<>> rhom_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* rhom = rhom_gko->get_values(); // marknum x 1 vector, Density of solid
std::unique_ptr<gko::matrix::Dense<>> etasm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* etasm = etasm_gko->get_values(); // marknum x 1 vector, Standard shear viscosity of bulk
std::unique_ptr<gko::matrix::Dense<>> etam_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* etam = etam_gko->get_values(); // marknum x 1 vector, Shear viscosity of bulk
std::unique_ptr<gko::matrix::Dense<>> cohescm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* cohescm = cohescm_gko->get_values(); // marknum x 1 vector
std::unique_ptr<gko::matrix::Dense<>> cohestm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* cohestm = cohestm_gko->get_values(); // marknumx1 vector
std::unique_ptr<gko::matrix::Dense<>> frictcm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* frictcm = frictcm_gko->get_values(); // marknumx1 vector,Marker rock type
std::unique_ptr<gko::matrix::Dense<>> dilatcm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* dilatcm = dilatcm_gko->get_values(); // marknumx1 vector, Marker rock type
std::unique_ptr<gko::matrix::Dense<>> fricttm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* fricttm = fricttm_gko->get_values(); // marknum x 1 vector Marker rock type
std::unique_ptr<gko::matrix::Dense<>> porm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* porm = porm_gko->get_values(); // marknum x 1 vector, Marker rock type
std::unique_ptr<gko::matrix::Dense<>> kkkm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* kkkm = kkkm_gko->get_values(); // // marknum x 1 vector, Marker rock type
std::unique_ptr<gko::matrix::Dense<>> xm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* xm = xm_gko->get_values(); // marknum x 1 vector, Horizontal coordinates of solid markers
std::unique_ptr<gko::matrix::Dense<>> ym_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* ym = ym_gko->get_values(); // marknum x 1 vector, Vertical coordinates of solid markers
std::unique_ptr<gko::matrix::Dense<>> sxxm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* sxxm = sxxm_gko->get_values(); // marknum x 1 vector, Marker SIGMAxx', Pa
std::unique_ptr<gko::matrix::Dense<>> syym_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* syym = syym_gko->get_values(); // marknum x 1 vector, Marker SIGMAyy', Pa
std::unique_ptr<gko::matrix::Dense<>> sxym_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum,1));
double* sxym = sxym_gko->get_values(); // marknum x 1 vector, Marker SIGMAxy', Pa



// (3) Defining global matrixes
// according to the global number of unknowns
// Sparse Matrix L is not yet defined as it will be built from a set of Triplets each step
std::unique_ptr<gko::matrix::Dense<>> R_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(N,1)); // N x 1 Vector of the right parts of equations
double* R = R_gko->get_values();// N x 1 Vector of the right parts of equations
std::unique_ptr<gko::matrix::Dense<>> S_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(N,1)); // N x 1 Vector
double* S = S_gko->get_values(); // N x 1 vector

// Initializing the input and solution vector on the gpu
std::unique_ptr<gko::matrix::Dense<>> R_gpu = gko::matrix::Dense<double>::create(gpu_exec, gko::dim<2>(N,1)); // N x 1 Vector
std::unique_ptr<gko::matrix::Dense<>> S_gpu = gko::matrix::Dense<double>::create(gpu_exec, gko::dim<2>(N, 1));

// variable type declaration
int ynlast, iterstep;

double timesum = 0;
double dt = dtelastic0;
bool yndtdecrease = true;

double dt00, dtx, dty, dtz, dtlapusta, Vmax, maxvxy;

//niterglobal x 1
std::unique_ptr<gko::matrix::Dense<>> DSYLSQ_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(niterglobal,1));
double* DSYLSQ = DSYLSQ_gko->get_values();

std::unique_ptr<gko::matrix::Dense<>> DVX0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* DVX0 = DVX0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> DVY0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* DVY0 = DVY0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> DVZ0_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
double* DVZ0 = DVZ0_gko->get_values(); // Ny1 x Nx1
std::unique_ptr<gko::matrix::Dense<>> DSY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* DSY = DSY_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> YNY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* YNY = YNY_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SigmaY_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SigmaY = SigmaY_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SII_fault_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SII_fault = SII_fault_gko->get_values(); // Ny x Nx
std::unique_ptr<gko::matrix::Dense<>> SIIB_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
double* SIIB = SIIB_gko->get_values(); // Ny x Nx


std::unique_ptr<gko::matrix::Dense<>> timesumcur_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1)); // num_timesteps x 1 Vector
double* timesumcur = timesumcur_gko->get_values();
std::unique_ptr<gko::matrix::Dense<>> dtcur_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1));// num_timesteps x 1 Vector
double* dtcur = dtcur_gko->get_values();// num_timesteps x 1 Vector

std::unique_ptr<gko::matrix::Dense<>> maxvxsmod_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1));// num_timesteps x 1 Vector
double* maxvxsmod = maxvxsmod_gko->get_values();// num_timesteps x 1 Vector
std::unique_ptr<gko::matrix::Dense<>> minvxsmod_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1));// num_timesteps x 1 Vector
double* minvxsmod = minvxsmod_gko->get_values();// num_timesteps x 1 Vector
std::unique_ptr<gko::matrix::Dense<>> maxvysmod_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1));// num_timesteps x 1 Vector
double* maxvysmod = maxvysmod_gko->get_values();// num_timesteps x 1 Vector
std::unique_ptr<gko::matrix::Dense<>> minvysmod_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1));// num_timesteps x 1 Vector
double* minvysmod = minvysmod_gko->get_values();// num_timesteps x 1 Vector
std::unique_ptr<gko::matrix::Dense<>> maxvzsmod_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1));// num_timesteps x 1 Vector
double* maxvzsmod = maxvzsmod_gko->get_values();// num_timesteps x 1 Vector
std::unique_ptr<gko::matrix::Dense<>> minvzsmod_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(num_timesteps, 1));// num_timesteps x 1 Vector
double* minvzsmod = minvzsmod_gko->get_values();// num_timesteps x 1 Vector












// variables for timer 
std::chrono::time_point<std::chrono::system_clock> time_start, time_end; // total runtime

int main() {
    t_marker_gko->fill(1.0);
    rhom_gko->fill(2800.0);
    etasm_gko->fill(1e21);
    cohescm_gko->fill(cohes);
    cohestm_gko->fill(cohes);
    frictcm_gko->fill(0.5);
    dilatcm_gko->fill(dilatation);
    fricttm_gko->fill(tensile);
    kkkm_gko->fill(2e-16);
    gm_gko->fill(shearmod);
    rhofm_gko->fill(1000);
    etafm_gko->fill(1.0e-3);


    time_start = std::chrono::system_clock::now();
    std::srand(time(nullptr));


    // read input
    // Load file
    std::string timestep_str;
    std::ifstream input_timestep("./input_data/StartingTimestep.txt");
    std::getline(input_timestep, timestep_str);
    int timestep = stoi(timestep_str);
    input_timestep.close();
    init_geometry();
    std::cout << "The A " << Nx << " x " << Ny << " grid is simulated, starting at timestep = " << timestep << "!" << std::endl;

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
