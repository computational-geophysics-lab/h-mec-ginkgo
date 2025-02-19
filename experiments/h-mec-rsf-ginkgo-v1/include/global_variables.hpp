#ifndef GLOBALVAR_H
#define GLOBALVAR_H

#include <string>
#include <ginkgo/ginkgo.hpp>

#include "constants.hpp"
#include "ginkgo_exec.hpp"

using namespace std;

// Basic nodes
extern std::unique_ptr<gko::matrix::Dense<double>> OM0_gko; //Ny x Nx, Old State parameter
extern double*  OM0; // Ny x Nx, Old State parameter
extern std::unique_ptr<gko::matrix::Dense<double>> OM_gko; // State parameter
extern double*  OM; // Ny x Nx, State parameter
extern std::unique_ptr<gko::matrix::Dense<double>> OM5_gko; // Ny x Nx
extern double*  OM5; // Ny x Nx
extern std::unique_ptr<gko::matrix::Dense<double>> ARSF_gko;// Ny x Nx, a - parameter of RSF
extern double*  ARSF;// Ny x Nx, a - parameter of RSF
extern std::unique_ptr<gko::matrix::Dense<double>> BRSF_gko;// Ny x Nx, b - parameter of RSF
extern double*  BRSF;// Ny x Nx, b - parameter of RSF
extern std::unique_ptr<gko::matrix::Dense<double>> LRSF_gko;// Ny x Nx, L - parameter of RSF
extern double*  LRSF;// Ny x Nx, L - parameter of RSF

// Unknown parameters
extern std::unique_ptr<gko::matrix::Dense<double>> pt_gko;   //Ny1 x Nx1, Total pressure
extern double*  pt;   // Ny1 x Nx1, Total pressure
extern std::unique_ptr<gko::matrix::Dense<double>> vxs_gko;  // Ny1 x Nx1, Solid vx - velocity
extern double*  vxs;  // Ny1 x Nx1, Solid vx - velocity
extern std::unique_ptr<gko::matrix::Dense<double>> vys_gko; // Ny1 x Nx1, Solid vy - velocity
extern double*  vys;  // Ny1 x Nx1, Solid vy - velocity
extern std::unique_ptr<gko::matrix::Dense<double>> vzs_gko;  // Ny1 x Nx1, Solid vz - velocity
extern double*  vzs;  // Ny1 x Nx1, Solid vz - velocity
extern std::unique_ptr<gko::matrix::Dense<double>> pf_gko;  // Ny1 x Nx1, Fluid pressure
extern double*  pf;   // Ny1 x Nx1, Fluid pressure
extern std::unique_ptr<gko::matrix::Dense<double>> vxD_gko;  // Ny1 x Nx1, Darsi vx - velocity
extern double*  vxD;  // Ny1 x Nx1, Darsi vx - velocity
extern std::unique_ptr<gko::matrix::Dense<double>> vyD_gko;// Ny1 x Nx1, Darsi vy - velocity
extern double*  vyD;  // Ny1 x Nx1, Darsi vy - velocity




// Nodal matrices
// Basic nodes
extern std::unique_ptr<gko::matrix::Dense<double>> RHO_gko, ETA_gko, ETA0_gko, ETA1_gko, ETA5_gko, ETA00_gko, IETAPLB_gko, SXY_gko, SXY0_gko, SZX0_gko, SZY0_gko, YNY0_gko, KKK_gko,
      GGG_gko, COHC_gko, COHT_gko, FRIC_gko, FRIT_gko, DILC_gko, TTT_gko, EIIB_gko, VSLIPB_gko, SZX_gko, SZY_gko;// Ny x Nx
extern double *RHO, *ETA, *ETA0, *ETA1, *ETA5, *ETA00, *IETAPLB, *SXY, *SXY0, *SZX0, *SZY0, *YNY0, *KKK,
      *GGG, *COHC, *COHT, *FRIC, *FRIT, *DILC, *TTT, *EIIB, *VSLIPB, *SZX, *SZY;// Ny x Nx

// Pressure nodes
extern std::unique_ptr<gko::matrix::Dense<double>> ETAB_gko, ETAB0_gko, ETAP_gko, ETAP0_gko, POR_gko, GGGP_gko, GGGB_gko, PTF0_gko, PT0_gko, PF0_gko, pt_ave_gko,
      pf_ave_gko, SXX_gko, SXX0_gko, SYY_gko, SYY0_gko, DILP_gko; // Ny1 x Nx1
extern double *ETAB, *ETAB0, *ETAP, *ETAP0, *POR, *GGGP, *GGGB, *PTF0, *PT0, *PF0, *pt_ave,
      *pf_ave, *SXX, *SXX0, *SYY, *SYY0, *DILP; // Ny1 x Nx1

// Vx nodes
extern std::unique_ptr<gko::matrix::Dense<double>> RHOX_gko, RHOFX_gko, ETADX_gko, PORX_gko, VX0_gko, VXF0_gko; // Ny1 x Nx1
extern double *RHOX, *RHOFX, *ETADX, *PORX, *VX0, *VXF0; // Ny1 x Nx1

// Vy nodes
extern std::unique_ptr<gko::matrix::Dense<double>> RHOY_gko, RHOFY_gko, ETADY_gko, PORY_gko, VY0_gko, VYF0_gko; // Ny1 x Nx1
extern double *RHOY, *RHOFY, *ETADY, *PORY, *VY0, *VYF0; // Ny1 x Nx1
// Vz nodes
extern std::unique_ptr<gko::matrix::Dense<double>> VZ0_gko; // Ny1 x Nx1
extern double*  VZ0; // Ny1 x Nx1

extern std::unique_ptr<gko::matrix::Dense<double>> ESP_gko, EXY_gko, EZX_gko, EZY_gko; // Ny x Nx
extern double *ESP, *EXY, *EZX, *EZY; // Ny x Nx
extern std::unique_ptr<gko::matrix::Dense<double>> EXX_gko, EYY_gko, EII_gko, EIIVP_gko, SII_gko, DSII_gko, DIS_gko; // Ny1 x Nx1
extern double *EXX, *EYY, *EII, *EIIVP, *SII, *DSII, *DIS; // Ny1 x Nx1

extern std::unique_ptr<gko::matrix::Dense<double>> EL_DECOM_gko; // Ny1 x Nx1, Elastic (de)compaction
extern double*  EL_DECOM;   // Ny1 x Nx1, Elastic (de)compaction
extern std::unique_ptr<gko::matrix::Dense<double>> VIS_COMP_gko;
extern double*  VIS_COMP;

// Lagrangian solid markers
extern std::unique_ptr<gko::matrix::Dense<double>> t_marker_gko; // marknum x 1 vector, Marker rock type
extern double*  t_marker; // marknum x 1 vector, Marker rock type
extern std::unique_ptr<gko::matrix::Dense<double>> rhom_gko; // marknum x 1 vector, Density of solid
extern double*  rhom; // marknum x 1 vector, Density of solid
extern std::unique_ptr<gko::matrix::Dense<double>> etasm_gko; // marknum x 1 vector, Standard shear viscosity of bulk
extern double*  etasm; // marknum x 1 vector, Standard shear viscosity of bulk
extern std::unique_ptr<gko::matrix::Dense<double>> etam_gko; // marknum x 1 vector, Shear viscosity of bulk
extern double*  etam; // marknum x 1 vector, Shear viscosity of bulk
extern std::unique_ptr<gko::matrix::Dense<double>> cohescm_gko; // marknum x 1 vector
extern double*  cohescm; // marknum x 1 vector
extern std::unique_ptr<gko::matrix::Dense<double>> cohestm_gko; //marknum x 1 vector
extern double*  cohestm; //marknum x 1 vector
extern std::unique_ptr<gko::matrix::Dense<>> frictcm_gko;//marknum x 1 vector
extern double*  frictcm; //marknum x 1 vector
extern std::unique_ptr<gko::matrix::Dense<double>> dilatcm_gko;//marknum x 1 vector
extern double*  dilatcm; //marknum x 1 vector
extern std::unique_ptr<gko::matrix::Dense<double>> fricttm_gko;//marknum x 1 vector
extern double*  fricttm; //marknum x 1 vector
extern std::unique_ptr<gko::matrix::Dense<double>> porm_gko;//marknum x 1 vector
extern double*  porm; //marknum x 1 vector
extern std::unique_ptr<gko::matrix::Dense<double>> kkkm_gko;//marknum x 1 vector
extern double* kkkm; //marknum x 1 vector
extern std::unique_ptr<gko::matrix::Dense<double>> xm_gko;//marknum x 1 vector, Horizontal coordinates of solid markers
extern double*  xm; //marknum x 1 vector, Horizontal coordinates of solid markers
extern std::unique_ptr<gko::matrix::Dense<double>> ym_gko;//marknum x 1 vector, Vertical coordinates of solid markers
extern double*  ym; //marknum x 1 vector, Vertical coordinates of solid markers
extern std::unique_ptr<gko::matrix::Dense<double>> sxxm_gko; // marknum x 1 vector, Marker SIGMAxx', Pa
extern double*  sxxm; //marknum x 1 vector, Marker SIGMAxx', Pa
extern std::unique_ptr<gko::matrix::Dense<double>> syym_gko; // marknum x 1 vector, Marker SIGMAyy', Pa
extern double*  syym; //marknum x 1 vector, Marker SIGMAyy', Pa
extern std::unique_ptr<gko::matrix::Dense<double>> sxym_gko; // marknum x 1 vector, Marker SIGMAxy', Pa
extern double*  sxym; // marknum x 1 vector, Marker SIGMAxy', Pa

// (3) Defining global matrixes
// according to the global number of unknowns
// Sparse Matrix L is not yet defined as it will be built from a set of Triplets each step
extern std::unique_ptr<gko::matrix::Dense<double>> R_gko; // N x 1 Vector of the right parts of equations
extern double*  R; // N x 1 Vector of the right parts of equations
extern std::unique_ptr<gko::matrix::Dense<double>> S_gko; // N x 1
extern double*  S; // N x 1
extern std::unique_ptr<gko::matrix::Dense<double>> L_dense_gko; // N x N, Dense matrix to assemble the stencil. Then converted to Csr called L that is used for the solver
extern double* L_dense; // N x N, Dense matrix to assemble the stencil. Then converted to Csr matrix called L that is used to build the solver
extern std::unique_ptr<gko::matrix::Dense<double>> R_gpu;
extern std::unique_ptr<gko::matrix::Dense<double>> S_gpu;

// variable type declaration
extern int ynlast, iterstep;

extern double timesum;
extern double dt;
extern bool yndtdecrease;

extern double dt00, dtx, dty, dtz, dtlapusta, Vmax, maxvxy;

extern std::unique_ptr<gko::matrix::Dense<double>> DSYLSQ_gko;
extern double *DSYLSQ;

extern std::unique_ptr<gko::matrix::Dense<double>> DVX0_gko, DVY0_gko, DVZ0_gko, DSY_gko; // Ny1 x Nx1
extern double *DVX0, *DVY0, *DVZ0, *DSY;
extern std::unique_ptr<gko::matrix::Dense<double>> YNY_gko, SigmaY_gko, SII_fault_gko, SIIB_gko; // Ny x Nx
extern double *YNY, *SigmaY, *SII_fault, *SIIB; // Ny x Nx

extern std::unique_ptr<gko::matrix::Dense<double>> timesumcur_gko, dtcur_gko; // num_timesteps x 1 vector
extern double *timesumcur, *dtcur; // num_timesteps x 1 vector
extern std::unique_ptr<gko::matrix::Dense<double>> maxvxsmod_gko, minvxsmod_gko, maxvysmod_gko, minvysmod_gko, maxvzsmod_gko, minvzsmod_gko;
extern double *maxvxsmod, *minvxsmod, *maxvysmod, *minvysmod, *maxvzsmod, *minvzsmod;


#endif
