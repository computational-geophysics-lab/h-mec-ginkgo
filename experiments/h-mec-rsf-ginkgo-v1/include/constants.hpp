#ifndef CONST_H
#define CONST_H

#include <string>
#include <ginkgo/ginkgo.hpp>

#include "ginkgo_exec.hpp"

using namespace std;


// set timestep limit
const int num_timesteps = 50000;    // max number of timesteps
const int savestep = 1000;           // storage periodicity
const double significant_digits = 1e-10; // Defines the relative number of digits for the safe_difference function

const bool antiplane = false;
const int Num_var = 6; // set to 7 in antiplane = true

// ========================================
// Define Numerical model
// Eulerian basic grid
const double xsize = 40000.; // size in horizontal direction, m
const double ysize = 10000.; // size in vertical direction, m
const int Nx = 401; //Originally Nx = 401, Ny = 51          // number of grid steps in horizontal directions
const int Ny = 51;           // number of grid steps in vertical direction

// Where to apply the transition on the left(1) and right(2)
const double TS_1 = 2e3;
const double TS_2 = 4e3;

const double TS_3 = 34e3;
const double TS_4 = 38e3;

// Eulerian Staggered Grid
const int Nx1 = Nx + 1;     // Number of horizontal lines for staggered grid
const int Ny1 = Ny + 1;     // Number of vertical lines for staggered grid

const int N = Nx1 * Ny1 * Num_var; // Global number of unknowns

// ========================================
// Output files
const int line_fault = (Ny - 1) / 2.; 

// ========================================
// Coordinates
const double dx = xsize / (Nx - 1); // grid step in horizontal direction, m
const double dy = ysize / (Ny - 1); // grid step in vertical direction, m
const double xbeg = 0;
const double xend = xsize;
const double ybeg = 0;
const double yend = ysize;

const double dx2 = pow(dx, 2), dy2 = pow(dy, 2);
const double dx_dy = dx * dy;

// Returns a std::__1::unique_ptr<gko::matrix::Dense<gko::default_precision>>, i.e. a gko::matrix::Dense<double> with linearly spaced entries in a interval
inline std::unique_ptr<gko::matrix::Dense<>> GinkgoLinSpaced(const int Nsteps, const double start, const double end) {
    std::unique_ptr<gko::matrix::Dense<>> ginkgo_vec=gko::matrix::Dense<double>::create(exec, gko::dim<2>(Nsteps,1));
    double* ginkgo_vec_data = ginkgo_vec->get_values();
    const double step = (end-start)/(Nsteps-1);
    for (int i = 0; i<Nsteps; i++){
        ginkgo_vec_data[i] = start + i*step;
    }
    return ginkgo_vec;
}


const double faultwidth = dx;    // Characteristic fault width, m

// Entering time cycle
const double dtelastic0 = 5e8; // elastic timestep
const double dtmin = 1e-4;

const double inertia = 1.;

const double lower_block = ysize / 2. + dy;
const double upper_block = ysize / 2. - dy;

extern const std::unique_ptr<gko::matrix::Dense<>> x_gko;
extern const double* x;
extern const std::unique_ptr<gko::matrix::Dense<>> y_gko;
extern const double* y;
extern const std::unique_ptr<gko::matrix::Dense<>> xvx_gko;
extern const double* xvx;
extern const std::unique_ptr<gko::matrix::Dense<>> yvx_gko;
extern const double* yvx;
extern const std::unique_ptr<gko::matrix::Dense<>> xvy_gko;
extern const double* xvy;
extern const std::unique_ptr<gko::matrix::Dense<>> yvy_gko;
extern const double* yvy;
extern const std::unique_ptr<gko::matrix::Dense<>> xp_gko;
extern const double* xp;
extern const std::unique_ptr<gko::matrix::Dense<>> yp_gko;
extern const double* yp;
// ========================
// Lagrangian solid markers
const int Nx_markers = (Nx - 1) * 4; // Marker resolution in x - dection
const int Ny_markers = (Ny - 1) * 4; // Marker resolution in y direction
const double dxms = xsize / (double)Nx_markers; // Standard marker horizontal step
const double dyms = ysize / (double)Ny_markers; // Standard marker vertical step
const int marknum = Nx_markers * Ny_markers; // Total number of markers

// Plastic Strain:
// const double gammap = plstrain * 200;

const double pi = M_PI;

//                      Block  Fault
inline double arsfm[2] = {0.025,0.006};// 2x1 vector, a - parameter of RSF
inline double brsfm[2] = {0.001, 0.015};// 2x1 vector, b - parameter of RSF
const double lrsfm[2] = {0.020, 0.0085};// 2x1 vector, L - parameter of RSF
const double omm[2] = {10, -10}; // 2x1 vector, State


const double V0 = 1.e-9;               // Reference slip velocity of RSF, m / s
const int alpha = 29;                  // Scaling factor for shear viscosity of porous matrix

// ========================================
const double POR0 = .01; // Standard porosity

// ========================================
// Brittle / plastic rheology
const double cohes = 0.;         // Cohesion, Pa
const double friction = .6;      // unused variable ???
const double dilatation = 0.;    // Dilatation coefficient confined
const double tensile = 1.;       // Internal friction coefficient tensile
const double shearmod = 3.e10;   // Shear modulus
const double BETAFLUID = 1.e-8;  // 4.0e-10; // Compressibility of fluid, 1 / Pa
const double BETASOLID = 2.e-11; // 2.5e-11; // Compressibility of solid, 1 / Pa

// ========================================
// Constants
const double gx = 0.;        // Horizontal gravity, m / s^2
const double gy = 0.;        // Vertical gravity, m / s^2
const double PCONF = 1.e7;   // Confining pressure
const double PTFDIFF = 3.e7; // Total - Fluid Pressure difference in the top row, Pa

// ========================================
// Limits
const double etamin = 1e-3;  // Lower shear viscosity cutoff
const double etamax = 1e50;  // Upper shear viscosity cutoff
const double kkkmin = 1e-22; // Lower Darsi viscosity cutoff
const double kkkmax = 1e-12; // Upper Darsi viscosity cutoff
const double stpmax = 2e-4;  // / dy * faultwidth; // Max gridstep fraction for marker displacement in the channel
const double stpmax1 = 6e-5; // / dy * faultwidth; // Max gridstep fraction for marker displacement

//Boundary conditions
const double bcupper = -1e-9;
const double bclower = 1e-9;
const double bcvyflower = 0; // -1e-12;

const string nname = "h_mec_";  // storage filename
const int niterglobal = 10000; // Max number of global iterations
const int ynlastmax = 499; // Max number of iterations without changing timestep
const int dtstep = 5; // Min number of plastic iterations before changing dt
const double vratiomax = 0.001; // Max change of velocity between iterations
const double dtkoef = 1.1; // Koefficient for decrease of previous timestep
const double dtkoefup = 1.2; // Koefficient for increase of previous timestep
const double dtkoefv = 1.001; // koef for velocity - based timestep reduction
const double errmin = 1e2; // min LSQ err for stoping iterations
const double etawt = 0.0;
const double syieldmin = 1e-3;
const double tyield = 1;
const double plstrain = 0; // Plastic strain for dilatancy

// These variables are only used once and could be inserted directly as numbers
const double aa1 = 20.93, aa2 = 35.28, aa3 = 2.34;
const double bb1 = 0.99, bb2 = 44.39, bb3 = 0.73;
const double cc1 = 0.37, cc2 = 3.54, cc3 = 0.47;


inline std::unique_ptr<gko::matrix::Dense<>> gm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum, 1));
inline std::unique_ptr<gko::matrix::Dense<>> rhofm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum, 1));
inline std::unique_ptr<gko::matrix::Dense<>> etafm_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(marknum, 1));

inline double* gm = gm_gko->get_values();
inline double* etafm = etafm_gko->get_values();
inline double* rhofm = rhofm_gko->get_values();

// Lagrangian solid markers
// moved to global_variables.hpp

const double dy_faultw = dy / faultwidth;

#endif