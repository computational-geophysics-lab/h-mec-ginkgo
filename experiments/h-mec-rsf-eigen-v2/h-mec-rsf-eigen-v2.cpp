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
Fluid X - Darsi:  - ETAfluid / K * VxD - dPf / dx = -RHOf * gx
Fluid Y - Darsi:  - ETAfluid / K * VyD - dPf / dy = -RHOf * gy
Fluid Continuity: dVxD / dx + dVyD / dy - (Pt - Pf) / ETAbulk = 0
+ staggered grid
+ P - v formulation
+ advecting material fileds by markers
========================================

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
to do:

1. Change variable names to useful / understandable names

2. Move all variable declarations outside of loops if possible

3. Do calculations outside of loops if possible
  a) Multiplications of every element can be done at once outside of loop
  
4. Remove unnessesary calculations if there are any
  a) This includes copying variables not needed
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse> // should be already included with the last include but somehow Euler needs this include too
#include <random>
#include <ctime>
#include <math.h>
#include <iomanip>
#include <chrono>

// using namespace Eigen;

// ====================================================================================
// global variable declarations
// ====================================================================================

// set timestep limit
const int num_timesteps = 1; // very small number for testing

// ========================================
// Define Numerical model
// Eulerian basic grid
/*
const double xsize = 40000.;      // size in horizontal direction, m
const double ysize = 10000.;      // size in vertical direction, m
const int Nx = 161;               // number of grid steps in horizontal directions
const int Ny = 41;                // number of grid steps in vertical direction

// Where to apply the transition on the left(1) and right(2)
const double TS_1 = 6e3;
const double TS_2 = 8e3;

const double TS_3 = 34e3;
const double TS_4 = 37e3;
*/

// very tiny test system
const double xsize = 10000.;      // size in horizontal direction, m
const double ysize = 2500.;      // size in vertical direction, m
const int Nx = 41;               // number of grid steps in horizontal directions
const int Ny = 11;

// Where to apply the transition on the left(1) and right(2)
const double TS_1 = 6e3 / 4.;
const double TS_2 = 8e3 / 4.;

const double TS_3 = 34e3 / 4.;
const double TS_4 = 37e3 / 4.;

// Eulerian Staggered Grid
const int Nx1 = Nx + 1;           // Number of horizontal lines for staggered grid
const int Ny1 = Ny + 1;           // Number of vertical lines for staggered grid

const int N = Nx1 * Ny1 * 6;      // Global number of unknowns

// ========================================
// Output files
const bool seismic_cycles_data = true; // unnessesary
int timesum_plus = 0;
const int line_fault = (Ny - 1) / 2.; 

// ========================================
// Coordinates
const double dx = xsize / (Nx - 1); // grid step in horizontal direction, m
const double dy = ysize / (Ny - 1); // grid step in vertical direction, m
const double xbeg = 0;
const double xend = xsize;
const double ybeg = 0;
const double yend = ysize;

// define type and replace : operator
Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Nx, xbeg, xend); // horizontal coordinates of basic grid points
Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(Ny, ybeg, yend); // vertical coordinates of basic grid points
Eigen::VectorXd xvx = Eigen::VectorXd::LinSpaced(Nx1, xbeg, xend + dx); // Horizontal coordinates of Vx - nodes
Eigen::VectorXd yvx = Eigen::VectorXd::LinSpaced(Ny1, ybeg - dy / 2., yend + dy / 2.); // Vertical coordinates of Vx - nodes
Eigen::VectorXd xvy = Eigen::VectorXd::LinSpaced(Nx1, xbeg - dx / 2., xend + dx / 2.); // Horizontal coordinates of Vy - nodes
Eigen::VectorXd yvy = Eigen::VectorXd::LinSpaced(Ny1, ybeg, yend + dy); // Vertical coordinates of Vy - nodes
Eigen::VectorXd xp = Eigen::VectorXd::LinSpaced(Nx1, xbeg - dx / 2., xend + dx / 2.); // Horizontal coordinates of P - nodes
Eigen::VectorXd yp = Eigen::VectorXd::LinSpaced(Ny1, ybeg - dy / 2., yend + dy / 2.); // Vertical coordinates of Vx - nodes

//               Block  Fault
Eigen::Vector2d arsfm = {.03,  .009}; // a - parameter of RSF
Eigen::Vector2d brsfm = {.001, .015}; // b - parameter of RSF
Eigen::Vector2d lrsfm = {.014, .014}; // L - parameter of RSF (characteristic slip distance)
Eigen::Vector2d omm   = {10,  -5   }; // State
double V0 = 1.e-9;                   // Reference slip velocity of RSF, m / s
int alpha = 29;                      // Scaling factor for shear viscosity of porous matrix

// ========================================
const double POR0 = .01; // Standard porosity

// ========================================
// Brittle / plastic rheology
const double cohes = 0.;     // Cohesion, Pa
const double friction = .6;     // Internal friction coefficient confined
const double dilatation = 0.;  // Dilatation coefficient confined
const double tensile = 1.;      // Internal friction coefficient tensile
const double shearmod = 3.e10; // Shear modulus
const double BETAFLUID = 1.e-8;  // 4.0e-10; // Compressibility of fluid, 1 / Pa
const double BETASOLID = 2.e-11; // 2.5e-11; // Compressibility of solid, 1 / Pa
const double faultwidth = dx;    // Characteristic fault width, m

// ========================================
// Constants
const double gx = 0.;             // Horizontal gravity, m / s^2
const double gy = 0.;             // Vertical gravity, m / s^2
const double PCONF = 1.e7;     // Confining pressure
const double PTFDIFF = 3.e7;    // Total - Fluid Pressure difference in the top row, Pa

// ========================================
// Limits
const double etamin = 1e-3; // Lower shear viscosity cutoff
const double etamax = 1e50; // Upper shear viscosity cutoff
const double kkkmin = 1e-22; // Lower Darsi viscosity cutoff
const double kkkmax = 1e-12; // Upper Darsi viscosity cutoff
const double stpmax = 2e-4; // / dy * faultwidth; // Max gridstep fraction for marker displacement in the channel
const double stpmax1 = 6e-5; // / dy * faultwidth; // Max gridstep fraction for marker displacement

//Boundary conditions
const double bcupper = -2e-9;
const double bclower = 2e-9;
const double bcvyflower = 0;// - 1e-12;
// bcvxfleft = bcvyflower * xsize / ysize;

// Entering time cycle
double timesum = 0;
const double dtelastic0 = 5e8; // elastic timestep
double dt = dtelastic0;
const double dtmin = 1e-4;

const double ascale = 1.;

// Timesteps between visualization frames
const int savematstep = 300;  // storage periodicity
std::string nname = "h_mec_";  // storage filename
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
bool yndtdecrease = true;
double plstrain = 0; // Plastic strain for dilatancy
int iterstep;

const double lower_block = ysize / 2. + dy;
const double upper_block = ysize / 2. - dy;

// ========================
// Lagrangian solid markers
const int Nx_markers = (Nx - 1) * 4; // Marker resolution in x - dection
const int Ny_markers = (Ny - 1) * 4; // Marker resolution in y direction
const double dxms = xsize / (double)Nx_markers; // Standard marker horizontal step
const double dyms = ysize / (double)Ny_markers; // Standard marker vertical step
const int marknum = Nx_markers * Ny_markers; // Total number of markers

// These variables are only used once and could be inserted directly as numbers
const double aa1 = 20.93, aa2 = 35.28, aa3 = 2.34;
const double bb1 = 0.99, bb2 = 44.39, bb3 = 0.73;
const double cc1 = 0.37, cc2 = 3.54, cc3 = 0.47;

double aa, bb, cc;
double ss3;

// Plastic Strain:
const double gammapij = plstrain * 200;
double dilij;

const double pi = M_PI;

// Basic nodes
Eigen::MatrixXd OM0(Ny, Nx);    // Old state parameter
Eigen::MatrixXd OM(Ny, Nx);     // State parameter
Eigen::MatrixXd OM5(Ny, Nx);
Eigen::MatrixXd ARSF(Ny, Nx); // a - parameter of RSF
Eigen::MatrixXd BRSF(Ny, Nx); // b - parameter of RSF
Eigen::MatrixXd LRSF(Ny, Nx); // L - parameter of RSF

// Unknown parameters
Eigen::MatrixXd pt(Ny1, Nx1);  // Total pressure
Eigen::MatrixXd vxs(Ny1, Nx1); // Solid vx - velocity
Eigen::MatrixXd vys(Ny1, Nx1); // Solid vy - velocity
Eigen::MatrixXd pf(Ny1, Nx1);  // Fluid pressure
Eigen::MatrixXd vxD(Ny1, Nx1); // Darsi vx - velocity
Eigen::MatrixXd vyD(Ny1, Nx1); // Darsi vy - velocity

// Nodal arrays
// Basic nodes
Eigen::MatrixXd RHO(Ny, Nx);
Eigen::MatrixXd ETA(Ny, Nx);
Eigen::MatrixXd ETA0(Ny, Nx);
Eigen::MatrixXd ETA1(Ny, Nx);
Eigen::MatrixXd ETA5(Ny, Nx);
Eigen::MatrixXd ETA50(Ny, Nx);
Eigen::MatrixXd ETA00(Ny, Nx);
Eigen::MatrixXd IETAPLB(Ny, Nx);
Eigen::MatrixXd SXY(Ny, Nx);
Eigen::MatrixXd SXY0(Ny, Nx);
Eigen::MatrixXd YNY0(Ny, Nx);
Eigen::MatrixXd YNY00(Ny, Nx);
Eigen::MatrixXd KKK(Ny, Nx);
Eigen::MatrixXd GGG(Ny, Nx);
Eigen::MatrixXd COHC(Ny, Nx);
Eigen::MatrixXd COHT(Ny, Nx);
Eigen::MatrixXd FRIC(Ny, Nx);
Eigen::MatrixXd FRIT(Ny, Nx);
Eigen::MatrixXd DILC(Ny, Nx);
Eigen::MatrixXd TTT(Ny, Nx);
Eigen::MatrixXd EIIB(Ny, Nx);
Eigen::MatrixXd STRPLB(Ny, Nx);

// Pressure nodes
Eigen::MatrixXd ETAB(Ny1, Nx1);
Eigen::MatrixXd ETAB0(Ny1, Nx1);
Eigen::MatrixXd ETAP(Ny1, Nx1);
Eigen::MatrixXd ETAP0(Ny1, Nx1);
Eigen::MatrixXd POR(Ny1, Nx1);
Eigen::MatrixXd GGGP(Ny1, Nx1);
Eigen::MatrixXd GGGB(Ny1, Nx1);
Eigen::MatrixXd PTF0(Ny1, Nx1);
Eigen::MatrixXd PT0(Ny1, Nx1);
Eigen::MatrixXd PF0(Ny1, Nx1);
Eigen::MatrixXd pt_ave(Ny1, Nx1);
Eigen::MatrixXd pf_ave(Ny1, Nx1);
Eigen::MatrixXd SXX(Ny1, Nx1);
Eigen::MatrixXd SXX0(Ny1, Nx1);
Eigen::MatrixXd SYY(Ny1, Nx1);
Eigen::MatrixXd SYY0(Ny1, Nx1);
Eigen::MatrixXd DILP(Ny1, Nx1);
// Vx nodes
Eigen::MatrixXd RHOX(Ny1, Nx1);
Eigen::MatrixXd RHOFX(Ny1, Nx1);
Eigen::MatrixXd ETADX(Ny1, Nx1);
Eigen::MatrixXd PORX(Ny1, Nx1);
Eigen::MatrixXd VX0(Ny1, Nx1);
Eigen::MatrixXd VXF0(Ny1, Nx1);
// Vy nodes
Eigen::MatrixXd RHOY(Ny1, Nx1);
Eigen::MatrixXd RHOFY(Ny1, Nx1);
Eigen::MatrixXd ETADY(Ny1, Nx1);
Eigen::MatrixXd PORY(Ny1, Nx1);
Eigen::MatrixXd VY0(Ny1, Nx1);
Eigen::MatrixXd VYF0(Ny1, Nx1);

Eigen::MatrixXd VSLIPB(Ny, Nx);

// Lagrangian solid markers
Eigen::VectorXd rhom(marknum);         // Density of solid
Eigen::VectorXd etasm(marknum);        // Standard shear viscosity of bulk
Eigen::VectorXd etam(marknum);         // Shear viscosity of bulk
Eigen::VectorXd etabm(marknum);        // Bulk viscosity of bulk
Eigen::VectorXd cohescm(marknum);      // Cohesion for confined fracture of solid
Eigen::VectorXd frictcm(marknum);      // friction for confined fracture of solid
Eigen::VectorXd dilatcm(marknum);      // dilatation for confined fracture of solid
Eigen::VectorXd cohestm(marknum);      // Cohesion for tensile fracture of solid
Eigen::VectorXd fricttm(marknum);      // friction for tensile fracture of solid
Eigen::VectorXd porm(marknum);         // Porosity of solid
Eigen::VectorXd kkkm(marknum);         // Standard permeability of solid
Eigen::VectorXd rhofm(marknum);        // Density of fluid
Eigen::VectorXd etafm(marknum);        // Viscosity of fluid
Eigen::VectorXd t_marker(marknum);           // Marker rock type
Eigen::VectorXd xm(marknum);           // Horizontal coordinates of solid markers
Eigen::VectorXd ym(marknum);           // Vertical coordinates of solid markers
Eigen::VectorXd sxxm(marknum);         // Marker SIGMAxx', Pa
Eigen::VectorXd syym(marknum);         // Marker SIGMAyy', Pa
Eigen::VectorXd sxym(marknum);         // Marker SIGMAxy', Pa
Eigen::VectorXd gsm(marknum);          // Standard shear modulus of bulk, Pa
Eigen::VectorXd gm(marknum);           // Shear modulus of bulk, Pa
Eigen::VectorXd vx0m(marknum);         // Marker horizontal velocity
Eigen::VectorXd vy0m(marknum);         // Marker vertical velocity
Eigen::VectorXd ptfm(marknum);         // Pt - Pf, Pa
Eigen::VectorXd amursfm(marknum);      // RSF a / mu parameter

// declaration of matrices that are set to zero each timestep
Eigen::MatrixXd ETA0SUM(Ny, Nx);
Eigen::MatrixXd COHCSUM(Ny, Nx);
Eigen::MatrixXd FRICSUM(Ny, Nx);
Eigen::MatrixXd DILCSUM(Ny, Nx);
Eigen::MatrixXd COHTSUM(Ny, Nx);
Eigen::MatrixXd FRITSUM(Ny, Nx);
Eigen::MatrixXd WTSUM(Ny, Nx);

// Interpolate ETA, RHO to nodal points
// Basic nodes
Eigen::MatrixXd RHOSUM(Ny, Nx);
Eigen::MatrixXd ETASUM(Ny, Nx);
Eigen::MatrixXd KKKSUM(Ny, Nx);
Eigen::MatrixXd TTTSUM(Ny, Nx);
Eigen::MatrixXd SXYSUM(Ny, Nx);
Eigen::MatrixXd GGGSUM(Ny, Nx);

//LDZ
Eigen::MatrixXd OM0SUM(Ny, Nx); // Old state parameter
Eigen::MatrixXd OMSUM(Ny, Nx); // State parameter
Eigen::MatrixXd ARSFSUM(Ny, Nx); // a - parameter of RSF
Eigen::MatrixXd BRSFSUM(Ny, Nx); // b - parameter of RSF
Eigen::MatrixXd LRSFSUM(Ny, Nx); // L - parameter of RSF

// Pressure nodes
Eigen::MatrixXd ETAPSUM(Ny1, Nx1);
Eigen::MatrixXd ETAP0SUM(Ny1, Nx1);
Eigen::MatrixXd ETAB0SUM(Ny1, Nx1);
Eigen::MatrixXd PORSUM(Ny1, Nx1);
Eigen::MatrixXd SXXSUM(Ny1, Nx1);
Eigen::MatrixXd SYYSUM(Ny1, Nx1);
Eigen::MatrixXd GGGPSUM(Ny1, Nx1);
Eigen::MatrixXd WTPSUM(Ny1, Nx1);
// Vx nodes
Eigen::MatrixXd RHOXSUM(Ny1, Nx1);
Eigen::MatrixXd RHOFXSUM(Ny1, Nx1);
Eigen::MatrixXd ETADXSUM(Ny1, Nx1);
Eigen::MatrixXd PORXSUM(Ny1, Nx1);
Eigen::MatrixXd VX0SUM(Ny1, Nx1);
Eigen::MatrixXd WTXSUM(Ny1, Nx1);
// Vy nodes
Eigen::MatrixXd RHOYSUM(Ny1, Nx1);
Eigen::MatrixXd RHOFYSUM(Ny1, Nx1);
Eigen::MatrixXd ETADYSUM(Ny1, Nx1);
Eigen::MatrixXd PORYSUM(Ny1, Nx1);
Eigen::MatrixXd VY0SUM(Ny1, Nx1);
Eigen::MatrixXd WTYSUM(Ny1, Nx1);

Eigen::MatrixXd ESP(Ny, Nx);
Eigen::MatrixXd EXY(Ny, Nx);
Eigen::MatrixXd DSXY(Ny, Nx);
Eigen::MatrixXd EXX(Ny1, Nx1);
Eigen::MatrixXd DSXX(Ny1, Nx1);
Eigen::MatrixXd EYY(Ny1, Nx1);
Eigen::MatrixXd DSYY(Ny1, Nx1);
Eigen::MatrixXd EII(Ny1, Nx1);
Eigen::MatrixXd EIIVP(Ny1, Nx1);
Eigen::MatrixXd SII(Ny1, Nx1);
Eigen::MatrixXd DSII(Ny1, Nx1);
Eigen::MatrixXd DIS(Ny1, Nx1);

Eigen::MatrixXd EL_DECOM(Ny1, Nx1);   // Elastic (de)compaction
Eigen::MatrixXd VIS_COMP(Ny1, Nx1);

// (3) Defining global matrixes
// according to the global number of unknowns
Eigen::SparseMatrix<double> L(N, N); // Matrix of coefficients in the left part
Eigen::VectorXd R(N); // Vector of the right parts of equations
Eigen::VectorXd S(N);

// variable type declaration
double cohescmm, cohestmm, frictcmm, dilatcmm, fricttmm, etasmm0, etamm0, etamm, rhomm, etadm, dxm, dym, wtmij, wtmi1j, wtmij1, wtmi1j1;

double dt00, etamincur;

int ynlast;

double KXX, IETAPL, BETADRAINED, KBW, KSK, KXY;

double pfscale, ptscale;

int kp, kx, ky, kpf, kxf, kyf; // this could maybe be made into a vector

double ETAXY1, ETAXY2, ETAXX1, ETAXX2, ETAYY1, ETAYY2, GXY1, GXY2, GXX1, GXX2, GYY1, GYY2, KXY1, KXY2, KXX1, KXX2, KYY1, KYY2, SXY1, SXY2, SXX1, SXX2, SYY1, SYY2;

double dRHOdx, dRHOdy;

double Vmax;

double SIIB1, SIIB2, SIIB3, SIIB4, SIIB5;

double avgpt, diffpt;

double dt0;

double EXY2, EXYVP2, DISXY;

double ptB, pfB, prB;
double kfxy, siiel, kfxy0, SIIB0;

double V, EIISLIP, ETAPL, ETAVP, kfxy1, DSIIB1, DSIIB2;

double syield, SIGMA2, etapl, dTETAmax, maxvxy;

double dtx, dty, dtlapusta;

double EXY1, PT0_ave, PF0_ave;

Eigen::VectorXd vxm(4), vym(4), spm(4);

Eigen::VectorXd DSYLSQ(niterglobal);

Eigen::MatrixXd DVX0(Ny1, Nx1);
Eigen::MatrixXd DVY0(Ny1, Nx1);

Eigen::MatrixXd AXY(Ny, Nx);
Eigen::MatrixXd DSY(Ny, Nx);
Eigen::MatrixXd YNY(Ny, Nx);
Eigen::MatrixXd SigmaY(Ny, Nx);
Eigen::MatrixXd SII_fault(Ny, Nx);

Eigen::MatrixXd SIIB(Ny, Nx);

Eigen::VectorXd timesumcur(num_timesteps);
Eigen::VectorXd dtcur(num_timesteps);
Eigen::VectorXd maxvxsmod(num_timesteps);
Eigen::VectorXd minvxsmod(num_timesteps);
Eigen::VectorXd maxvysmod(num_timesteps);
Eigen::VectorXd minvysmod(num_timesteps);

// ====================================================================================
// end global variable declarations
// ====================================================================================

// lambda function that rounds towards 0
auto fix = [](double x) {
  return x < 0. ? (int)ceil(x) : (int)floor(x);
};

int check_bounds(int k, int bound) {
    if (k < 0) {
        k = 0;
    } else if (k > bound - 2) {
        k = bound - 2;
    }
    return k;
}

int main() {
    // ====================================================================================
    // initialize random generator with set seed for testing purpose
    std::srand(42);
    // initialize random generator with current time as seed
    // std::srand(time(nullptr));
    // ====================================================================================

    // read input
    // Load file
    std::string timestep_str;
    std::ifstream input_timestep("data/file.txt");
    std::getline(input_timestep, timestep_str);
    int timestep = std::stoi(timestep_str);
    input_timestep.close();

    // no. time step
    std::cout << "Total no. time steps: num_timesteps == " << num_timesteps << std::endl;

    if (timestep > 0) {
        std::string filename = "h_mec_" + std::to_string(timestep) + ".txt";
        std::ifstream input_file(filename);
        // read input file here
        timestep++;
    } else {
        timestep = 1;
        // Basic nodes
        OM0 = Eigen::MatrixXd::Constant(Ny, Nx, omm(0));    // Old state parameter
        ARSF = Eigen::MatrixXd::Constant(Ny, Nx, arsfm(0)); // a - parameter of RSF
        BRSF = Eigen::MatrixXd::Constant(Ny, Nx, brsfm(0)); // b - parameter of RSF
        LRSF = Eigen::MatrixXd::Constant(Ny, Nx, lrsfm(0)); // L - parameter of RSF
        
        // Unknown parameters
        pt.setZero();  // Total pressure
        vxs.setZero(); // Solid vx - velocity
        vys.setZero(); // Solid vy - velocity
        pf.setZero();  // Fluid pressure
        vxD.setZero(); // Darsi vx - velocity
        vyD.setZero(); // Darsi vy - velocity
        
        // Nodal arrays
        // Basic nodes
        RHO.setZero();
        ETA.setZero();
        ETA0.setZero();
        IETAPLB.setZero();
        SXY.setZero();
        SXY0.setZero();
        YNY0.setZero();
        KKK.setZero();
        GGG.setZero();
        COHC.setZero();
        COHT.setZero();
        FRIC.setZero();
        FRIT.setZero();
        DILC.setZero();
        TTT.setZero();
        EIIB.setZero();
        STRPLB.setZero();
        
        // Define Fault
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (y(i) > upper_block && y(i) < lower_block) {
                    OM0(i, j) = omm(1);
                
                    if (x(j) >= TS_1 && x(j) < TS_2) {
                        BRSF(i, j) = brsfm(0) - (brsfm(0) - brsfm(1)) * ((x(j) - TS_1) / (TS_2 - TS_1));
                        ARSF(i, j) = arsfm(0) - (arsfm(0) - arsfm(1)) * ((x(j) - TS_1) / (TS_2 - TS_1));
                    }
                    if (x(j) >= TS_2 && x(j) <= TS_3) {
                        BRSF(i, j) = brsfm(1);
                        ARSF(i, j) = arsfm(1);
                    }
                    if (x(j) > TS_3 && x(j) <= TS_4) {
                        BRSF(i, j) = brsfm(1) - (brsfm(1) - brsfm(0)) * ((x(j) - TS_3) / (TS_4 - TS_3));
                        ARSF(i, j) = arsfm(1) - (arsfm(1) - arsfm(0)) * ((x(j) - TS_3) / (TS_4 - TS_3));
                    }
                }
            }
        }
        
        OM = OM0;
        
        // Pressure nodes
        ETAB.setZero();
        ETAB0.setZero();
        ETAP.setZero();
        ETAP0.setZero();
        POR.setZero();
        GGGP.setZero();
        GGGB.setZero();
        PTF0.setZero();
        PT0 = Eigen::MatrixXd::Constant(Ny1, Nx1, PCONF + PTFDIFF);
        PF0 = Eigen::MatrixXd::Constant(Ny1, Nx1, PCONF);
        pt_ave.setZero();
        pf_ave.setZero();
        SXX.setZero();
        SXX0.setZero();
        SYY.setZero();
        SYY0.setZero();
        DILP.setZero();
        // Vx nodes
        RHOX.setZero();
        RHOFX.setZero();
        ETADX.setZero();
        PORX.setZero();
        VX0.setZero();
        VXF0.setZero();
        // Vy nodes
        RHOY.setZero();
        RHOFY.setZero();
        ETADY.setZero();
        PORY.setZero();
        VY0.setZero();
        VYF0.setZero();
        
        VSLIPB.setZero();
        
        // Lagrangian solid markers
        rhom.setZero();        // Density of solid
        etasm.setZero();       // Standard shear viscosity of bulk
        etam.setZero();        // Shear viscosity of bulk
        etabm.setZero();       // Bulk viscosity of bulk
        cohescm.setZero();     // Cohesion for confined fracture of solid
        frictcm.setZero();     // friction for confined fracture of solid
        dilatcm.setZero();     // dilatation for confined fracture of solid
        cohestm.setZero();     // Cohesion for tensile fracture of solid
        fricttm.setZero();     // friction for tensile fracture of solid
        porm.setZero();        // Porosity of solid
        kkkm.setZero();        // Standard permeability of solid
        rhofm.setZero();       // Density of fluid
        etafm.setZero();       // Viscosity of fluid
        t_marker.setZero();    // Marker rock type
        xm.setZero();          // Horizontal coordinates of solid markers
        ym.setZero();          // Vertical coordinates of solid markers
        sxxm.setZero();        // Marker SIGMAxx', Pa
        syym.setZero();        // Marker SIGMAyy', Pa
        sxym.setZero();        // Marker SIGMAxy', Pa
        gsm.setZero();         // Standard shear modulus of bulk, Pa
        gm.setZero();          // Shear modulus of bulk, Pa
        vx0m.setZero();        // Marker horizontal velocity
        vy0m.setZero();        // Marker vertical velocity
        ptfm.setZero();        // Pt - Pf, Pa
        amursfm.setZero();     // RSF a / mu parameter
        
        int m = 0;
        for (int jm = 0; jm < Nx_markers; jm++) {
            for (int im = 0; im < Ny_markers; im++, m++) { // Update marker index each time
                // Define randomized regular coordinates
                xm(m) = xbeg + jm * dxms + (rand() % 1) * dxms;
                ym(m) = ybeg + im * dyms + (rand() % 1) * dyms;
                //Matrix
                t_marker(m) = 1;
                rhom(m) = 3000;
                etasm(m) = 1e22;
                gsm(m) = shearmod;
                cohescm(m) = cohes;
                cohestm(m) = cohes;
                frictcm(m) = friction;
                dilatcm(m) = dilatation;
                fricttm(m) = tensile;
                porm(m) = .01 * (1 + .0 * (rand() % 1 - .5));
                kkkm(m) = 5e-16; // * (dy / faultwidth)^2;
                rhofm(m) = 1000;
                etafm(m) = 1e-3;
                etam(m) = etasm(m) * exp(-alpha * porm(m));
                etabm(m) = etam(m) / porm(m); // / (1 - porm(m));
                gm(m) = gsm(m); // * (1 - porm(m));
                
                // Air, wedge, slab
                if (ym(m) < upper_block || ym(m) > lower_block) {
                    t_marker(m) = -1;
                    etam(m) = 1e23;
                    etasm(m) = 1e23;
                    rhom(m) = 3000;
                    kkkm(m) = 5e-16; // * (dy / faultwidth)^2;  // does not change the value of that variable
                    cohescm(m) = cohes * 1e3;
                    cohestm(m) = cohes * 1e3;
                    frictcm(m) = friction;
                    dilatcm(m) = dilatation;
                    fricttm(m) = tensile;
                }
            }
        }
        L.setZero();
        R.setZero();
    }
    
    // /////////////////////////////////////////////////////////////////////////////////////// 
    // actual computations start here
    // /////////////////////////////////////////////////////////////////////////////////////// 

    for (; timestep <= num_timesteps; timestep++) {

        // Interpolate ETA, RHO to nodal points
        // Basic nodes
        RHOSUM.setZero();
        ETASUM.setZero();
        KKKSUM.setZero();
        TTTSUM.setZero();
        SXYSUM.setZero();
        GGGSUM.setZero();
        ETA.setZero();
        ETA0SUM.setZero();
        COHCSUM.setZero();
        FRICSUM.setZero();
        DILCSUM.setZero();
        COHTSUM.setZero();
        FRITSUM.setZero();
        WTSUM.setZero();

        //LDZ
        OM0SUM.setZero();  // Old state parameter
        OMSUM.setZero();   // State parameter
        ARSFSUM.setZero(); // a - parameter of RSF
        BRSFSUM.setZero(); // b - parameter of RSF
        LRSFSUM.setZero(); // L - parameter of RSF

        // Pressure nodes
        ETAPSUM.setZero();
        ETAP0SUM.setZero();
        ETAB0SUM.setZero();
        PORSUM.setZero();
        SXXSUM.setZero();
        SYYSUM.setZero();
        GGGPSUM.setZero();
        WTPSUM.setZero();
        // Vx nodes
        RHOXSUM.setZero();
        RHOFXSUM.setZero();
        ETADXSUM.setZero();
        PORXSUM.setZero();
        VX0SUM.setZero();
        WTXSUM.setZero();
        // Vy nodes
        RHOYSUM.setZero();
        RHOFYSUM.setZero();
        ETADYSUM.setZero();
        PORYSUM.setZero();
        VY0SUM.setZero();
        WTYSUM.setZero();

        // Cycle on markers
        #pragma omp parallel for // about 3-4x faster with n = 4
        for (int m = 0; m < marknum; m++) {
            double cohescmm, cohestmm, frictcmm, dilatcmm, fricttmm, etasmm0, etamm0, etamm, rhomm, etadm;

            // Marker properties
            double kkkmm = kkkm(m) * pow(porm(m) / POR0, 3);

            // Viscosity of porous matrix
            if (t_marker(m) != 0) {
                cohescmm = cohescm(m) * (1 - porm(m)); // * exp(-alpha * porm(m));
                cohestmm = cohestm(m) * (1 - porm(m)); // * exp(-alpha * porm(m));
                frictcmm = frictcm(m);
                dilatcmm = dilatcm(m);
                fricttmm = fricttm(m);
                etasmm0 = etasm(m);
                etamm0 = etasm(m) * exp(-alpha * porm(m));
                etamm = etam(m);
                // total density
                rhomm = rhom(m) * (1 - porm(m)) + rhofm(m) * porm(m);
            } else {
                cohescmm = cohescm(m);
                cohestmm = cohestm(m);
                frictcmm = frictcm(m);
                dilatcmm = 0;
                fricttmm = fricttm(m);
                etasmm0 = etasm(m);
                etamm0 = etasm(m);
                etamm = etam(m);
                // total density
                rhomm = rhom(m);
            }

            // Darsi "viscosity"
            etadm = etafm(m) / kkkmm;
            
            // Interpolate to basic nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            int j = check_bounds(fix(xm(m) / dx), Nx);
            int i = check_bounds(fix(ym(m) / dy), Ny);
            
            double dxm = (xm(m) - x(j)) / dx;
            double dym = (ym(m) - y(i)) / dy;
            
            // merged similar computations of matlab version
            // The computations on the 4 elements are now done on a 2x2-sub-block of the matrix
            // Similarly wtm is a 2x2 matrix that stores the corresponding value
            // Sub-block += variable * wtm

            Eigen::Matrix2d wtm;
            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;

            ETASUM.block(i, j, 2, 2) += etamm * wtm;
            RHOSUM.block(i, j, 2, 2) += rhomm * wtm;
            KKKSUM.block(i, j, 2, 2) += kkkmm * wtm;
            TTTSUM.block(i, j, 2, 2) += t_marker(m) * wtm;
            SXYSUM.block(i, j, 2, 2) += sxym(m) * wtm;
            GGGSUM.block(i, j, 2, 2) += 1. / gm(m) * wtm;
            ETA0SUM.block(i, j, 2, 2) += etamm0 * wtm;
            COHTSUM.block(i, j, 2, 2) += cohestmm * wtm;
            FRITSUM.block(i, j, 2, 2) += fricttmm * wtm;
            COHCSUM.block(i, j, 2, 2) += cohescmm * wtm;
            FRICSUM.block(i, j, 2, 2) += frictcmm * wtm;
            DILCSUM.block(i, j, 2, 2) += dilatcmm * wtm;

            WTSUM.block(i, j, 2, 2) += wtm;
            
            // Interpolate to pressure nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix((xm(m) + dx / 2.) / dx), Nx);
            i = check_bounds(fix((ym(m) + dy / 2.) / dy), Ny);
            
            dxm = (xm(m) - xp(j)) / dx;
            dym = (ym(m) - yp(i)) / dy;

            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;

            ETAPSUM.block(i, j, 2, 2) += etamm * wtm;
            PORSUM.block(i, j, 2, 2) += porm(m) * wtm;
            SXXSUM.block(i, j, 2, 2) += sxxm(m) * wtm;
            SYYSUM.block(i, j, 2, 2) += syym(m) * wtm;
            GGGPSUM.block(i, j, 2, 2) += 1. / gm(m) * wtm;
            ETAP0SUM.block(i, j, 2, 2) += etamm0 * wtm;
            ETAB0SUM.block(i, j, 2, 2) += etasmm0 * wtm;

            WTPSUM.block(i, j, 2, 2) += wtm;

            // Interpolate to Vx nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix((xm(m)) / dx), Nx);
            i = check_bounds(fix((ym(m) + dy / 2.) / dy), Ny);
            
            dxm = (xm(m) - xvx(j)) / dx;
            dym = (ym(m) - yvx(i)) / dy;
            
            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;

            RHOXSUM.block(i, j, 2, 2) += rhomm * wtm;
            ETADXSUM.block(i, j, 2, 2) += 1. / etadm * wtm;
            RHOFXSUM.block(i, j, 2, 2) += rhofm(m) * wtm;
            PORXSUM.block(i, j, 2, 2) += porm(m) * wtm;
            
            WTXSUM.block(i, j, 2, 2) += wtm;
            
            // Interpolate to Vy nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix((xm(m) + dx / 2.) / dx), Nx);
            i = check_bounds(fix((ym(m)) / dy), Ny);
        
            dxm = (xm(m) - xvy(j)) / dx;
            dym = (ym(m) - yvy(i)) / dy;

            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;
            
            RHOYSUM.block(i, j, 2, 2) += rhomm * wtm;
            ETADYSUM.block(i, j, 2, 2) += 1. / etadm * wtm;
            RHOFYSUM.block(i, j, 2, 2) += rhofm(m) * wtm;
            PORYSUM.block(i, j, 2, 2) += porm(m) * wtm;

            WTYSUM.block(i, j, 2, 2) += wtm;
        } // ends loop through markers

        // Computing ETA and RHO
        for (int i = 0; i < Ny; i++) {
            // Basic nodes
            for (int j = 0; j < Nx; j++) {
                double wtsum = WTSUM(i, j);
                if (wtsum > 0) {
                    RHO(i, j) = RHOSUM(i, j) / wtsum;
                    KKK(i, j) = KKKSUM(i, j) / wtsum;
                    TTT(i, j) = TTTSUM(i, j) / wtsum;
                    GGG(i, j) = 1. / (GGGSUM(i, j) / wtsum);
                    ETA0(i, j) = ETA0SUM(i, j) / wtsum;
                    COHT(i, j) = COHTSUM(i, j) / wtsum;
                    FRIT(i, j) = FRITSUM(i, j) / wtsum;
                    COHC(i, j) = COHCSUM(i, j) / wtsum;
                    FRIC(i, j) = FRICSUM(i, j) / wtsum;
                    DILC(i, j) = DILCSUM(i, j) / wtsum;
                }
            }
            // Vy nodes
            for (int j = 0; j < Nx1; j++) {
                double wtysum = WTYSUM(i, j);
                if (wtysum > 0) {
                    RHOY(i, j) = RHOYSUM(i, j) / wtysum;
                    ETADY(i, j) = 1. / (ETADYSUM(i, j) / wtysum);
                    RHOFY(i, j) = RHOFYSUM(i, j) / wtysum;
                    PORY(i, j) = PORYSUM(i, j) / wtysum;
                }
            }
        }
        
        for (int i = 0; i < Ny1; i++) {
            // Vx nodes
            for (int j = 0; j < Nx; j++) {
                double wtxsum = WTXSUM(i, j);
                if (wtxsum > 0) {
                    RHOX(i, j) = RHOXSUM(i, j) / wtxsum;
                    ETADX(i, j) = 1. / (ETADXSUM(i, j) / wtxsum);
                    RHOFX(i, j) = RHOFXSUM(i, j) / wtxsum;
                    PORX(i, j) = PORXSUM(i, j) / wtxsum;
                }
            }
            //Pressure nodes
            for (int j = 0; j < Nx1; j++) {
                double wtpsum = WTPSUM(i, j);
                if (wtpsum > 0) {
                    ETAP(i, j) = ETAPSUM(i, j) / wtpsum;
                    POR(i, j) = PORSUM(i, j) / wtpsum;
                    ETAB0(i, j) = ETAB0SUM(i, j) / wtpsum / POR(i, j);
                    ETAP0(i, j) = ETAP0SUM(i, j) / wtpsum;
                    GGGP(i, j) = 1. / (GGGPSUM(i, j) / wtpsum);
                }
            }
        }
        
        // Save viscosity
        if (timestep == 1) {
            ETA1 = ETA0;
            ETA = ETA0;
            ETA50 = ETA0;
        } else {
            ETA1 = ETA00;
            ETA = ETA00;
            ETA50 = ETA00;
        }
        YNY00 = YNY0;
        
        OM = OM0;
        
        // Multiple solving of equations
        if (!yndtdecrease) {
            dt = std::max(std::min(dt * dtkoefup, dtelastic0), dtmin);
        } else {
            dt = std::max(std::min(dt, dtelastic0), dtmin);
        }
        
        yndtdecrease = false;
        dt00 = dt;
        DSYLSQ.setZero();
        ynlast = 0;
        
        for (iterstep = 0; iterstep < niterglobal; iterstep++) {
            
            // Limiting viscosity
            etamincur = shearmod * dt * 1e-4;
            
            // External P - nodes: symmetry
            pt.col(0) = pt.col(1);
            pt.col(Nx) = pt.col(Nx - 1);
            pt.row(0) = pt.row(1);
            pt.row(Ny) = pt.row(Ny - 1);
            pf.col(0) = pf.col(1);
            pf.col(Nx) = pf.col(Nx - 1);
            pf.row(0) = pf.row(1);
            pf.row(Ny) = pf.row(Ny - 1); 
            
            // Basic nodes
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    if (ETA(i, j) < etamincur) {
                        ETA(i, j) = etamincur;
                    }
                    // Compute plastic strain rate
                    if (ETA(i, j) < ETA0(i, j)) {
                        double SXX_tot_4 = (SXX(i, j) + SXX(i + 1, j) + SXX(i, j + 1) + SXX(i + 1, j + 1)) / 4.;
                        double SYY_tot_4 = (SYY(i, j) + SYY(i + 1, j) + SYY(i, j + 1) + SYY(i + 1, j + 1)) / 4.;
                        SIIB(i, j) = sqrt(pow(SXY(i, j), 2) + .5 * pow(SXX_tot_4, 2) + .5 * pow(SYY_tot_4, 2) + .5 * pow(-SXX_tot_4 - SYY_tot_4, 2));
                        EIIB(i, j) = dy / faultwidth * (SIIB(i, j) / 2. / ETA(i, j) - SIIB(i, j) / 2. / ETA0(i, j));
                        IETAPLB(i, j) = (1. / ETA(i, j) - 1. / ETA0(i, j));
                    } else {
                        EIIB(i, j) = 0;
                        IETAPLB(i, j) = 0;
                    }
                }
            }
            
            // Computing viscosity and dilatation in pressure nodes
            for (int i = 1; i < Ny; i++) {
                for (int j = 1; j < Nx; j++) {
                    // Compute viscoplastic viscosity
                    IETAPL = (IETAPLB(i - 1, j - 1) + IETAPLB(i, j - 1) + IETAPLB(i - 1, j) + IETAPLB(i, j)) / 4.;
                    if (YNY0(i - 1, j - 1) > 0 || YNY0(i, j - 1) > 0 || YNY0(i - 1, j) > 0 || YNY0(i, j) > 0) {
                        ETAP(i, j) = 1. / (1. / ETAP0(i, j) + IETAPL);
                        ETAB(i, j) = 1. / (1. / ETAB0(i, j) + dy / faultwidth * IETAPL * POR(i, j));
                    } else {
                        ETAP(i, j) = ETAP0(i, j);
                        ETAB(i, j) = ETAB0(i, j);
                    }
                    // Check viscosity
                    if (ETAP(i, j) < etamincur) {
                        ETAP(i, j) = etamincur;
                    }
                    if (ETAB(i, j) * POR(i, j) < etamincur) {
                        ETAB(i, j) = etamincur / POR(i, j);
                    }
                    // Pores compressibility
                    GGGB(i, j) = GGGP(i, j) / POR(i, j);
                    // Dilation
                    // Zhao and Cai, International Journal of Rock Mechanics & Mining Sciences 47 (2010) 368–384
                    // Weak sandstone parameters 
                    
                    ss3 = std::min(std::max((pt(i, j) - pf(i, j)) * 1e-6, 0.), 100.); // SIGMA3, MPa
                    aa = aa1 + aa2 * exp(-ss3 / aa3);
                    bb = bb1 + bb2 * exp(-ss3 / bb3);
                    cc = cc1 + cc2 / 100. * pow(ss3, cc3);
                    dilij = sin(aa * bb * (exp(-bb * gammapij) - exp(-cc * gammapij)) / (cc - bb) / 180. * pi);
                    DILP(i, j) = 0; //2 * (dili1j1 * EIIB(i - 1, j - 1) + dilij1 * EIIB(i, j - 1) + dili1j * EIIB(i - 1, j) + dilij * EIIB(i, j)) / 4;
                }
            }
            
            if (dt > 2e5) {
                pfscale = ETADX.block(1, 0, Ny - 1, Nx).minCoeff() * dx * 1e17 / pow(dt, 2);
            } else {
                pfscale = ETADX.block(1, 0, Ny - 1, Nx).minCoeff() * dx;
            }
            
            ptscale = pfscale;
            
            // correct ->
            // 5)Composing global matrixes L(), R()
            for (int j = 0; j < Nx1; j++) {
                for (int i = 0; i < Ny1; i++) {
                    // Computing global indexes for vx, vy, p
                    kp = (j * Ny1 + i) * 6;
                    kx = kp + 1;
                    ky = kp + 2;
                    kpf = kp + 3;
                    kxf = kp + 4;
                    kyf = kp + 5;
                    
                    // 5a) Composing equation for vxs
                    if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (j == Nx) {
                            L.coeffRef(kx, kx) = 1;
                            R(kx) = 0;
                        }
                        
                        // Upper boundary
                        // prescribed velocity
                        if (i == 0 && j < Nx) {
                            L.coeffRef(kx, kx) = 1;
                            L.coeffRef(kx, kx + 6) = 1;
                            R(kx) = 2 * bcupper;
                        }
                        
                        // Lower boundary
                        // prescribed velocity
                        if (i == Ny && j < Nx) {
                            L.coeffRef(kx, kx) = 1;
                            L.coeffRef(kx, kx - 6) = 1;
                            R(kx) = 2 * bclower;
                        }
                        
                        // Left boundary:
                        if (j == 0 && i > 0 && i < Ny) {
                            L.coeffRef(kx, kx) = 1;
                            L.coeffRef(kx, kx + 6 * Ny1) = -1;
                            R(kx) = 0;
                        }
                        
                        // Right boundary
                        if (j == Nx - 1 && i > 0 && i < Ny) {
                            L.coeffRef(kx, kx) = 1;
                            L.coeffRef(kx, kx - 6 * Ny1) = -1;
                            R(kx) = 0;
                        }
                        
                    } else {
                        // Total X - Stokes: dSIGMAxxt' / dx + dSIGMAxyt' / dy - dPt / dx = -RHOt * gx
                        // SIGMAijt = 2 * ETA * EPSILONijs * K + SIGMAijt0 * (1 - K)
                        //             vxs2
                        //        vys1  |    vys3
                        //              |
                        //  vxs1 -- Pt1 -- vxs3 -- Pt2 -- vxs5
                        //              |
                        //        vys2  |    vys4
                        //             vxs4
                        // Viscosity
                        ETAXY1 = ETA(i - 1, j);
                        ETAXY2 = ETA(i, j);
                        ETAXX1 = ETAP(i, j);
                        ETAXX2 = ETAP(i, j + 1);
                        // Shear modulus
                        GXY1 = GGG(i - 1, j);
                        GXY2 = GGG(i, j);
                        GXX1 = GGGP(i, j);
                        GXX2 = GGGP(i, j + 1);
                        // Viscoelasticity factor
                        KXY1 = dt * GXY1 / (dt * GXY1 + ETAXY1);
                        KXY2 = dt * GXY2 / (dt * GXY2 + ETAXY2);
                        KXX1 = dt * GXX1 / (dt * GXX1 + ETAXX1);
                        KXX2 = dt * GXX2 / (dt * GXX2 + ETAXX2);
                        // Numerical viscosity
                        ETAXY1 *= KXY1;
                        ETAXY2 *= KXY2;
                        ETAXX1 *= KXX1;
                        ETAXX2 *= KXX2;
                        // Numerical stresses
                        SXY1 = SXY0(i - 1, j) * (1 - KXY1);
                        SXY2 = SXY0(i, j) * (1 - KXY2);
                        SXX1 = SXX0(i, j) * (1 - KXX1);
                        SXX2 = SXX0(i, j + 1) * (1 - KXX2);
                        // Density derivatives
                        dRHOdx = (RHOX(i, j + 1) - RHOX(i, j - 1)) / (2. * dx);
                        dRHOdy = (RHO(i, j) - RHO(i - 1, j)) / dy;
                        // Left part
                        double dx2 = pow(dx, 2), dy2 = pow(dy, 2);
                        L.coeffRef(kx, kx) = -4. / 3. * (ETAXX1 + ETAXX2) / dx2 - (ETAXY1 + ETAXY2) / dy2 - gx * dt * dRHOdx - ascale * RHOX(i, j) / dt; //vxs3
                        L.coeffRef(kx, kx - Ny1 * 6) = 4. / 3. * ETAXX1 / dx2; //vxs1
                        L.coeffRef(kx, kx + Ny1 * 6) = 4. / 3. * ETAXX2 / dx2; //vxs5
                        L.coeffRef(kx, kx - 6) = ETAXY1 / dy2; //vxs2
                        L.coeffRef(kx, kx + 6) = ETAXY2 / dy2; //vxs4
                        double gx_dt_dRHOdy = gx * dt * dRHOdy / 4.;
                        double dx_dy = dx * dy;
                        L.coeffRef(kx, ky - 6) = ETAXY1 / dx_dy - 2. / 3. * ETAXX1 / dx_dy - gx_dt_dRHOdy; //vys1
                        L.coeffRef(kx, ky) = -ETAXY2 / dx_dy + 2. / 3. * ETAXX1 / dx_dy - gx_dt_dRHOdy; //vys2
                        L.coeffRef(kx, ky - 6 + Ny1 * 6) = -ETAXY1 / dx_dy + 2. / 3. * ETAXX2 / dx_dy - gx_dt_dRHOdy; //vys3
                        L.coeffRef(kx, ky + Ny1 * 6) = ETAXY2 / dx_dy - 2. / 3. * ETAXX2 / dx_dy - gx_dt_dRHOdy; //vys4
                        L.coeffRef(kx, kp) = ptscale / dx; //Pt1'
                        L.coeffRef(kx, kp + Ny1 * 6) = -ptscale / dx; //Pt2'
                        // Right part
                        R(kx) = -RHOX(i, j) * (ascale * VX0(i, j) / dt + gx) - (SXX2 - SXX1) / dx - (SXY2 - SXY1) / dy;
                    }
                    
                    // 5b) Composing equation for vys
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vys = 0
                        if (i == Ny) {
                            L.coeffRef(ky, ky) = 1;
                            R(ky) = 0;
                        }
                        
                        // Left boundary
                        // Free Slip
                        if (j == 0) {
                            L.coeffRef(ky, ky) = 1;
                            L.coeffRef(ky, ky + Ny1 * 6) = 1;
                            R(ky) = 0;
                        }
                        
                        // Right boundary
                        // Free Slip
                        if (j == Nx) {
                            L.coeffRef(ky, ky) = 1;
                            L.coeffRef(ky, ky - Ny1 * 6) = 1;
                            R(ky) = 0;
                        }
                        
                        // Upper boundary: no penetration
                        if (i == 0 && j > 0 && j < Nx) {
                            L.coeffRef(ky, ky) = 1;
                            R(ky) = 0;
                        }
                        
                        // Lower boundary: no penetration
                        if (i == Ny - 1 && j > 0 && j < Nx) {
                            L.coeffRef(ky, ky) = 1;
                            R(ky) = 0;
                        }
                        
                    } else {
                        // Total Y - Stokes: dSIGMAyxt' / dx + dSIGMAyyt' / dy - dPt / dy = -RHOt * gy
                        // y - Stokes equation: dSIGMA'yx / dx + dSIGMA'yy / dy - dP / dy = -RHO * gy
                        //
                        //               vys2
                        //                |
                        //         vxs1  Pt1  vxs3
                        //                |
                        //   vys1 --------- vys3 -------- vys5
                        //                |
                        //         vxs2  Pt2  vxs4
                        //                |
                        //               vys4
                        // Viscosity
                        ETAXY1 = ETA(i, j - 1);
                        ETAXY2 = ETA(i, j);
                        ETAYY1 = ETAP(i, j);
                        ETAYY2 = ETAP(i + 1, j);
                        // Shear modulus
                        GXY1 = GGG(i, j - 1);
                        GXY2 = GGG(i, j);
                        GYY1 = GGGP(i, j);
                        GYY2 = GGGP(i + 1, j);
                        // Viscoelasticity factor
                        KXY1 = dt * GXY1 / (dt * GXY1 + ETAXY1);
                        KXY2 = dt * GXY2 / (dt * GXY2 + ETAXY2);
                        KYY1 = dt * GYY1 / (dt * GYY1 + ETAYY1);
                        KYY2 = dt * GYY2 / (dt * GYY2 + ETAYY2);
                        // Numerical viscosity
                        ETAXY1 *= KXY1;
                        ETAXY2 *= KXY2;
                        ETAYY1 *= KYY1;
                        ETAYY2 *= KYY2;
                        // Numerical stresses
                        SXY1 = SXY0(i, j - 1) * (1 - KXY1);
                        SXY2 = SXY0(i, j) * (1 - KXY2);
                        SYY1 = SYY0(i, j) * (1 - KYY1);
                        SYY2 = SYY0(i + 1, j) * (1 - KYY2);
                        // Density derivatives
                        dRHOdy = (RHOY(i + 1, j) - RHOY(i - 1, j)) / 2. / dy;
                        dRHOdx = (RHO(i, j) - RHO(i, j - 1)) / dx;
                        // Left part
                        double dx2 = pow(dx, 2), dy2 = pow(dy, 2);
                        L.coeffRef(ky, ky) = -4. / 3. * (ETAYY1 + ETAYY2) / dy2 - (ETAXY1 + ETAXY2) / dx2 - gy * dt * dRHOdy - ascale * RHOY(i, j) / dt; //vys3
                        L.coeffRef(ky, ky - Ny1 * 6) = ETAXY1 / dx2; //vys1
                        L.coeffRef(ky, ky + Ny1 * 6) = ETAXY2 / dx2; //vys5
                        L.coeffRef(ky, ky - 6) = 4. / 3. * ETAYY1 / dy2; //vys2
                        L.coeffRef(ky, ky + 6) = 4. / 3. * ETAYY2 / dy2; //vys4
                        double gy_dt_dRHOdx = gy * dt * dRHOdx / 4.;
                        double dx_dy = dx * dy;
                        L.coeffRef(ky, kx - Ny1 * 6) = ETAXY1 / dx_dy - 2. / 3. * ETAYY1 / dx_dy - gy_dt_dRHOdx; //vxs1
                        L.coeffRef(ky, kx + 6 - Ny1 * 6) = -ETAXY1 / dx_dy + 2. / 3. * ETAYY2 / dx_dy - gy_dt_dRHOdx; //vxs2
                        L.coeffRef(ky, kx) = -ETAXY2 / dx_dy + 2. / 3. * ETAYY1 / dx_dy - gy_dt_dRHOdx; //vxs3
                        L.coeffRef(ky, kx + 6) = ETAXY2 / dx_dy - 2. / 3. * ETAYY2 / dx_dy - gy_dt_dRHOdx; //vxs4
                        L.coeffRef(ky, kp) = ptscale / dy; //Pt1'
                        L.coeffRef(ky, kp + 6) = -ptscale / dy; //Pt2'
                        // Right part
                        R(ky) = -RHOY(i, j) * (ascale * VY0(i, j) / dt + gy) - (SYY2 - SYY1) / dy - (SXY2 - SXY1) / dx;
                    }
                    
                    // 5c) Composing equation for Pt
                    if (i == 0 || j == 0 || i == Ny || j == Nx) { // || (i == 2 && j == 2))
                        // BC equation: 1 * Pt = 0
                        L.coeffRef(kp, kp) = 1;
                        R(kp) = 0;
                    } else {
                        // Solid Continuity: dVxs / dx + dVys / dy + (Pt - Pf) / ETAbulk = 0
                        //              vys1
                        //               |
                        //        vxs1 -- Pt, Pf -- vxs2
                        //               |
                        //              vys2
                        // Drained compressibility
                        BETADRAINED = (1. / GGGB(i, j) + BETASOLID) / (1 - POR(i, j));
                        // Biott - Willis koefficient
                        KBW = 1 - BETASOLID / BETADRAINED;
                        // Left part
                        L.coeffRef(kp, kx - Ny1 * 6) = -1. / dx; //vxs1
                        L.coeffRef(kp, kx) = 1. / dx; //vxs2
                        L.coeffRef(kp, ky - 6) = -1. / dy; //vys1
                        L.coeffRef(kp, ky) = 1. / dy; //vys2
                        L.coeffRef(kp, kp) = ptscale * (1. / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED / dt); //Pt
                        L.coeffRef(kp, kpf) = -pfscale * (1. / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / dt); //Pf
                        // Right part
                        R(kp) = BETADRAINED * (PT0(i, j) - KBW * PF0(i, j)) / dt + DILP(i, j);
                    }
                    
                    // 5d) Composing equation for vxD
                    if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (j == Nx) {
                            L.coeffRef(kxf, kxf) = 1;
                            R(kxf) = 0;
                        }
                        
                        // Upper boundary: symmetry
                        if (i == 0 && j < Nx) {
                            L.coeffRef(kxf, kxf) = 1;
                            L.coeffRef(kxf, kxf + 6) = -1;
                            R(kxf) = 0;
                        }
                        
                        // Lower boundary: symmetry
                        if (i == Ny && j < Nx) {
                            L.coeffRef(kxf, kxf) = 1;
                            L.coeffRef(kxf, kxf - 6) = -1;
                            R(kxf) = 0;
                        }
                        
                        // Left boundary
                        // no penetration
                        if (j == 0) {
                            L.coeffRef(kxf, kxf) = 1;
                            R(kxf) = 0; //bcvxfleft;
                        }
                        
                        // Right boundary
                        // no penetration
                        if (j == Nx - 1) {
                            L.coeffRef(kxf, kxf) = 1;
                            R(kxf) = 0;
                        }
                        
                    } else {
                        // Fluid X - Darsi:  - ETAfluid / K * VxD - dPf / dx = -RHOf * gx + RHOf * DVxs / Dt
                        //
                        //  Pf1 --- vxD, vxs --- Pf2
                        //
                        // Left part
                        L.coeffRef(kxf, kxf) = -ETADX(i, j) - RHOFX(i, j) / PORX(i, j) * ascale / dt; //vxD
                        L.coeffRef(kxf, kx) = -RHOFX(i, j) * ascale / dt; //vxs
                        L.coeffRef(kxf, kpf) = pfscale / dx; //Pf1'
                        L.coeffRef(kxf, kpf + Ny1 * 6) = -pfscale / dx; //Pf2'
                        // Right part
                        R(kxf) = -RHOFX(i, j) * (ascale * VXF0(i, j) / dt + gx);
                    }
                    
                    // 5e) Composing equation for vyD
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (i == Ny) {
                            L.coeffRef(kyf, kyf) = 1;
                            R(kyf) = 0;
                        }
                        
                        // Left boundary
                        // symmetry
                        if (j == 0 && i > 0 && i < Ny - 1) {
                            L.coeffRef(kyf, kyf) = 1;
                            L.coeffRef(kyf, kyf + Ny1 * 6) = -1;
                            R(kyf) = 0;
                        }
                        
                        // Right boundary
                        // symmetry
                        if (j == Nx && i > 0 && i < Ny - 1) {
                            L.coeffRef(kyf, kyf) = 1;
                            L.coeffRef(kyf, kyf - Ny1 * 6) = -1;
                            R(kyf) = 0;
                        }
                        
                        // Upper boundary: no penetration
                        if (i == 0) {
                            L.coeffRef(kyf, kyf) = 1;
                            R(kyf) = bcvyflower;
                        }
                        
                        // Lower boundary: no penetration
                        if (i == Ny - 1) {
                            L.coeffRef(kyf, kyf) = 1;
                            R(kyf) = bcvyflower;
                        }
                    } else {
                        // Fluid Y - Darsi:  - ETAfluid / K * VyD - dPf / dy = -RHOf * gy + RHOf * DVys / Dt
                        //
                        //   Pf1
                        //    |
                        //   vyD, vy
                        //    |
                        //   Pf2
                        //
                        // Left part
                        L.coeffRef(kyf, kyf) = -ETADY(i, j) - RHOFY(i, j) / PORY(i, j) * ascale / dt; //vyD
                        L.coeffRef(kyf, ky) = -RHOFY(i, j) * ascale / dt; //vys
                        L.coeffRef(kyf, kpf) = pfscale / dy; //Pf1'
                        L.coeffRef(kyf, kpf + 6) = -pfscale / dy; //Pf2'
                        // Right part
                        R(kyf) = -RHOFY(i, j) * (ascale * VYF0(i, j) / dt + gy);
                    }
                                        
                    // 5f) Composing equation for Pf
                    if (j == 0 || j == Nx || i <= 1 || i >= Ny - 1) { //same if clause but more compact
                        // BC equation: 1 * Pf = 0
                        L.coeffRef(kpf, kpf) = 1;
                        R(kpf) = 0;
                        
                        // Real BC
                        if (i == 1) {
                            L.coeffRef(kpf, kpf) = pfscale;
                            L.coeffRef(kpf, kp) = -ptscale;
                            R(kpf) = -PTFDIFF;
                        }
                        
                        // Real BC
                        if (i == Ny - 1) {
                            L.coeffRef(kpf, kpf) = pfscale;
                            L.coeffRef(kpf, kp) = -ptscale;
                            R(kpf) = -PTFDIFF;
                        }
                        
                    } else {
                        // Fluid Continuity: dVxD / dx + dVyD / dy - (Pt - Pf) / ETAbulk = 0
                        //              vyD1
                        //               |
                        //        vxD1 -- Pt, Pf -- vxD2
                        //               |
                        //              vyD2
                        // Compute elastic coefficients
                        // Drained compressibility
                        BETADRAINED = (1 / GGGB(i, j) + BETASOLID) / (1 - POR(i, j));
                        // Biott - Willis koefficient
                        KBW = 1 - BETASOLID / BETADRAINED;
                        // Skempton koefficient
                        KSK = (BETADRAINED - BETASOLID) / (BETADRAINED - BETASOLID + POR(i, j) * (BETAFLUID - BETASOLID));
                        // Left part
                        L.coeffRef(kpf, kxf - Ny1 * 6) = -1. / dx; //vxs1
                        L.coeffRef(kpf, kxf) = 1. / dx; //vxs2
                        L.coeffRef(kpf, kyf - 6) = -1. / dy; //vys1
                        L.coeffRef(kpf, kyf) = 1. / dy; //vys2
                        L.coeffRef(kpf, kp) = -ptscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / dt); //Pt
                        L.coeffRef(kpf, kpf) = pfscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / KSK / dt); //Pf
                        // Right part
                        R(kpf) = -BETADRAINED * KBW * (PT0(i, j) - 1. / KSK * PF0(i, j)) / dt - DILP(i, j);
                    }
                }
            }

            std::cout << "\nstart solving LSE" << std::endl; // testtesttest
            
            auto start = std::chrono::high_resolution_clock::now();

            // 6) Solving matrix
            Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
            L.makeCompressed();
            solver.analyzePattern(L);
            auto stop_fact_3 = std::chrono::high_resolution_clock::now();
            solver.factorize(L); // Computationally most costly line in the code

            auto stop_fact = std::chrono::high_resolution_clock::now();

            S = solver.solve(R);

            auto stop = std::chrono::high_resolution_clock::now();

            auto duration_fact_2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop_fact_3 - start);
            auto duration_fact_3 = std::chrono::duration_cast<std::chrono::milliseconds>(stop_fact - start);
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "factorisation        " << (duration_fact_3.count() - duration_fact_2.count()) / 1000. << " seconds" << std::endl;
            std::cout << "solving              " << (duration.count() - duration_fact_3.count()) / 1000. << " seconds" << std::endl;

            std::cout << "LSE success\n" << std::endl;
            
            // 7) Reload solution
            // pfavr = 0;
            // pcount = 0;
            for (int j = 0; j < Nx1; j++) {
                for (int i = 0; i < Ny1; i++) {
                    // Global indexes for vx, vy, P
                    kp = (j * Ny1 + i) * 6;
                    kx = kp + 1;
                    ky = kp + 2;
                    kpf = kp + 3;
                    kxf = kp + 4;
                    kyf = kp + 5;
                    // Reload solution
                    pt(i, j) = S(kp) * ptscale;
                    vxs(i, j) = S(kx);
                    vys(i, j) = S(ky);
                    pf(i, j) = S(kpf) * pfscale;
                    vxD(i, j) = S(kxf);
                    vyD(i, j) = S(kyf);
                }
            }

            Vmax = VSLIPB.maxCoeff();
            
            if (dt > 1e4 && Vmax < 1e-7) {
                avgpt = pt.sum() / (double)(pt.rows() * pt.cols()); //calculate average total pressure
                diffpt = (PCONF + PTFDIFF) - avgpt;
                pt += Eigen::MatrixXd::Constant(Ny1, Nx1, diffpt);
            }
            
            // Velocity change
            DVX0 = vxs - VX0;
            DVY0 = vys - VY0;
            
            // Define timestep
            dt0 = dt;
            bool yn = false;
            
            // Plastic iterations
            // Compute strain rate, stress and stress change
            ESP.setZero();
            EXY.setZero();
            SXY.setZero();
            DSXY.setZero();
            EXX.setZero();
            SXX.setZero();
            DSXX.setZero();
            EYY.setZero();
            SYY.setZero();
            DSYY.setZero();
            EII.setZero();
            EIIVP.setZero();
            SII.setZero();
            DIS.setZero();
            
            EL_DECOM.setZero();   // Elastic (de)compaction
            VIS_COMP.setZero();   // Viscous compaction
            
            // Process internal basic nodes
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    // ESP = .5 *(dVy / dx - dVx / dy), EXY, SXY, DSXY
                    ESP(i, j) = .5 * ((vys(i, j + 1) - vys(i, j)) / dx - (vxs(i + 1, j) - vxs(i, j)) / dy);
                    EXY(i, j) = .5 * ((vxs(i + 1, j) - vxs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dx);
                    KXY = dt * GGG(i, j) / (dt * GGG(i, j) + ETA(i, j));
                    SXY(i, j) = 2 * ETA(i, j) * EXY(i, j) * KXY + SXY0(i, j) * (1 - KXY);
                    DSXY(i, j) = SXY(i, j) - SXY0(i, j);
                }
            }

            // Process pressure cells
            for (int i = 1; i < Ny; i++) {
                for (int j = 1; j < Nx; j++) {
                    // EXX, SXX, DSXX
                    EXX(i, j) = (2 * (vxs(i, j) - vxs(i, j - 1)) / dx - (vys(i, j) - vys(i - 1, j)) / dy) / 3.;
                    EYY(i, j) = (2 * (vys(i, j) - vys(i - 1, j)) / dy - (vxs(i, j) - vxs(i, j - 1)) / dx) / 3.;
                    KXX = dt * GGGP(i, j) / (dt * GGGP(i, j) + ETAP(i, j));
                    SXX(i, j) = 2 * ETAP(i, j) * EXX(i, j) * KXX + SXX0(i, j) * (1 - KXX);
                    SYY(i, j) = 2 * ETAP(i, j) * EYY(i, j) * KXX + SYY0(i, j) * (1 - KXX);
                    DSXX(i, j) = SXX(i, j) - SXX0(i, j);
                    DSYY(i, j) = SYY(i, j) - SYY0(i, j);
                }
            }
            
            // External P - nodes: symmetry
            pt.col(0) = pt.col(1);
            pt.col(Nx) = pt.col(Nx - 1);
            pt.row(0) = pt.row(1);
            pt.row(Ny) = pt.row(Ny - 1);
            pf.col(0) = pf.col(1);
            pf.col(Nx) = pf.col(Nx - 1);
            pf.row(0) = pf.row(1);
            pf.row(Ny) = pf.row(Ny - 1);
            EXX.col(0) = EXX.col(1);
            EXX.col(Nx) = EXX.col(Nx - 1);
            EXX.row(0) = EXX.row(1);
            EXX.row(Ny) = EXX.row(Ny - 1);
            SXX.col(0) = SXX.col(1);
            SXX.col(Nx) = SXX.col(Nx - 1);
            SXX.row(0) = SXX.row(1);
            SXX.row(Ny) = SXX.row(Ny - 1);
            SXX0.col(0) = SXX0.col(1);
            SXX0.col(Nx) = SXX0.col(Nx - 1);
            SXX0.row(0) = SXX0.row(1);
            SXX0.row(Ny) = SXX0.row(Ny - 1);
            EYY.col(0) = EYY.col(1);
            EYY.col(Nx) = EYY.col(Nx - 1);
            EYY.row(0) = EYY.row(1);
            EYY.row(Ny) = EYY.row(Ny - 1);
            SYY.col(0) = SYY.col(1);
            SYY.col(Nx) = SYY.col(Nx - 1);
            SYY.row(0) = SYY.row(1);
            SYY.row(Ny) = SYY.row(Ny - 1);
            SYY0.col(0) = SYY0.col(1);
            SYY0.col(Nx) = SYY0.col(Nx - 1);
            SYY0.row(0) = SYY0.row(1);
            SYY0.row(Ny) = SYY0.row(Ny - 1);
            ETAP.col(0) = ETAP.col(1);
            ETAP.col(Nx) = ETAP.col(Nx - 1);
            ETAP.row(0) = ETAP.row(1);
            ETAP.row(Ny) = ETAP.row(Ny - 1);
            ETAB.col(0) = ETAB.col(1);
            ETAB.col(Nx) = ETAB.col(Nx - 1);
            ETAB.row(0) = ETAB.row(1);
            ETAB.row(Ny) = ETAB.row(Ny - 1);
            GGGP.col(0) = GGGP.col(1);
            GGGP.col(Nx) = GGGP.col(Nx - 1);
            GGGP.row(0) = GGGP.row(1);
            GGGP.row(Ny) = GGGP.row(Ny - 1);
            GGGB.col(0) = GGGB.col(1);
            GGGB.col(Nx) = GGGB.col(Nx - 1);
            GGGB.row(0) = GGGB.row(1);
            GGGB.row(Ny) = GGGB.row(Ny - 1);

            // Compute stress and strain rate invariants and dissipation
            // Process pressure cells
            for (int i = 1; i < Ny; i++) {
                for (int j = 1; j < Nx; j++) {
                    double sky_ij = SXY(i, j), sky_i1j = SXY(i - 1, j), sky_ij1 = SXY(i, j - 1), sky_i1j1 = SXY(i - 1, j - 1);
                    double eta_ij_2 = 2 * ETA(i, j), eta_i1j_2 = 2 * ETA(i - 1, j), eta_ij1_2 = 2 * ETA(i, j - 1), eta_i1j1_2 = 2 * ETA(i - 1, j - 1);
                    double sxx_ij = SXX(i, j), syy_ij = SYY(i, j);
                    double etap_2 = ETAP(i, j);

                    // EXY term is averaged from four surrounding basic nodes
                    EXY2 = (pow(EXY(i, j), 2) + pow(EXY(i - 1, j), 2) + pow(EXY(i, j - 1), 2) + pow(EXY(i - 1, j - 1), 2)) / 4.;
                    EII(i, j) = sqrt(0.5 * (pow(EXX(i, j), 2) + pow(EYY(i, j), 2)) + EXY2);
                    EXYVP2 = (pow(sky_ij / eta_ij_2, 2) + pow(sky_i1j / eta_i1j_2, 2) + pow(sky_ij1 / eta_ij1_2, 2) + pow(sky_i1j1 / eta_i1j1_2, 2)) / 4.;
                    EIIVP(i, j) = sqrt(0.5 * (pow(sxx_ij / etap_2, 2) + pow(syy_ij / etap_2, 2)) + EXYVP2);
                    // Second strain rate invariant SII
                    // SXY term is averaged from four surrounding basic nodes
                    SXY2 = (pow(sky_ij, 2) + pow(sky_i1j, 2) + pow(sky_ij1, 2) + pow(sky_i1j1, 2)) / 4.;
                    SII(i, j) = sqrt(0.5 * (pow(sxx_ij, 2) + pow(syy_ij, 2)) + SXY2);
                    
                    // Dissipation
                    DISXY = (pow(sky_ij, 2) / eta_ij_2 + pow(sky_i1j, 2) / eta_i1j_2 + pow(sky_ij1, 2) / eta_ij1_2 + pow(sky_i1j1, 2) / eta_i1j1_2) / 4.;
                    DIS(i, j) = pow(sxx_ij, 2) / etap_2 + pow(syy_ij, 2) / etap_2 + 2 * DISXY;
                }
            }
                        
            // Update viscosity for yielding
            AXY.setZero();
            // dt0 = dt;dt = dt * 1.1;
            ETA5 = ETA0;
            // Basic nodes
            DSY.setZero();
            YNY.setZero();
            SigmaY.setZero();
            SII_fault.setZero();
            int ynpl = 0;
            double ddd = 0;
            dtlapusta = 1e7;
            OM5 = OM;

            // Power law plasticity model Yi et al., 2018
            // Journal of Offshore Mechanics and Arctic Engineering
            double dtslip = 1e30;
            if (timestep > tyield) {
                for (int i = 0; i < Ny; i++) {
                    if (y(i)  >= upper_block && y(i)  <= lower_block) {
                        for (int j = 0; j < Nx; j++) {
                            // reducing matrix calls
                            double arsf_temp = ARSF(i, j);
                            double brsf_temp = BRSF(i, j);
                            double lrsf_temp = LRSF(i, j);
                            double fric_temp = FRIC(i, j);
                            double eta0_temp = ETA0(i, j);

                            double sxx_temp_4 = SXX(i, j) + SXX(i + 1, j) + SXX(i, j + 1) + SXX(i + 1, j + 1) / 4.;
                            double syy_temp_4 = SYY(i, j) + SYY(i + 1, j) + SYY(i, j + 1) + SYY(i + 1, j + 1) / 4.;

                            // SXX, pt are averaged from four surrounding pressure nodes
                            SIIB(i, j) = sqrt(pow(SXY(i, j), 2) + .5 * pow(sxx_temp_4, 2) + .5 * pow(syy_temp_4, 2) + .5 * pow(-sxx_temp_4 - syy_temp_4, 2)); // - can be changed to + as term is squared
                            ptB = (pt(i, j) + pt(i + 1, j) + pt(i, j + 1) + pt(i + 1, j + 1)) / 4.;
                            pfB = (pf(i, j) + pf(i + 1, j) + pf(i, j + 1) + pf(i + 1, j + 1)) / 4.;
                            // Computing "elastic" stress invariant
                            kfxy = ETA(i, j) / (GGG(i, j) * dt + ETA(i, j));
                            siiel = SIIB(i, j) / kfxy;
                            // Compute maximal stress invariant with viscous viscosity
                            kfxy0 = eta0_temp / (GGG(i, j) * dt + eta0_temp);
                            SIIB0 = siiel * kfxy0;
                            
                            // Compute old viscoplastic slip rate
                            // Compute PEFF
                            prB = (ptB - pfB);
                            
                            if (prB < 1e3) {
                                prB = 1e3;
                            }
                            // Compute old power law strain rate
                            SIIB1 = SIIB(i, j);
                            
                            //Compute slip velocity for current stress invariant and state
                            V = 2 * V0 * sinh(std::max(SIIB1, 0.) / arsf_temp / prB) * exp(-(brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                            
                            EIISLIP = V / dx / 2.;
                            
                            // Compute new ETAVP
                            ETAPL = SIIB1 / 2. / EIISLIP;
                            ETAVP = 1. / (1. / eta0_temp + 1. / ETAPL);
                            // Compute new stress invariant
                            kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                            SIIB2 = siiel * kfxy1;
                            DSIIB1 = SIIB2 - SIIB1;
                            
                            //Compute slip velocity for current stress invariant and state
                            V = 2 * V0 * sinh(std::max(SIIB2, 0.) / arsf_temp / prB) * exp(-(brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                            
                            EIISLIP = V / dx / 2.;
                            
                            // Compute new ETAVP
                            ETAPL = SIIB2 / 2. / EIISLIP;
                            ETAVP = 1. / (1. / eta0_temp + 1. / ETAPL);
                            // Compute new stress invariant
                            kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                            SIIB3 = siiel * kfxy1;
                            DSIIB2 = SIIB3 - SIIB2;
                            
                            if ((DSIIB1 >= 0 && DSIIB2 <= 0) || (DSIIB1 <= 0 && DSIIB2 >= 0)) {
                                double DSIIB = 1e9;
                                
                                while(abs(DSIIB) > 1e-3) {
                                    SIIB4 = (SIIB1 + SIIB2) / 2.;
                                    
                                    //Compute slip velocity for current stress invariant and state
                                    V = 2 * V0 * sinh(std::max((SIIB4), 0.) / arsf_temp / prB) * exp( - (brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                                    
                                    EIISLIP = V / dx / 2.;
                                    
                                    // Compute new ETAVP
                                    ETAPL = SIIB4 / 2. / EIISLIP;
                                    ETAVP = 1. / (1. / eta0_temp + 1. / ETAPL);
                                    // Compute new stress invariant
                                    kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                                    SIIB5 = siiel * kfxy1;
                                    DSIIB = SIIB5 - SIIB4;
                                    if ((DSIIB >= 0 && DSIIB1 >= 0) || (DSIIB <= 0 && DSIIB1 <= 0)) {
                                        SIIB1 = SIIB4;
                                    } else {
                                        SIIB2 = SIIB4;
                                    }
                                }
                            }
                            
                            if (V * dt / lrsf_temp > 1e-6) {
                                OM5(i, j) = log(V0 / V + (exp(OM0(i, j)) - V0 / V) * exp( - V * dt / lrsf_temp));
                            } else {
                                OM5(i, j) = log(exp(OM0(i, j)) * (1 - V * dt / lrsf_temp) + V0 * dt / lrsf_temp);
                            }
                            
                            
                            kfxy1 = ETAVP / (GGG(i, j) * dt + ETAVP);
                            SIGMA2 = siiel * kfxy1;
                            
                            // Compute yielding stress
                            syield = std::max(syieldmin, (ptB - pfB) * arsf_temp * asinh(V / 2. / V0 * exp((brsf_temp * OM5(i, j) + fric_temp) / arsf_temp)));
                            
                            // Compute visco - plastic viscosity
                            etapl = eta0_temp * syield / (eta0_temp * V + syield);
                            
                            //LDZ: save syield
                            SigmaY(i, j) = syield;
                            VSLIPB(i, j) = V;
                            SII_fault(i, j) = SIIB4;

                            // reduces calls on matrix
                            double g_temp = 2 * GGG(i, j); // reduces calls on matrix
                            
                            // "/ BETASOLID" -> Timestep criterion, Lapusta et al., 2000; Lapusta and Liu, 2009
                            double vi = (3. / BETASOLID - g_temp) / (6. / BETASOLID + g_temp);
                            double k = g_temp / (pi * (1 - vi) * dx);
                            double xi = .25 * pow(k * lrsf_temp / arsf_temp / prB - (brsf_temp - arsf_temp) / arsf_temp, 2) - k * lrsf_temp / arsf_temp / prB;
                            if (xi < 0) {
                                dTETAmax = std::min(1. - (brsf_temp - arsf_temp) * prB / (k * lrsf_temp), .2);
                            } else {
                                dTETAmax = std::min(arsf_temp * prB / (k * lrsf_temp - (brsf_temp - arsf_temp) * prB), .2);
                            }
                            dtlapusta = std::min(dtlapusta, dTETAmax * lrsf_temp / V);
                            
                            
                            double A = syield / siiel;
                            AXY(i, j) = A;
                            // Count old yelding nodes
                            bool ynn = false;
                            if (YNY0(i, j) > 0) {
                                ynn = true;
                                DSY(i, j) = SIIB(i, j) - syield;
                                ddd += pow(DSY(i, j), 2);
                                ynpl++;
                            }
                            // Update viscosity
                            if (A < 1) {
                                // New viscosity for the basic node
                                etapl = dt * GGG(i, j) * A / (1 - A);
                                if (etapl < eta0_temp) {
                                    // Update plastic nodes
                                    ETA5(i, j) = pow(etapl, 1 - etawt) * pow(ETA(i, j), etawt);
                                    YNY(i, j) = 1;
                                    // Count yelding nodes
                                    if (!ynn) {
                                        DSY(i, j) = SIIB(i, j) - syield;
                                        ddd += pow(DSY(i, j), 2);
                                        ynpl++;
                                    }
                                } else {
                                    ETA5(i, j) = eta0_temp;
                                }
                            } else {
                                ETA5(i, j) = eta0_temp;
                            }
                        }
                    }
                }
            }
            
            // Compute Error
            DSYLSQ(iterstep) = 0;
            if (ynpl > 0) {
                DSYLSQ(iterstep) = sqrt(ddd / ynpl);
            }
            if (ynpl == 0) {
                ETA = ETA0;
            }

            // connot calculate DSYLSQ(iterstep - 1) if iterstep = 0
            // need to know what to do when iterstep = 0
            double D_iter = DSYLSQ(iterstep);
            double D_iter_quot = 0;
            if (iterstep != 0) {
                D_iter_quot = D_iter / DSYLSQ(iterstep - 1);
            }
            
            
            // Adjust timestep
            double dtpl = dt;
            // if (ynlast >= dtstep && ynpl > 0 && DSYLSQ(iterstep) > errmax && iterstep < niterglobal)
            //     dtpl = dt / dtkoef
            //     yn = true;
            // end
            if (ynpl > 0 && iterstep < niterglobal && ynlast >= dtstep && (ynlast > ynlastmax || log10(D_iter_quot) >= 0 || log10(D_iter_quot) > log10(errmin / D_iter) / (ynlastmax - ynlast))) {
                dtpl = dt / dtkoef;
                yn = true;
            }
            
            double maxvxy0;
            // Define displacement timesteps
            if (iterstep > 0) {
                maxvxy0 = maxvxy;
            }
            double maxvxs = std::max(vxs.maxCoeff(), abs(vxs.minCoeff()));
            double maxvys = std::max(vys.maxCoeff(), abs(vys.minCoeff()));
            maxvxy = sqrt(pow(vxs.maxCoeff() - vxs.minCoeff(), 2) + pow(vys.maxCoeff() - vys.minCoeff(), 2));
            double stpmaxcur = stpmax1;
            dtx = dt;
            if (dt > dx * stpmaxcur / maxvxs) {
                dtx = dx / dtkoefv * stpmaxcur / maxvxs;
                yn = true;
            }
            dty = dt;
            if (dt > dy * stpmaxcur / maxvys) {
                dty = dy / dtkoefv * stpmaxcur / maxvys;
                yn = true;
            }
            maxvxs = 0;
            maxvys = 0;

            for (int i = 0; i < Ny1; i++) {
                    for (int j = 0; j < Nx1; j++) {
                    if (yvx(i) >= upper_block && yvx(i) <= lower_block) {
                        maxvxs = std::max(maxvxs, abs(vxs(i, j)));
                    }
                    if (yvy(i) >= upper_block && yvy(i) <= lower_block) {
                        maxvys = std::max(maxvys, abs(vys(i, j)));
                    }
                }
            }
            
            stpmaxcur = stpmax;
            if (dt > dx * stpmaxcur / maxvxs) {
                dtx = dx / dtkoefv * stpmaxcur / maxvxs;
                yn = true;
            }
            if (dt > dy * stpmaxcur / maxvys) {
                dty = dy / dtkoefv * stpmaxcur / maxvys;
                yn = true;
            }
            
            dtslip = 1e30;
            for (int i = 0; i < Ny; i++) {
                    for (int j = 0; j < Nx; j++) {
                    if (VSLIPB(i, j) > 0) {
                        dtslip = std::min(dtslip, dx * stpmax / VSLIPB(i, j));
                    }
                }
            }
            
            if (ynpl > 0 && dtslip < dt) {
                yn = true;
                dtslip = dtslip / dtkoefv;
            }
            
            
            // Chose minimal timestep
            if (yn && dt > dtmin) {
                double dtold = dt;
                dt = std::max(std::min(std::min(dtx, dty), std::min(std::min(dtpl, dtslip), dtlapusta)), dtmin);
                if (dt < dtold) {
                    ynlast = 0;
                }
            } else {
                yn = false;
            }

            // Exit iterations
            bool ynstop = false;
            double vratio;
            // Velocity change ratio
            if (iterstep > 0) {
                vratio = log10(maxvxy / maxvxy0);
            }
            if (!yn && (ynpl == 0 || (DSYLSQ(iterstep) < errmin && iterstep > 0 && abs(vratio) < vratiomax))) {
                ynstop = true;
            } else {
                // Recomputing ETA
                for (int i = 0; i < Ny; i++) {
                    for (int j = 0; j < Nx; j++) {
                        ETA(i, j) = std::max(std::min(ETA5(i, j), ETA0(i, j)), etamin);
                    }
                }
                // Save current viscosity
                ETA50 = ETA;
                YNY0 = YNY;
                OM = OM5;
            }
            
            // Exit iteration
            if (ynstop) {
                break;
            }

            ynlast++;
        }

        // /////////////////////////////////////////////////////////////////////////////////////// 
        // end of loop through global iterations
        // /////////////////////////////////////////////////////////////////////////////////////// 
        
        // Mark dt decrease
        if (dt00 > dt) {
            yndtdecrease = true;
        }
        
        // Save current viscosity
        ETA00 = ETA50;
        OM0 = OM;
        // Recheck displacement timestep
        dt = std::min(std::min(dtx, dty), std::min(dt, dtlapusta));
        
        // Compute strain rate, stress and stress change
        ESP.setZero();
        EXY.setZero();
        SXY.setZero();
        DSXY.setZero();
        EXX.setZero();
        SXX.setZero();
        DSXX.setZero();
        EYY.setZero();
        SYY.setZero();
        DSYY.setZero();
        EII.setZero();
        EIIVP.setZero();
        SII.setZero();
        DSII.setZero();
        DIS.setZero();

        // Process internal basic nodes
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                // ESP = .5 *(dVy / dx - dVx / dy), EXY, SXY, DSXY
                ESP(i, j) = .5 * ((vys(i, j + 1) - vys(i, j)) / dx - (vxs(i + 1, j) - vxs(i, j)) / dy);
                EXY(i, j) = .5 * ((vxs(i + 1, j) - vxs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dx);
                KXY = dt * GGG(i, j) / (dt * GGG(i, j) + ETA(i, j));
                SXY(i, j) = 2 * ETA(i, j) * EXY(i, j) * KXY + SXY0(i, j) * (1 - KXY);
                DSXY(i, j) = SXY(i, j) - SXY0(i, j);
            }
        }

        // Process pressure cells
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx; j++) {
                // EXX, SXX, DSXX
                EXX(i, j) = (2 * (vxs(i, j) - vxs(i, j - 1)) / dx - (vys(i, j) - vys(i - 1, j)) / dy) / 3.;
                EYY(i, j) = (2 * (vys(i, j) - vys(i - 1, j)) / dy - (vxs(i, j) - vxs(i, j - 1)) / dx) / 3.;
                KXX = dt * GGGP(i, j) / (dt * GGGP(i, j) + ETAP(i, j));
                SXX(i, j) = 2 * ETAP(i, j) * EXX(i, j) * KXX + SXX0(i, j) * (1 - KXX);
                SYY(i, j) = 2 * ETAP(i, j) * EYY(i, j) * KXX + SYY0(i, j) * (1 - KXX);
                DSXX(i, j) = SXX(i, j) - SXX0(i, j);
                DSYY(i, j) = SYY(i, j) - SYY0(i, j);
        
                // Compute stress and strain rate invariants and dissipation

                double sky_ij = SXY(i, j), sky_i1j = SXY(i - 1, j), sky_ij1 = SXY(i, j - 1), sky_i1j1 = SXY(i - 1, j - 1);
                double eta_ij_2 = 2 * ETA(i, j), eta_i1j_2 = 2 * ETA(i - 1, j), eta_ij1_2 = 2 * ETA(i, j - 1), eta_i1j1_2 = 2 * ETA(i - 1, j - 1);
            
                // EXY term is averaged from four surrounding basic nodes
                EXY2 = (pow(EXY(i, j), 2) + pow(EXY(i - 1, j), 2) + pow(EXY(i, j - 1), 2) + pow(EXY(i - 1, j - 1), 2)) / 4.;
                EII(i, j) = sqrt(pow(EXX(i, j), 2) + EXY2);
                EXYVP2 = (pow(sky_ij / eta_ij_2, 2) + pow(sky_i1j / eta_i1j_2, 2) + pow(sky_ij1 / eta_ij1_2, 2) + pow(sky_i1j1 / eta_i1j1_2, 2)) / 4.;
                EIIVP(i, j) = sqrt(.5 * (pow(SXX(i, j) / (2 * ETAP(i, j)), 2) + pow(SYY(i, j) / (2 * ETAP(i, j)), 2)) + EXYVP2);
                // Second strain rate invariant SII
                // SXY term is averaged from four surrounding basic nodes
                SXY2 = (pow(sky_ij, 2) + pow(sky_i1j, 2) + pow(sky_ij1, 2) + pow(sky_i1j1, 2)) / 4.;
                SII(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + pow(SYY(i, j), 2)) + SXY2);
                
                // Dissipation
                DISXY = (pow(sky_ij, 2) / eta_ij_2 + pow(sky_i1j, 2) /  eta_i1j_2 + pow(sky_ij1, 2) / eta_ij1_2 + pow(sky_i1j1, 2) / eta_i1j1_2) / 4.;
                DIS(i, j) = pow(SXX(i, j), 2) / (2 * ETAP(i, j)) + pow(SYY(i, j), 2) / (2 * ETAP(i, j)) + 2 * DISXY;
                
                
                if (i < Ny - 1) {
                    pt_ave(i, j) = (pt(i, j) + pt(i + 1, j)) / 2.;
                    pf_ave(i, j) = (pf(i, j) + pf(i + 1, j)) / 2.;
                    PT0_ave = (PT0(i, j) + PT0(i + 1, j)) / 2.;
                    PF0_ave = (PF0(i, j) + PF0(i + 1, j)) / 2.;
                } else {
                    pt_ave(i, j) = pt(i, j);
                    pf_ave(i, j) = pf(i, j);
                    PT0_ave = PT0(i, j);
                    PF0_ave = PF0(i, j);
                }
                
                // Compute elastic and viscous compaction
                VIS_COMP(i, j) = (pt_ave(i, j) - pf_ave(i, j)) / (ETAB(i, j) * (1 - POR(i, j)));
                // Drained compressibility
                BETADRAINED = (1 / GGGB(i, j) + BETASOLID) / (1 - POR(i, j));
                // Biott - Willis koefficient
                KBW = 1 - BETASOLID / BETADRAINED;
                EL_DECOM(i, j) = BETADRAINED * (pt_ave(i, j) - PT0_ave - KBW * pf_ave(i, j) + KBW * PF0_ave) / dt;
            }
        }

        // Runge-Kutta velocity, spin array
        vxm.setZero();
        vym.setZero();
        spm.setZero();

        // Move markers by nodal velocity field
        for (int m = 0; m < marknum; m++) {
            // Save marker position
            double xold = xm(m);
            double yold = ym(m);
            for (int rk = 0; rk < 4; rk++) {
                // vx - velocity interpolation
                // [i, j] -------- [i, j + 1]
                //   |                |
                //   |    o m         |
                //   |                |
                // [i + 1, j] ------- [i + 1, j + 1]
                // Indexes and distances
                int j = fix(xm(m) / dx);
                int i = fix((ym(m) + dy / 2.) / dy);
                
                j = check_bounds(j, Nx);
                i = check_bounds(i, Ny);

                //Distances
                dxm = (xm(m) - xvx(j)) / dx;
                dym = (ym(m) - yvx(i)) / dy;
                // Weights
                wtmij = (1 - dxm) * (1 - dym);
                wtmi1j = (1 - dxm) * dym;
                wtmij1 = dxm * (1 - dym);
                wtmi1j1 = dxm * dym;
                // Interpolation
                vxm(rk) = vxs(i, j) * wtmij + vxs(i + 1, j) * wtmi1j + vxs(i, j + 1) * wtmij1 + vxs(i + 1, j + 1) * wtmi1j1;
                
                // vy - velocity interpolation
                // [i, j] -------- [i, j + 1]
                //   |                |
                //   |    o m         |
                //   |                |
                // [i + 1, j] ------- [i + 1, j + 1]
                // Indexes and distances
                j = fix((xm(m) + dx / 2.) / dx);
                i = fix(ym(m) / dy);

                j = check_bounds(j, Nx);
                i = check_bounds(i, Ny);

                //Distances
                dxm = (xm(m) - xvy(j)) / dx;
                dym = (ym(m) - yvy(i)) / dy;
                // Weights
                wtmij = (1 - dxm) * (1 - dym);
                wtmi1j = (1 - dxm) * dym;
                wtmij1 = dxm * (1 - dym);
                wtmi1j1 = dxm * dym;
                // Interpolation
                vym(rk) = vys(i, j) * wtmij + vys(i + 1, j) * wtmi1j + vys(i, j + 1) * wtmij1 + vys(i + 1, j + 1) * wtmi1j1;
                
                
                // ESP = .5 *(dVy / dx - dVx / dy) interpolation
                // [i, j] -------- [i, j + 1]
                //   |                |
                //   |    o m         |
                //   |                |
                // [i + 1, j] ------- [i + 1, j + 1]
                // Indexes and distances
                j = fix((xm(m)) / dx);
                i = fix((ym(m)) / dy);

                j = check_bounds(j, Nx);
                i = check_bounds(i, Ny);
                
                //Distances
                dxm = (xm(m) - x(j)) / dx;
                dym = (ym(m) - y(i)) / dy;
                // Weights
                wtmij = (1 - dxm) * (1 - dym);
                wtmi1j = (1 - dxm) * dym;
                wtmij1 = dxm * (1 - dym);
                wtmi1j1 = dxm * dym;
                // Interpolation ESP = .5 *(dVy / dx - dVx / dy) for the marker
                spm(rk) = ESP(i, j) * wtmij + ESP(i + 1, j) * wtmi1j + ESP(i, j + 1) * wtmij1 + ESP(i + 1, j + 1) * wtmi1j1;
                
                // Moving between A, B, C, D points
                if (rk < 2) {
                    // Moving A -> B and A -> C
                    xm(m) = xold + vxm(rk) * dt / 2.;
                    ym(m) = yold + vym(rk) * dt / 2.;
                } else if (rk == 3) {
                    // Moving A -> D
                    xm(m) = xold + vxm(rk) * dt;
                    ym(m) = yold + vym(rk) * dt;
                }
            }
            // Compute effective velocity, rotation rate
            double vxeff = (vxm(0) + 2 * vxm(1) + 2 * vxm(2) + vxm(3)) / 6.;
            double vyeff = (vym(0) + 2 * vym(1) + 2 * vym(2) + vym(3)) / 6.;
            double speff = spm(0);
            
            // Rotate stress on marker according to its spin
            // Compute amount of rotation from spin rate:
            // Espin = .5 *(dvy / dx - dvx / dy) i.e. positive for clockwise rotation
            // (when x axis is directed rightward and y axis is directed downward)
            double dspeff = speff * dt;
            // Save old stresses
            double msxxold = sxxm(m);
            double msyyold = syym(m);
            double msxyold = sxym(m);
            sxym(m) = .5 * (msxxold - msyyold) * sin(2 * dspeff) + msxyold * cos(2 * dspeff);
            sxxm(m) = msxxold * pow(cos(dspeff), 2) + msyyold * pow(sin(dspeff), 2) - msxyold * sin(2 * dspeff);
            syym(m) = msxxold * pow(sin(dspeff), 2) + msyyold * pow(cos(dspeff), 2) + msxyold * sin(2 * dspeff);
            
            // Move markers
            xm(m) = xold + vxeff * dt;
            ym(m) = yold + vyeff * dt;
            
            // Recycling
            if (xm(m) < 0) {
                xm(m) = xm(m) + xsize;
            }
            if (xm(m) > xsize) {
                xm(m) = xm(m) - xsize;
            }
        }

        int t_minus1 = timestep - 1;
        
        // Update timesum
        timesum = timesum + dt;
        timesumcur(t_minus1) = timesum;
        dtcur(t_minus1) = dt;
        
        maxvxsmod(t_minus1) = -1e30;
        minvxsmod(t_minus1) = 1e30;
        maxvysmod(t_minus1) = -1e30;
        minvysmod(t_minus1) = 1e30;

        VX0 = vxs;
        VY0 = vys;
        
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                // Vx
                if (RHOX(i, j) > 2000 && i != 0) {
                    maxvxsmod(t_minus1) = std::max(maxvxsmod(t_minus1), vxs(i, j));
                    minvxsmod(t_minus1) = std::min(minvxsmod(t_minus1), vxs(i, j));
                }
                // Vy
                if (RHOY(i, j) > 2000 && j != 0) {
                    maxvysmod(t_minus1) = std::max(maxvysmod(t_minus1), vys(i, j));
                    minvysmod(t_minus1) = std::min(minvysmod(t_minus1), vys(i, j));
                }
            }
        }

        for (int i = 0; i < Ny1; i++) {
            for (int j = 0; j < Nx1; j++) {
                // Update VX0
                if (PORX(i, j) > 0) {
                    VXF0(i, j) = vxs(i, j) + vxD(i, j) / PORX(i, j);
                }
                // Update VY0
                if (PORY(i, j) > 0) {
                    VYF0(i, j) = vys(i, j) + vyD(i, j) / PORY(i, j);
                }
            }
        }

        // Update SXX0
        SXX0 = SXX;
        // Update SYY0
        SYY0 = SYY;
        // Update SXY0
        SXY0 = SXY;
        // Update PTF0
        PTF0 = (pt - pf);
        PT0 = pt;
        PF0 = pf;
        
        ///////////////////////////////////////////////////////////////////////////////////////// 
        // output
        /////////////////////////////////////////////////////////////////////////////////////// 
        
        std::cout << "====================================" << std::endl;
        std::cout << "total time:        " << timesum << " sec" << std::endl;
        std::cout << "time step:         " << dt << " sec" << std::endl;
        std::cout << "Vslip max:         " << Vmax << std::endl;
        std::cout << "iter - iterations:   " << iterstep << std::endl;
        std::cout << "global - iterations: " << ynlast << std::endl;
        
        if (seismic_cycles_data) {
            
            if (timesum == dt) {
                std::ofstream out_fault;
                out_fault.open("output/x_fault.txt");
                out_fault << x << std::endl;
                out_fault.close();
            }
            
            if (timesum_plus < timesum) {
                // ========== save slip rate
                std::ofstream out_Vslip("output/EVO_Vslip.txt", std::ios_base::app | std::ios_base::out);
                out_Vslip << timesum << "    " << dt << "    " << VSLIPB.row(line_fault) << std::endl; // might be .row() ???
                out_Vslip.close();
                
                // ========== save viscosity
                std::ofstream out_viscosity("output/EVO_viscosity.txt", std::ios_base::app | std::ios_base::out);
                out_viscosity << timesum << "    " << dt << "    " << ETA.row(line_fault) << std::endl; // might be .row() ???
                out_viscosity.close();
                
                // ========== save fluid 
                std::ofstream out_press_flu("output/EVO_press_flu.txt", std::ios_base::app | std::ios_base::out);
                out_press_flu << timesum << "    " << dt << "    " << pf.row(line_fault) << std::endl; // might be .row() ???
                out_press_flu.close();
                
                // ========== save effective pressure
                std::ofstream out_press_eff("output/EVO_press_eff.txt", std::ios_base::app | std::ios_base::out);
                Eigen::MatrixXd P_diff = pt - pf;
                out_press_eff << timesum << "    " << dt << "    " << P_diff.row(line_fault) << std::endl; // might be .row() ???
                out_press_eff.close();
                
                // ========== save SigmaY
                std::ofstream out_SigmaY("output/EVO_SigmaY.txt", std::ios_base::app | std::ios_base::out);
                out_SigmaY << timesum << "    " << dt << "    " << SigmaY.row(line_fault) << std::endl; // might be .row() ???
                out_SigmaY.close();
                
                // ========== save SII
                std::ofstream out_Sii("output/EVO_Sii.txt", std::ios_base::app | std::ios_base::out);
                out_Sii << timesum << "    " << dt << "    " << SII_fault.row(line_fault) << std::endl; // might be .row() ???
                out_Sii.close();
                
                // ========== save Theta
                std::ofstream out_Theta("output/EVO_Theta.txt", std::ios_base::app | std::ios_base::out);
                out_Theta << timesum << "    " << dt << "    " << OM.row(line_fault) << std::endl; // might be .row() ???
                out_Theta.close();
                
                // ========== save viscous compaction
                std::ofstream out_Visc_comp("output/EVO_Visc_comp.txt", std::ios_base::app | std::ios_base::out);
                out_Visc_comp << timesum << "    " << dt << "    " << VIS_COMP.row(line_fault) << std::endl; // might be .row() ???
                out_Visc_comp.close();
                
                // ========== save elastic compaction
                std::ofstream out_Elast_comp("output/EVO_Elast_comp.txt", std::ios_base::app | std::ios_base::out);
                out_Elast_comp << timesum << "    " << dt << "    " << EL_DECOM.row(line_fault) << std::endl; // might be .row() ???
                out_Elast_comp.close();

                // ========== save vx Darcy
                std::ofstream out_EVO_vxD("output/EVO_vxD.txt", std::ios_base::app | std::ios_base::out);
                out_EVO_vxD << timesum << "    " << dt << "    " << vxD.row(line_fault) << std::endl; // might be .row() ???
                out_EVO_vxD.close();
                
                // ========== save time, dt, vmax
                std::ofstream out_data("output/EVO_data.txt", std::ios_base::app | std::ios_base::out);
                Vmax = VSLIPB.maxCoeff();
                out_data << std::setw(20) << timesum << std::setw(20) << dt << std::setw(20) << Vmax << std::setw(20) << ynlast << std::setw(20) << iterstep << std::endl;
                out_data.close();
                
                timesum_plus = timesum;
            }
        }
        
        if (fix(timestep / savematstep) * savematstep == timestep) {
            std::ofstream out_file;
            out_file.open("data/file.txt");
            out_file << timestep << std::endl;
            out_file.close();

            std::ofstream out_save_timestep;
            std::string save_file_name = nname + std::to_string(timestep);
            out_save_timestep.open(save_file_name);
            // std::cout << everything << std::endl;
            out_save_timestep.close();
        }
        
    }

    return 0;
}