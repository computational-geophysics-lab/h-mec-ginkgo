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
#include <math.h>
#include <iomanip>
#include <chrono>
#include <H5Cpp.h>

#include "hdf5.hpp"
#include "constants.hpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;
using namespace H5;

typedef Triplet<double> Trp;
typedef high_resolution_clock hrc;
typedef VectorXd VecXd;
typedef MatrixXd MatXd;

// ====================================================================================
// global variable declarations
// ====================================================================================

// Basic nodes
MatXd OM0(Ny, Nx);  // Old state parameter
MatXd OM(Ny, Nx);   // State parameter
MatXd OM5(Ny, Nx);
MatXd ARSF(Ny, Nx); // a - parameter of RSF
MatXd BRSF(Ny, Nx); // b - parameter of RSF
MatXd LRSF(Ny, Nx); // L - parameter of RSF

// Unknown parameters
MatXd pt(Ny1, Nx1);  // Total pressure
MatXd vxs(Ny1, Nx1); // Solid vx - velocity
MatXd vys(Ny1, Nx1); // Solid vy - velocity
MatXd vzs(Ny1, Nx1); // Solid vz - velocity
MatXd pf(Ny1, Nx1);  // Fluid pressure
MatXd vxD(Ny1, Nx1); // Darsi vx - velocity
MatXd vyD(Ny1, Nx1); // Darsi vy - velocity

// Nodal matrices
// Basic nodes
MatXd RHO(Ny, Nx), ETA(Ny, Nx), ETA0(Ny, Nx), ETA1(Ny, Nx), ETA5(Ny, Nx), ETA00(Ny, Nx), IETAPLB(Ny, Nx), SXY(Ny, Nx), SXY0(Ny, Nx), SZX0(Ny, Nx), SZY0(Ny, Nx), YNY0(Ny, Nx), KKK(Ny, Nx),
      GGG(Ny, Nx), COHC(Ny, Nx), COHT(Ny, Nx), FRIC(Ny, Nx), FRIT(Ny, Nx), DILC(Ny, Nx), TTT(Ny, Nx), EIIB(Ny, Nx), VSLIPB(Ny, Nx), SZX(Ny, Nx), SZY(Ny, Nx);

// Pressure nodes
MatXd ETAB(Ny1, Nx1), ETAB0(Ny1, Nx1), ETAP(Ny1, Nx1), ETAP0(Ny1, Nx1), POR(Ny1, Nx1), GGGP(Ny1, Nx1), GGGB(Ny1, Nx1), PTF0(Ny1, Nx1), PT0(Ny1, Nx1), PF0(Ny1, Nx1), pt_ave(Ny1, Nx1),
      pf_ave(Ny1, Nx1), SXX(Ny1, Nx1), SXX0(Ny1, Nx1), SYY(Ny1, Nx1), SYY0(Ny1, Nx1), DILP(Ny1, Nx1);
// Vx nodes
MatXd RHOX(Ny1, Nx1), RHOFX(Ny1, Nx1), ETADX(Ny1, Nx1), PORX(Ny1, Nx1), VX0(Ny1, Nx1), VXF0(Ny1, Nx1);
// Vy nodes
MatXd RHOY(Ny1, Nx1), RHOFY(Ny1, Nx1), ETADY(Ny1, Nx1), PORY(Ny1, Nx1), VY0(Ny1, Nx1), VYF0(Ny1, Nx1);
// Vz nodes
MatXd VZ0(Ny1, Nx1);

MatXd ESP(Ny, Nx), EXY(Ny, Nx), EZX(Ny, Nx), EZY(Ny, Nx), EXX(Ny1, Nx1), EYY(Ny1, Nx1), EII(Ny1, Nx1), EIIVP(Ny1, Nx1), SII(Ny1, Nx1), DSII(Ny1, Nx1), DIS(Ny1, Nx1);

MatXd EL_DECOM(Ny1, Nx1);   // Elastic (de)compaction
MatXd VIS_COMP(Ny1, Nx1);

// (3) Defining global matrixes
// according to the global number of unknowns
// Sparse Matrix L is not yet defined as it will be built from a set of Triplets each step
VecXd R(N); // Vector of the right parts of equations
VecXd S(N);
PardisoLU<SparseMatrix<double>> solver;

// variable type declaration
int ynlast, iterstep;

double timesum = 0;
double dt = dtelastic0;
bool yndtdecrease = true;

double dt00, dtx, dty, dtz, dtlapusta, Vmax, maxvxy;

VecXd DSYLSQ(niterglobal);

MatXd DVX0(Ny1, Nx1), DVY0(Ny1, Nx1), DVZ0(Ny1, Nx1), DSY(Ny, Nx), YNY(Ny, Nx), SigmaY(Ny, Nx), SII_fault(Ny, Nx), SIIB(Ny, Nx); // DVX0, DVY0, DVZ0 unused

VecXd timesumcur(num_timesteps), dtcur(num_timesteps);
VecXd maxvxsmod(num_timesteps), minvxsmod(num_timesteps), maxvysmod(num_timesteps), minvysmod(num_timesteps), maxvzsmod(num_timesteps), minvzsmod(num_timesteps);

// ====================================================================================
// end global variable declarations
// ====================================================================================

// lambda function that rounds towards 0
auto fix = [](const double temp) {
    return temp < 0. ? (int)ceil(temp) : (int)floor(temp);
};

// lambda function that checks if a value is between a min and max value
auto enforce_bounds = [](double val, const double min, const double max) {
    if (val < min) {
        val = min;
    } else if (val > max) {
        val = max;
    }
    return val;
};

// lambda function that squares each element of a 2x2 block
auto square_block = [](Matrix2d mat) {
    return (pow(mat(0, 0), 2) + pow(mat(0, 1), 2) + pow(mat(1, 0), 2) + pow(mat(1, 1), 2));
};

// lambda function that computes: x / (x + y)
auto divplus = [](const double x, const double y) {
    return (x / (x + y));
};

// function that checks if a integer value is bigger than 0 and lower than a set bound
int check_bounds(int k, const int bound) {
    if (k < 0) {
        k = 0;
    } else if (k > bound - 2) {
        k = bound - 2;
    }
    return k;
}

// function that sets the boundary values to values of one column/row further in
void copy_bounds(MatXd& Temp) {
    Temp.col(0) = Temp.col(1);
    Temp.col(Nx) = Temp.col(Nx - 1);
    Temp.row(0) = Temp.row(1);
    Temp.row(Ny) = Temp.row(Ny - 1);
}

// function that removes already calculated results in some output files so that these are in order (of timesteps)
void remove_double_output(const string file_name, const int num) {
    ifstream file(file_name);
    string new_content = "";
    string line;
    for (int i = 0; i < num; i++) {
        getline(file, line);
        new_content += line + '\n';
    }
    file.close();
    ofstream file_out(file_name);
    file_out << new_content;
    file_out.close();
}

int main() {
    Eigen::initParallel();
    // ====================================================================================
    // initialize random generator with set seed for testing purpose
    // srand(42);
    // initialize random generator with current time as seed
    srand(time(nullptr));
    // ====================================================================================

    // read input
    // Load file
    string timestep_str;
    ifstream input_timestep("file.txt");
    getline(input_timestep, timestep_str);
    int timestep = stoi(timestep_str);
    input_timestep.close();

    for (int m = 0; m < marknum; m++) {
        // Define randomized regular coordinates
        xm(m) = xbeg + floor(m / Ny_markers) * dxms + (rand() % 1) * dxms;
        ym(m) = ybeg + m % Ny_markers * dyms + (rand() % 1) * dyms;
        // Matrix
        porm(m) = .01 * (1 + .0 * (rand() % 1 - .5));
        etam(m) = etasm(m) * exp(-alpha * porm(m));

        // Air, wedge, slab
        if (ym(m) < upper_block || ym(m) > lower_block) {
            t_marker(m) = -1;
            etam(m) = 1e23;
            etasm(m) = 1e23;
            rhom(m) = 2800;
            kkkm(m) = 2e-16; // * (dy_faultw)^2;
            cohescm(m) = cohes * 1e3;
            cohestm(m) = cohes * 1e3;
            frictcm(m) = .8;
            dilatcm(m) = dilatation;
            fricttm(m) = tensile;
        }
    }

    if (timestep > 0) {
        const string filename = nname + to_string(timestep) + ".h5";
        const string group_matrix = "Matrix";
        const string group_vector = "Vector";
        const string group_values = "Value";

        const hsize_t dims1[2] = {Ny, Nx};
        const string matrix_names[24] = {"SIIB", "OM0", "OM", "ARSF", "BRSF", "LRSF", "RHO", "ETA0", "ETA1", "ETA5", "ETA00", "IETAPLB", "SXY0", "YNY0", "KKK",
                                   "GGG", "COHC", "COHT", "FRIC", "FRIT", "DILC", "TTT", "EIIB", "VSLIPB"}; // {"names"} has to be the same as in *matrix
        int j = 0;
        for (auto i : {&SIIB, &OM0, &OM, &ARSF, &BRSF, &LRSF, &RHO, &ETA0, &ETA1, &ETA5, &ETA00, &IETAPLB, &SXY0, &YNY0, &KKK, &GGG, &COHC, &COHT, &FRIC, &FRIT,
                       &DILC, &TTT, &EIIB, &VSLIPB}) { // {names} *matrix
            *i = read_matrix(filename, group_matrix, matrix_names[j], dims1);
            j++;
        }

        const hsize_t dims2[2] = {Ny1, Nx1};
        const string matrix_names_plus[30] = {"pt", "vxs", "vys", "pf", "vxD", "vyD", "ETAB", "ETAB0", "ETAP", "ETAP0", "POR", "GGGP", "GGGB", "PTF0", "PT0", "PF0",
                                        "SXX0", "SYY0", "RHOX", "RHOFX", "ETADX", "PORX", "VX0", "VXF0", "RHOY", "RHOFY", "ETADY", "PORY", "VY0", "VYF0"}; // {"names"} has to be the same as in *matrix_plus
        j = 0;
        for (auto i : {&pt, &vxs, &vys, &pf, &vxD, &vyD, &ETAB, &ETAB0, &ETAP, &ETAP0, &POR, &GGGP, &GGGB, &PTF0, &PT0, &PF0, &SXX0, &SYY0,
                       &RHOX, &RHOFX, &ETADX, &PORX, &VX0, &VXF0, &RHOY, &RHOFY, &ETADY, &PORY, &VY0, &VYF0}) { // {names} *matrix_plus
            *i = read_matrix(filename, group_matrix, matrix_names_plus[j], dims2);
            j++;
        }
        
        const hsize_t dim1[1] = {num_timesteps};
        const string vector_names[6] = {"timesumcur", "dtcur", "maxvxsmod", "minvxsmod", "maxvysmod", "minvysmod"}; // {"names"} has to be the same as in *vec
        j = 0;
        for (auto i : {&timesumcur, &dtcur, &maxvxsmod, &minvxsmod, &maxvysmod, &minvysmod}) { // {names} *vec
            *i = read_vector(filename, group_vector, vector_names[j], dim1);
            j++;
        }

        const hsize_t dim3[1] = {marknum};
        const string vector_names_marker[5] = {"xm", "ym", "sxxm", "syym", "sxym"}; // {"names"} has to be the same as in *vec2
        j = 0;
        for (auto i : {&xm, &ym, &sxxm, &syym, &sxym}) { // {names} *vec2
            *i = read_vector(filename, group_vector, vector_names_marker[j], dim3);
            j++;
        }

        const hsize_t dim2[1] = {9};
        const VecXd temp = read_vector(filename, group_values, "values", dim2);
        timesum = temp(0); dt00 = temp(1); dtx = temp(2); dty = temp(3); dtlapusta = temp(4); Vmax = temp(5); maxvxy = temp(6); dt = temp(7), yndtdecrease = temp(8);

        if (antiplane) {
            const string group_antiplane = "Antiplane";
            SZX0 = read_matrix(filename, group_antiplane, "SZX0", dims1);
            SZY0 = read_matrix(filename, group_antiplane, "SZY0", dims1);

            VZ0 = read_matrix(filename, group_antiplane, "VZ0", dims2);

            maxvzsmod = read_vector(filename, group_antiplane, "maxvzsmod", dim1);
            minvzsmod = read_vector(filename, group_antiplane, "minvzsmod", dim1);
        }

        timestep++;

        int num_lines = 0;
        string line;
        ifstream myfile("EVO_data.txt");
        while (getline(myfile, line)) {++num_lines;}
        myfile.close();
        
        if (num_lines >= timestep) {
            const string files[11] = {"EVO_Vslip.txt", "EVO_viscosity.txt", "EVO_press_flu.txt", "EVO_press_eff.txt", "EVO_SigmaY.txt", "EVO_Sii.txt", "EVO_Theta.txt", 
                                      "EVO_Visc_comp.txt", "EVO_Elast_comp.txt", "EVO_VxD.txt", "EVO_data.txt"};
            for (int i = 0; i < 11; i++) {
                remove_double_output(files[i], timestep - 1);
            }
        }
    } else {
        timestep = 1;
        // Basic nodes
        OM0 = MatXd::Constant(Ny, Nx, omm(0));    // Old state parameter
        ARSF = MatXd::Constant(Ny, Nx, arsfm(1)); // a - parameter of RSF
        BRSF = MatXd::Constant(Ny, Nx, brsfm(1)); // b - parameter of RSF
        LRSF = MatXd::Constant(Ny, Nx, lrsfm(0)); // L - parameter of RSF
        
        // set matrices to 0
        for (auto i : {pt, vxs, vys, pf, vxD, vyD, RHO, ETA0, IETAPLB, SXY, SXY0, SZX0, SZY0, YNY0, KKK, GGG, COHC, COHT, FRIC, FRIT, DILC, TTT, EIIB, ETAB, ETAB0, ETAP, ETAP0, POR, GGGP,
                       GGGB, PTF0, pt_ave, pf_ave, SXX, SXX0, SYY, SYY0, RHOX, RHOFX, ETADX, PORX, VX0, VXF0, RHOY, RHOFY, ETADY, PORY, VY0, VYF0, VSLIPB, VZ0}) {
            i.setZero();
        }

        // Define Fault
        // #pragma omp parallel for (slower with n = 4)
        for (int i = 0; i < Ny; i++) {
            if (y(i) > upper_block && y(i) < lower_block) {
                for (int j = 0; j < Nx; j++) {
                    OM0(i, j) = omm(1);
                
                    if (x(j) < TS_1) {
                        BRSF(i, j) = brsfm(0);
                        ARSF(i, j) = arsfm(0);
                    }
                    if (x(j) >= TS_1 && x(j) < TS_2) {
                        BRSF(i, j) = brsfm(0) - (brsfm(0) - brsfm(1)) * ((x(j) - TS_1) / (TS_2 - TS_1));
                        ARSF(i, j) = arsfm(0) - (arsfm(0) - arsfm(1)) * ((x(j) - TS_1) / (TS_2 - TS_1));
                        LRSF(i, j) = lrsfm(0) - (lrsfm(0) - lrsfm(1)) * ((x(j) - TS_1) / (TS_2 - TS_1));
                    }
                    if (x(j) >= TS_2 && x(j) <= TS_3) {
                        LRSF(i, j) = lrsfm(1);
                    }
                    if (x(j) > TS_3 && x(j) <= TS_4) {
                        BRSF(i, j) = brsfm(1) - (brsfm(1) - brsfm(0)) * ((x(j) - TS_3) / (TS_4 - TS_3));
                        ARSF(i, j) = arsfm(1) - (arsfm(1) - arsfm(0)) * ((x(j) - TS_3) / (TS_4 - TS_3));
                        LRSF(i, j) = lrsfm(1) - (lrsfm(1) - lrsfm(0)) * ((x(j) - TS_3) / (TS_4 - TS_3));
                    }
                    if (x(j) > TS_4) {
                        BRSF(i, j) = brsfm(0);
                        ARSF(i, j) = arsfm(0);
                    }
                }
            }
        }
        
        OM = OM0;

        // set vectors to 0
        for (auto i : {sxxm, syym, sxym}) {
            i.setZero();
        }

        PT0 = MatXd::Constant(Ny1, Nx1, PCONF + PTFDIFF);
        PF0 = MatXd::Constant(Ny1, Nx1, PCONF);
    }

    // declaration of matrices that are set to zero each timestep
    MatXd ETA0SUM(Ny, Nx), COHCSUM(Ny, Nx), FRICSUM(Ny, Nx), DILCSUM(Ny, Nx), COHTSUM(Ny, Nx), FRITSUM(Ny, Nx), WTSUM(Ny, Nx);

    // Interpolate ETA, RHO to nodal points
    // Basic nodes
    MatXd RHOSUM(Ny, Nx), ETASUM(Ny, Nx), KKKSUM(Ny, Nx), TTTSUM(Ny, Nx), SXYSUM(Ny, Nx), GGGSUM(Ny, Nx);

    //LDZ
    MatXd OM0SUM(Ny, Nx);  // Old state parameter
    MatXd OMSUM(Ny, Nx);   // State parameter
    MatXd ARSFSUM(Ny, Nx); // a - parameter of RSF
    MatXd BRSFSUM(Ny, Nx); // b - parameter of RSF
    MatXd LRSFSUM(Ny, Nx); // L - parameter of RSF

    // Pressure nodes
    MatXd ETAPSUM(Ny1, Nx1), ETAP0SUM(Ny1, Nx1), ETAB0SUM(Ny1, Nx1), PORSUM(Ny1, Nx1), SXXSUM(Ny1, Nx1), SYYSUM(Ny1, Nx1), GGGPSUM(Ny1, Nx1), WTPSUM(Ny1, Nx1);
    // Vx nodes
    MatXd RHOXSUM(Ny1, Nx1), RHOFXSUM(Ny1, Nx1), ETADXSUM(Ny1, Nx1), PORXSUM(Ny1, Nx1), WTXSUM(Ny1, Nx1);
    // Vy nodes
    MatXd RHOYSUM(Ny1, Nx1), RHOFYSUM(Ny1, Nx1), ETADYSUM(Ny1, Nx1), PORYSUM(Ny1, Nx1), WTYSUM(Ny1, Nx1);

    MatXd ETA50(Ny, Nx);

    // nucleation size
    const double hstar = pi / 2. * shearmod * lrsfm(1) * brsfm(1) / pow(brsfm(1) - arsfm(1), 2) / PTFDIFF;
    // cohesive zone size
    const double coh = 9. / 32. * pi * shearmod * lrsfm(1) / brsfm(1) / PTFDIFF;
    // Print information about discretization
    cout << ">> VW width = " << (((TS_4 + TS_3) / 2) - ((TS_2 + TS_1) / 2)) / 1e3 << " (km)" << endl;
    cout << ">> Critical nucleation size = " << hstar / 1e3 << " (km)" << endl;
    cout << ">> Cohesive zone = " << coh << " (m)" << endl;

    // /////////////////////////////////////////////////////////////////////////////////////// 
    // actual computations start here
    // /////////////////////////////////////////////////////////////////////////////////////// 

    for (; timestep <= num_timesteps; timestep++) {
        for (auto i : {RHOSUM, ETASUM, KKKSUM, TTTSUM, SXYSUM, GGGSUM, ETA, ETA0SUM, COHCSUM, FRICSUM, DILCSUM, COHTSUM, FRITSUM, WTSUM, 
                       OM0SUM, OMSUM, ARSFSUM, BRSFSUM, LRSFSUM, ETAPSUM, ETAP0SUM, ETAB0SUM, PORSUM, SXXSUM, SYYSUM, GGGPSUM, WTPSUM,
                       RHOXSUM, RHOFXSUM, ETADXSUM, PORXSUM, WTXSUM, RHOYSUM, RHOFYSUM, ETADYSUM, PORYSUM, WTYSUM}) {
            i.setZero();
        }
        
        // Cycle on markers
        #pragma omp parallel for // about 3-4x faster with n = 4
        for (int m = 0; m < marknum; m++) {
            double cohescmm, cohestmm, frictcmm, dilatcmm, fricttmm, etasmm0, etamm0, etamm, rhomm, etadm;

            // Marker properties
            double kkkmm = kkkm(m) * pow(porm(m) / POR0, 3);
            // Checking permeability limits
            kkkmm = enforce_bounds(kkkmm, kkkmin, kkkmax);

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
            // Matrix viscosity
            etamm0 = enforce_bounds(etamm0, etamin, etamax);
            // Effective viscosity
            etamm = enforce_bounds(etamm, etamin, etamax);
            
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

            Matrix2d wtm;
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
        
        OM = OM0;
        
        // Multiple solving of equations
        if (!yndtdecrease) {
            dt = max(min(dt * dtkoefup, dtelastic0), dtmin);
        } else {
            dt = max(min(dt, dtelastic0), dtmin);
        }
        
        yndtdecrease = false;
        dt00 = dt;
        ynlast = 0;
        DSYLSQ.setZero();
        
        for (iterstep = 0; iterstep < niterglobal; iterstep++) {
            // Limiting viscosity
            double etamincur = dt * shearmod * 1e-4;
            double ptscale, pfscale;
            
            // External P - nodes: symmetry
            copy_bounds(pt);
            copy_bounds(pf);
            
            // Basic nodes
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    if (ETA(i, j) < etamincur) {
                        ETA(i, j) = etamincur;
                    }
                    // Compute plastic strain rate
                    if (ETA(i, j) < ETA0(i, j)) {
                        const double SXX_temp = SXX.block(i, j, 2, 2).sum(), SYY_temp = SYY.block(i, j, 2, 2).sum();
                        SIIB(i, j) = sqrt(pow(SXY(i, j), 2) + (pow(SXX_temp, 2) + pow(SYY_temp, 2) + pow(SXX_temp + SYY_temp, 2)) / 32.);
                        IETAPLB(i, j) = (1. / ETA(i, j) - 1. / ETA0(i, j));
                        EIIB(i, j) = dy_faultw * SIIB(i, j) * IETAPLB(i, j) / 2.;
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
                    const double IETAPL = IETAPLB.block(i - 1, j - 1, 2, 2).sum() / 4.;
                    if (YNY0(i - 1, j - 1) > 0 || YNY0(i, j - 1) > 0 || YNY0(i - 1, j) > 0 || YNY0(i, j) > 0) {
                        ETAP(i, j) = 1. / (1. / ETAP0(i, j) + IETAPL);
                        ETAB(i, j) = 1. / (1. / ETAB0(i, j) + dy_faultw * IETAPL * POR(i, j));
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
                    // double ss3 = min(max((pt(i, j) - pf(i, j)) * 1e-6, 0.), 100.); // SIGMA3, MPa
                    // double aa = aa1 + aa2 * exp(-ss3 / aa3);
                    // double bb = bb1 + bb2 * exp(-ss3 / bb3);
                    // double cc = cc1 + cc2 / 100. * pow(ss3, cc3);
                    // double dil = sin(aa * bb * (exp(-bb * gammap) - exp(-cc * gammap)) / (cc - bb) / 180. * pi);
                    DILP(i, j) = 0; // 2 * (dil * EIIB(i - 1, j - 1) + dil * EIIB(i, j - 1) + dil * EIIB(i - 1, j) + dil * EIIB(i, j)) / 4;
                }
            }
            
            if (dt > 2e5) {
                pfscale = ETADX.block(1, 0, Ny - 1, Nx).minCoeff() * dx * 1e17 / pow(dt, 2);
            } else {
                pfscale = ETADX.block(1, 0, Ny - 1, Nx).minCoeff() * dx;
            }
            
            ptscale = pfscale;
            
            // This loop here is now the most timeconsuming part -> next step: try to replace coeffRef if possible, maybe use Triplet format
            // 5)Composing global matrixes L(), R()

            vector<Trp> Trip; // define triplet list to build sparse matrix
            Trip.reserve(N * 10); // reserving memory space for all the entries to be stored with some reserve

            R.setZero();

            for (int j = 0; j < Nx1; j++) {
                for (int i = 0; i < Ny1; i++) {
                    // Computing global indexes for vx, vy, p
                    const int kp = (j * Ny1 + i) * Num_var;
                    const int kx = kp + 1;
                    const int ky = kp + 2;
                    const int kpf = kp + 3;
                    const int kxf = kp + 4;
                    const int kyf = kp + 5;
                    const int kz = kp + 6;
                    
                    // 5a) Composing equation for vxs
                    if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (j == Nx) {
                            Trip.push_back(Trp(kx, kx, 1));
                        }
                        
                        // Upper boundary
                        // prescribed velocity
                        if (i == 0 && j < Nx) {
                            Trip.insert(Trip.end(), {Trp(kx, kx, 1), Trp(kx, kx + Num_var, 1)});
                            R(kx) = 2 * bcupper;
                        }
                        
                        // Lower boundary
                        // prescribed velocity
                        if (i == Ny && j < Nx) {
                            Trip.insert(Trip.end(), {Trp(kx, kx, 1), Trp(kx, kx - Num_var, 1)});
                            R(kx) = 2 * bclower;
                        }
                        
                        // Left boundary
                        if (j == 0 && i > 0 && i < Ny) {
                            Trip.insert(Trip.end(), {Trp(kx, kx, 1), Trp(kx, kx + Num_var * Ny1, -1)});
                        }
                        
                        // Right boundary
                        if (j == Nx - 1 && i > 0 && i < Ny) {
                            Trip.insert(Trip.end(), {Trp(kx, kx, 1), Trp(kx, kx - Num_var * Ny1, -1)});
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
                        double ETAXY1 = ETA(i - 1, j), ETAXY2 = ETA(i, j);
                        double ETAXX1 = ETAP(i, j), ETAXX2 = ETAP(i, j + 1);
                        // Shear modulus
                        const double GXY1 = GGG(i - 1, j), GXY2 = GGG(i, j);
                        const double GXX1 = GGGP(i, j), GXX2 = GGGP(i, j + 1);
                        // Viscoelasticity factor
                        const double KXY1 = divplus(dt * GXY1, ETAXY1), KXY2 = divplus(dt * GXY2, ETAXY2);
                        const double KXX1 = divplus(dt * GXX1, ETAXX1), KXX2 = divplus(dt * GXX2, ETAXX2);
                        // Numerical viscosity
                        ETAXY1 *= KXY1;
                        ETAXY2 *= KXY2;
                        ETAXX1 *= KXX1;
                        ETAXX2 *= KXX2;
                        // Numerical stresses
                        const double SXY1 = SXY0(i - 1, j) * (1 - KXY1), SXY2 = SXY0(i, j) * (1 - KXY2);
                        const double SXX1 = SXX0(i, j) * (1 - KXX1), SXX2 = SXX0(i, j + 1) * (1 - KXX2);
                        // Density derivatives
                        const double dRHOdx = (RHOX(i, j + 1) - RHOX(i, j - 1)) / (2. * dx);
                        const double dRHOdy = (RHO(i, j) - RHO(i - 1, j)) / dy;
                        // Left part
                        const double temp = gx * dt * dRHOdy / 4.;
                        Trip.insert(Trip.end(), {Trp(kx, kx, -(ETAXX1 + ETAXX2) / dx2 - (ETAXY1 + ETAXY2) / dy2 - gx * dt * dRHOdx - inertia * RHOX(i, j) / dt), Trp(kx, kx - Ny1 * Num_var, ETAXX1 / dx2), 
                                    Trp(kx, kx + Ny1 * Num_var, ETAXX2 / dx2), Trp(kx, kx - Num_var, ETAXY1 / dy2), Trp(kx, kx + Num_var, ETAXY2 / dy2), /* vxs3, vxs1, vxs5, vxs2, vxs4 */
                                    Trp(kx, ky - Num_var, (ETAXY1 - ETAXX1) / dx_dy - temp), Trp(kx, ky, (ETAXX1 - ETAXY2) / dx_dy - temp),
                                    Trp(kx, ky - Num_var + Ny1 * Num_var, (ETAXX2 - ETAXY1) / dx_dy - temp), Trp(kx, ky + Ny1 * Num_var, (ETAXY2 - ETAXX2) / dx_dy - temp), /* vys1, vys2, vys3, vys4 */
                                    Trp(kx, kp, ptscale / dx), Trp(kx, kp + Ny1 * Num_var, -ptscale / dx)}); /* Pt1', Pt2' */
                        // Right part
                        R(kx) = -RHOX(i, j) * (inertia * VX0(i, j) / dt + gx) - (SXX2 - SXX1) / dx - (SXY2 - SXY1) / dy;
                    }
                    
                    // 5b) Composing equation for vys
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vys = 0
                        if (i == Ny) {
                            Trip.push_back(Trp(ky, ky, 1));
                        }
                        
                        // Left boundary
                        // Free Slip
                        if (j == 0) {
                            Trip.insert(Trip.end(), {Trp(ky, ky, 1), Trp(ky, ky + Ny1 * Num_var, 1)});
                        }
                        
                        // Right boundary
                        // Free Slip
                        if (j == Nx) {
                            Trip.insert(Trip.end(), {Trp(ky, ky, 1), Trp(ky, ky - Ny1 * Num_var, 1)});
                        }
                        
                        // Upper boundary: no penetration
                        if (i == 0 && j > 0 && j < Nx) {
                            Trip.push_back(Trp(ky, ky, 1));
                        }
                        
                        // Lower boundary: no penetration
                        if (i == Ny - 1 && j > 0 && j < Nx) {
                            Trip.push_back(Trp(ky, ky, 1));
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
                        double ETAXY1 = ETA(i, j - 1), ETAXY2 = ETA(i, j);
                        double ETAYY1 = ETAP(i, j), ETAYY2 = ETAP(i + 1, j);
                        // Shear modulus
                        const double GXY1 = GGG(i, j - 1), GXY2 = GGG(i, j);
                        const double GYY1 = GGGP(i, j), GYY2 = GGGP(i + 1, j);
                        // Viscoelasticity factor
                        const double KXY1 = divplus(dt * GXY1, ETAXY1), KXY2 = divplus(dt * GXY2, ETAXY2);
                        const double KYY1 = divplus(dt * GYY1, ETAYY1), KYY2 = divplus(dt * GYY2, ETAYY2);
                        // Numerical viscosity
                        ETAXY1 *= KXY1;
                        ETAXY2 *= KXY2;
                        ETAYY1 *= KYY1;
                        ETAYY2 *= KYY2;
                        // Numerical stresses
                        const double SXY1 = SXY0(i, j - 1) * (1 - KXY1), SXY2 = SXY0(i, j) * (1 - KXY2);
                        const double SYY1 = SYY0(i, j) * (1 - KYY1), SYY2 = SYY0(i + 1, j) * (1 - KYY2);
                        // Density derivatives
                        const double dRHOdy = (RHOY(i + 1, j) - RHOY(i - 1, j)) / 2. / dy;
                        const double dRHOdx = (RHO(i, j) - RHO(i, j - 1)) / dx;
                        // Left part
                        const double temp = gy * dt * dRHOdx / 4.;
                        Trip.insert(Trip.end(), {Trp(ky, ky, -(ETAYY1 + ETAYY2) / dy2 - (ETAXY1 + ETAXY2) / dx2 - gy * dt * dRHOdy - inertia * RHOY(i, j) / dt), Trp(ky, ky - Ny1 * Num_var, ETAXY1 / dx2),
                                    Trp(ky, ky + Ny1 * Num_var, ETAXY2 / dx2), Trp(ky, ky - Num_var, ETAYY1 / dy2), Trp(ky, ky + Num_var, ETAYY2 / dy2), /* vys3, vys1, vys5, vys2, vys4 */
                                    Trp(ky, kx - Ny1 * Num_var, (ETAXY1 - ETAYY1) / dx_dy - temp), Trp(ky, kx + Num_var - Ny1 * Num_var, (ETAYY2 - ETAXY1) / dx_dy - temp),
                                    Trp(ky, kx, (ETAYY1 - ETAXY2) / dx_dy - temp), Trp(ky, kx + Num_var, (ETAXY2 - ETAYY2) / dx_dy - temp), /* vxs1, vxs2, vxs3, vxs4 */
                                    Trp(ky, kp, ptscale / dy), Trp(ky, kp + Num_var, -ptscale / dy)}); /* Pt1', Pt2' */
                        // Right part
                        R(ky) = -RHOY(i, j) * (inertia * VY0(i, j) / dt + gy) - (SYY2 - SYY1) / dy - (SXY2 - SXY1) / dx;
                    }

                    // 5c) Composing equation for Pt
                    if (i == 0 || j == 0 || i == Ny || j == Nx) { // || (i == 2 && j == 2))
                        // BC equation: 1 * Pt = 0
                        Trip.push_back(Trp(kp, kp, 1));
                    } else {
                        // Solid Continuity: dVxs / dx + dVys / dy + (Pt - Pf) / ETAbulk = 0
                        //              vys1
                        //               |
                        //        vxs1 -- Pt, Pf -- vxs2
                        //               |
                        //              vys2
                        // Drained compressibility
                        const double BETADRAINED = (1. / GGGB(i, j) + BETASOLID) / (1 - POR(i, j));
                        // Biott - Willis koefficient
                        const double KBW = 1 - BETASOLID / BETADRAINED;
                        // Left part
                        Trip.insert(Trip.end(), {Trp(kp, kx - Ny1 * Num_var, -1. / dx), Trp(kp, kx, 1. / dx), /* vxs1, vxs2 */ Trp(kp, ky - Num_var, -1. / dy), Trp(kp, ky, 1. / dy), /* vys1, vys2 */
                                    Trp(kp, kp, ptscale * (1. / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED / dt)), Trp(kp, kpf, -pfscale * (1. / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / dt))}); /* Pt, Pf */
                        // Right part
                        R(kp) = BETADRAINED * (PT0(i, j) - KBW * PF0(i, j)) / dt + DILP(i, j);
                    }
                    
                    // 5d) Composing equation for vxD
                    if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (j == Nx) {
                            Trip.push_back(Trp(kxf, kxf, 1));
                        }
                        
                        // Upper boundary: symmetry
                        if (i == 0 && j < Nx) {
                            Trip.insert(Trip.end(), {Trp(kxf, kxf, 1), Trp(kxf, kxf + Num_var, -1)});
                        }
                        
                        // Lower boundary: symmetry
                        if (i == Ny && j < Nx) {
                            Trip.insert(Trip.end(), {Trp(kxf, kxf, 1), Trp(kxf, kxf - Num_var, -1)});
                        }
                        
                        // Left boundary
                        // no penetration
                        if (j == 0) {
                            Trip.push_back(Trp(kxf, kxf, 1));
                        }
                        
                        // Right boundary
                        // no penetration
                        if (j == Nx - 1) {
                            Trip.push_back(Trp(kxf, kxf, 1));
                        }
                        
                    } else {
                        // Fluid X - Darsi: - ETAfluid / K * VxD - dPf / dx = -RHOf * gx + RHOf * DVxs / Dt
                        //
                        //  Pf1 --- vxD, vxs --- Pf2
                        //
                        // Left part
                        Trip.insert(Trip.end(), {Trp(kxf, kxf, -ETADX(i, j) - RHOFX(i, j) / PORX(i, j) * inertia / dt), Trp(kxf, kx, -RHOFX(i, j) * inertia / dt), /* vxD, vxs */
                                    Trp(kxf, kpf, pfscale / dx), Trp(kxf, kpf + Ny1 * Num_var, -pfscale / dx)}); /* Pf1', Pf2' */
                        // Right part
                        R(kxf) = -RHOFX(i, j) * (inertia * VXF0(i, j) / dt + gx);
                    }
                    
                    // 5e) Composing equation for vyD
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (i == Ny) {
                            Trip.push_back(Trp(kyf, kyf, 1));
                        }
                        
                        // Left boundary
                        // symmetry
                        if (j == 0 && i > 0 && i < Ny - 1) {
                            Trip.insert(Trip.end(), {Trp(kyf, kyf, 1), Trp(kyf, kyf + Ny1 * Num_var, -1)});
                        }
                        
                        // Right boundary
                        // symmetry
                        if (j == Nx && i > 0 && i < Ny - 1) {
                            Trip.insert(Trip.end(), {Trp(kyf, kyf, 1), Trp(kyf, kyf - Ny1 * Num_var, -1)});
                        }
                        
                        // Upper boundary: no penetration
                        if (i == 0) {
                            Trip.push_back(Trp(kyf, kyf, 1));
                            R(kyf) = bcvyflower;
                        }
                        
                        // Lower boundary: no penetration
                        if (i == Ny - 1) {
                            Trip.push_back(Trp(kyf, kyf, 1));
                            R(kyf) = bcvyflower;
                        }
                    } else {
                        // Fluid Y - Darsi: - ETAfluid / K * VyD - dPf / dy = -RHOf * gy + RHOf * DVys / Dt
                        //
                        //   Pf1
                        //    |
                        //   vyD, vy
                        //    |
                        //   Pf2
                        //
                        // Left part
                        Trip.insert(Trip.end(), {Trp(kyf, kyf, -ETADY(i, j) - RHOFY(i, j) / PORY(i, j) * inertia / dt), Trp(kyf, ky, -RHOFY(i, j) * inertia / dt), /* vyD, vys */
                                    Trp(kyf, kpf, pfscale / dy), Trp(kyf, kpf + Num_var, -pfscale / dy)}); /* Pf1', Pf2' */
                        // Right part
                        R(kyf) = -RHOFY(i, j) * (inertia * VYF0(i, j) / dt + gy);
                    }
                                        
                    // 5f) Composing equation for Pf
                    if (j == 0 || j == Nx || i <= 1 || i >= Ny - 1) {
                        // BC equation: 1 * Pf = 0
                        // Real BC
                        if (i == 1 || i == Ny - 1) {
                            Trip.insert(Trip.end(), {Trp(kpf, kpf, pfscale), Trp(kpf, kp, -ptscale)});
                            R(kpf) = -PTFDIFF;
                        } else {
                            Trip.push_back(Trp(kpf, kpf, 1));
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
                        const double BETADRAINED = (1 / GGGB(i, j) + BETASOLID) / (1 - POR(i, j));
                        // Biott - Willis koefficient
                        const double KBW = 1 - BETASOLID / BETADRAINED;
                        // Skempton koefficient
                        const double KSK = (BETADRAINED - BETASOLID) / (BETADRAINED - BETASOLID + POR(i, j) * (BETAFLUID - BETASOLID));
                        // Left part
                        Trip.insert(Trip.end(), {Trp(kpf, kxf - Ny1 * Num_var, -1. / dx), Trp(kpf, kxf, 1. / dx), /* vxs1, vxs2 */ Trp(kpf, kyf - Num_var, -1. / dy), Trp(kpf, kyf, 1. / dy), /* vys1, vys2 */
                                    Trp(kpf, kp, -ptscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / dt)), Trp(kpf, kpf, pfscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / KSK / dt))}); /* Pt, Pf */
                        // Right part
                        R(kpf) = -BETADRAINED * KBW * (PT0(i, j) - 1. / KSK * PF0(i, j)) / dt - DILP(i, j);
                    }
                    if (antiplane) {
                        // 5o) Composing equation for vzs (out-of-plane component)
                        if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                            // Ghost nodes: 1*vzs=0
                            if (j == Nx1) {
                                Trip.push_back(Trp(kz, kz, 1));
                            }
                            
                            // Upper boundary
                            if (i == 1 && j < Nx1) {
                                Trip.insert(Trip.end(), {Trp(kz, kz, 1), Trp(kz, kz + Num_var, 1)});
                                R(kz) = 2 * bclower;
                            }
                            
                            // Lower boundary
                            if (i == Ny1 && j < Nx1) {
                                Trip.insert(Trip.end(), {Trp(kz, kz, 1), Trp(kz, kz - Num_var, 1)});
                                R(kz) = -2 * bclower;
                            }
                            
                            // Left boundary
                            if (j == 1 && i > 1 && i < Ny1) {
                                Trip.insert(Trip.end(), {Trp(kz, kz, 1), Trp(kz, kz + Num_var * Ny1, -1)});
                            }
                            
                            // Right boundary
                            if (j == Nx && i > 1 && i < Ny1) {
                                Trip.insert(Trip.end(), {Trp(kz, kz, 1), Trp(kz, kz - Num_var * Ny1, -1)});
                                R(kz)=0;
                            }
                        } else {
                            // Total Z-Stokes: dSIGMAzxt'/dx+dSIGMAzyt'/dy-dPt/dx=-RHOt*gz
                            // SIGMAijt=2*ETA*EPSILONijs*K+SIGMAijt0*(1-K)
                            //             vzs2
                            //              |
                            //        (P)  SZY   (P)
                            //              |
                            //  vzs1--SZX--vzs3--SZX--vzs5
                            //              |
                            //        (P)  SZY   (P)
                            //              |
                            //             vzs4
                            
                            double ETAXY = ETA(i, j);
                            // Shear modulus
                            const double GXY = GGG(i, j);
                            // Viscoelasticity factor
                            const double KXY = dt * GXY / (dt * GXY + ETAXY);
                            // Numerical viscosity
                            ETAXY = ETAXY * KXY;
                            // Numerical stresses
                            const double SZY1 = SZY0(i - 1, j) * (1 - KXY);
                            const double SZY2 = SZY0(i, j) * (1 - KXY);
                            const double SZX1 = SZX0(i, j - 1) * (1 - KXY);
                            const double SZX2 = SZX0(i, j) * (1 - KXY);
                            // Left part
                            Trip.insert(Trip.end(), {Trp(kz, kz, -2 * (ETAXY / dx2 + ETAXY / dy2) - RHOX(i, j) / dt), Trp(kz, kz - Ny1 * Num_var, ETAXY / dx2), Trp(kz, kz + Ny1 * Num_var, ETAXY / dx2),
                                                     Trp(kz, kz - Num_var, ETAXY / dy2), Trp(kz, kz + Num_var, ETAXY / dy2)}); // vzs3, vzs3, vzs3, vzs2, vzs4
                            // Right part
                            R(kz) = -RHOX(i, j) * (VZ0(i, j) / dt) - (SZX2 - SZX1) / dx - (SZY2 - SZY1) / dy;
                        }
                    }
                }
            }

            SparseMatrix<double> L(N, N); // Matrix of coefficients in the left part
            L.setFromTriplets(Trip.begin(), Trip.end()); // Build Sparse Matrix
            // 6) Solving matrix
            L.makeCompressed();
            solver.analyzePattern(L);
            solver.factorize(L);
            S = solver.solve(R);

            // 7) Reload solution
            // pfavr = 0;
            // pcount = 0;
            
            // slightly slower in parallel
            // #pragma omp parallel for collapse(2)
            for (int j = 0; j < Nx1; j++) {
                for (int i = 0; i < Ny1; i++) {
                    // Global indexes for vx, vy, P
                    const int kp = (j * Ny1 + i) * Num_var;
                    // Reload solution
                    pt(i, j) = S(kp) * ptscale;
                    vxs(i, j) = S(kp + 1);
                    vys(i, j) = S(kp + 2);
                    pf(i, j) = S(kp + 3) * pfscale;
                    vxD(i, j) = S(kp + 4);
                    vyD(i, j) = S(kp + 5);
                    if (antiplane) {
                        vzs(i, j) = S(kp + 6);
                    }
                }
            }

            Vmax = VSLIPB.maxCoeff();
            
            /*
            if (dt > 1e4 && Vmax < 1e-7) {
                double avgpt = pt.sum() / (double)(pt.rows() * pt.cols()); //calculate average total pressure
                double diffpt = (PCONF + PTFDIFF) - avgpt;
                pt += MatXd::Constant(Ny1, Nx1, diffpt);
            }
            */
            
            // Velocity change
            if (antiplane) {
                vzs.col(Nx) = vzs.col(Nx - 1);
                DVZ0 = vzs - VZ0;
            }

            DVX0 = vxs - VX0;
            DVY0 = vys - VY0;
            
            // Define timestep
            bool yn = false;

            // Plastic iterations
            // Compute strain rate, stress and stress change
            // Process internal basic nodes
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    EXY(i, j) = .5 * ((vxs(i + 1, j) - vxs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dx);
                    const double KXY = divplus(dt * GGG(i, j), ETA(i, j));
                    SXY(i, j) = 2 * ETA(i, j) * EXY(i, j) * KXY + SXY0(i, j) * (1 - KXY);

                    if (antiplane) {
                        EZX(i, j) = .5 * ((vzs(i, j + 1) - vzs(i, j)) / dx + (vxs(i + 1, j) - vxs(i, j)) / dx);
                        EZY(i, j) = .5 * ((vzs(i + 1, j) - vzs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dy);
                        SZX(i, j) = 2 * ETA(i, j) * EZX(i, j) * KXY + SZX0(i, j) * (1 - KXY);
                        SZY(i, j) = 2 * ETA(i, j) * EZY(i, j) * KXY + SZY0(i, j) * (1 - KXY);
                    }
                }
            }

            // Process pressure cells
            for (int i = 1; i < Ny; i++) {
                for (int j = 1; j < Nx; j++) {
                    // EXX, SXX
                    EXX(i, j) = (2 * (vxs(i, j) - vxs(i, j - 1)) / dx - (vys(i, j) - vys(i - 1, j)) / dy) / 3.;
                    EYY(i, j) = (2 * (vys(i, j) - vys(i - 1, j)) / dy - (vxs(i, j) - vxs(i, j - 1)) / dx) / 3.;
                    const double KXX = divplus(dt * GGGP(i, j), ETAP(i, j));
                    SXX(i, j) = 2 * ETAP(i, j) * EXX(i, j) * KXX + SXX0(i, j) * (1 - KXX);
                    SYY(i, j) = 2 * ETAP(i, j) * EYY(i, j) * KXX + SYY0(i, j) * (1 - KXX);
                }
            }
            
            // External P - nodes: symmetry
            for (auto i : {pt, pf, EXX, SXX, SXX0, EYY, SYY, SYY0, ETAP, ETAB, GGGP, GGGB}) {
                copy_bounds(i);
            }

            if (antiplane) {
                SZX.col(Nx - 1) = SZX.col(Nx - 2);
                EZX.col(Nx - 1) = EZX.col(Nx - 2);

                // Compute stress and strain rate invariants and dissipation
                // Process pressure cells
                for (int i = 1; i < Ny; i++) {
                    for (int j = 1; j < Nx; j++) {
                        // EXY term is averaged from four surrounding basic nodes
                        EII(i, j) = sqrt(.5 * (pow(EXX(i, j), 2) + pow(EYY(i, j), 2)) + (square_block(EXY.block(i - 1, j - 1, 2, 2)) + square_block(EZX.block(i - 1, j - 1, 2, 2)) + square_block(EZY.block(i - 1, j - 1, 2, 2))) / 4.);
                        
                        // Second strain rate invariant SII
                        // SXY term is averaged from four surrounding basic nodes
                        SII(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + pow(SYY(i, j), 2)) + (square_block(SXY.block(i - 1, j - 1, 2, 2)) + square_block(SZX.block(i - 1, j - 1, 2, 2)) + square_block(SZY.block(i - 1, j - 1, 2, 2))) / 4.);
                    }
                }
            }

            // Update viscosity for yielding
            // dt0 = dt; dt = dt * 1.1;
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
                    // if (y(i)  >= upper_block && y(i)  <= lower_block) {
                    if (i == line_fault) {
                        for (int j = 0; j < Nx; j++) {
                            // reducing matrix calls
                            const double arsf_temp = ARSF(i, j), brsf_temp = BRSF(i, j), lrsf_temp = LRSF(i, j), fric_temp = FRIC(i, j), eta0_temp = ETA0(i, j);

                            const double sxx_temp = SXX.block(i, j, 2, 2).sum();
                            const double syy_temp = SYY.block(i, j, 2, 2).sum();

                            // SXX, pt are averaged from four surrounding pressure nodes
                            if (antiplane) {
                                SIIB(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + pow(SYY(i, j), 2)) + pow(SXY(i, j), 2) + pow(SZX(i, j), 2) + pow(SZY(i, j), 2));                        
                            } else {
                                SIIB(i, j) = sqrt(pow(SXY(i, j), 2) + (pow(sxx_temp, 2) + pow(syy_temp, 2) + pow(sxx_temp + syy_temp, 2))/ 32.);
                            }
                            // Computing "elastic" stress invariant
                            const double siiel = SIIB(i, j) / (ETA(i, j) / (GGG(i, j) * dt + ETA(i, j)));
                            
                            // Compute old viscoplastic slip rate
                            // Compute PEFF
                            const double prB_temp = (pt.block(i, j, 2, 2).sum() - pf.block(i, j, 2, 2).sum()) / 4.;
                            double prB = prB_temp;
                            
                            if (prB < 1e3) {
                                prB = 1e3;
                            }
                            // Compute old power law strain rate
                            double SIIB1 = SIIB(i, j);
                            
                            // Compute slip velocity for current stress invariant and state
                            double EIISLIP = V0 * sinh(max(SIIB1, 0.) / arsf_temp / prB) * exp(-(brsf_temp * OM(i, j) + fric_temp) / arsf_temp) / dx;
                            
                            // Compute new ETAVP
                            double ETAPL = SIIB1 / 2. / EIISLIP;
                            double ETAVP = 1. / (1. / eta0_temp + 1. / ETAPL);
                            // Compute new stress invariant
                            double SIIB2 = siiel * ETAVP / (GGG(i, j) * dt + ETAVP);
                            const double DSIIB1 = SIIB2 - SIIB1;
                            
                            // Compute slip velocity for current stress invariant and state
                            double V = 2 * V0 * sinh(max(SIIB2, 0.) / arsf_temp / prB) * exp(-(brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                            
                            EIISLIP = V / dx / 2.;
                            
                            // Compute new ETAVP
                            ETAPL = SIIB2 / 2. / EIISLIP;
                            ETAVP = 1. / (1. / eta0_temp + 1. / ETAPL);
                            // Compute new stress invariant
                            const double DSIIB2 = siiel * ETAVP / (GGG(i, j) * dt + ETAVP) - SIIB2;
                            double SIIB4 = 0.;
                            
                            if ((DSIIB1 >= 0 && DSIIB2 <= 0) || (DSIIB1 <= 0 && DSIIB2 >= 0)) {
                                double DSIIB = 1e9;
                                
                                while(abs(DSIIB) > 1e-3) {
                                    SIIB4 = (SIIB1 + SIIB2) / 2.;
                                    
                                    //Compute slip velocity for current stress invariant and state
                                    V = 2 * V0 * sinh(max((SIIB4), 0.) / arsf_temp / prB) * exp( - (brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                                    
                                    EIISLIP = V / dx / 2.;
                                    
                                    // Compute new ETAVP
                                    ETAPL = SIIB4 / 2. / EIISLIP;
                                    ETAVP = 1. / (1. / eta0_temp + 1. / ETAPL);
                                    // Compute new stress invariant
                                    DSIIB = siiel * ETAVP / (GGG(i, j) * dt + ETAVP) - SIIB4;
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
                            
                            // Compute yielding stress
                            const double syield = max(syieldmin, prB_temp * arsf_temp * asinh(V / 2. / V0 * exp((brsf_temp * OM5(i, j) + fric_temp) / arsf_temp)));
                            
                            // Compute visco - plastic viscosity
                            double etapl = eta0_temp * syield / (eta0_temp * V + syield);
                            
                            // Save syield
                            SigmaY(i, j) = syield;
                            VSLIPB(i, j) = V;
                            SII_fault(i, j) = SIIB4;

                            // reduces calls on matrix
                            const double g_temp = 2 * GGG(i, j); // reduces calls on matrix
                            
                            // "/ BETASOLID" -> Timestep criterion, Lapusta et al., 2000; Lapusta and Liu, 2009
                            const double k = g_temp / (pi * (1 - ((3. / BETASOLID - g_temp) / (6. / BETASOLID + g_temp))) * dx);
                            double dTETAmax;
                            if ((pow((k * lrsf_temp / prB - brsf_temp) / arsf_temp - 1, 2) / 4. - k * lrsf_temp / arsf_temp / prB) < 0) {
                                dTETAmax = min(1. - (brsf_temp - arsf_temp) * prB / (k * lrsf_temp), .2);
                            } else {
                                dTETAmax = min(arsf_temp * prB / (k * lrsf_temp - (brsf_temp - arsf_temp) * prB), .2);
                            }
                            dtlapusta = min(dtlapusta, dTETAmax * lrsf_temp / V);
                            
                            // Count old yelding nodes
                            bool ynn = false;
                            if (YNY0(i, j) > 0) {
                                ynn = true;
                                DSY(i, j) = SIIB(i, j) - syield;
                                ddd += pow(DSY(i, j), 2);
                                ynpl++;
                            }

                            // Update viscosity
                            const double A = syield / siiel;
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
            const double D_iter = DSYLSQ(iterstep);
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
            double maxvxs = max(vxs.maxCoeff(), abs(vxs.minCoeff()));
            double maxvys = max(vys.maxCoeff(), abs(vys.minCoeff()));
            double maxvzs;
            if (antiplane) {
                maxvzs = max(vzs.maxCoeff(), abs(vzs.minCoeff()));
                maxvxy = sqrt(pow(vxs.maxCoeff() - vxs.minCoeff(), 2) + pow(vys.maxCoeff() - vys.minCoeff(), 2) + pow(vzs.maxCoeff() - vzs.minCoeff(), 2));
            } else {
                maxvxy = sqrt(pow(vxs.maxCoeff() - vxs.minCoeff(), 2) + pow(vys.maxCoeff() - vys.minCoeff(), 2));
            }
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
            if (antiplane) {
                dtz = dt;
                if(dt > dy * stpmaxcur / maxvzs) {
                    dtz = dx / dtkoefv * stpmaxcur / maxvzs;
                    yn = true;
                }
                maxvzs = 0;
            }
            maxvxs = 0;
            maxvys = 0;

            for (int i = 0; i < Ny1; i++) {
                for (int j = 0; j < Nx1; j++) {
                    if (yvx(i) >= upper_block && yvx(i) <= lower_block) {
                        maxvxs = max(maxvxs, abs(vxs(i, j)));
                    }
                    if (yvy(i) >= upper_block && yvy(i) <= lower_block) {
                        maxvys = max(maxvys, abs(vys(i, j)));
                    }
                    if (antiplane && yvy(i) >= upper_block && yvy(i) <= lower_block) {
                        maxvzs = max(maxvzs, abs(vzs(i, j)));
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
            if (antiplane && dt > dy * stpmaxcur / maxvzs) {
                dtz = dy / dtkoefv * stpmaxcur / maxvzs;
                yn = true;
            }

            dtslip = 1e30;
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    if (VSLIPB(i, j) > 0) {
                        dtslip = min(dtslip, dx * stpmax / VSLIPB(i, j));
                    }
                }
            }
            
            if (ynpl > 0 && dtslip < dt) {
                yn = true;
                dtslip = dtslip / dtkoefv;
            }
            
            
            // Chose minimal timestep
            if (yn && dt > dtmin) {
                const double dtold = dt;
                dt = max(min(min(dtx, dty), min(min(dtpl, dtslip), dtlapusta)), dtmin);
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
                        ETA(i, j) = max(min(ETA5(i, j), ETA0(i, j)), etamin);
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
        dt = min(min(dtx, dty), min(dt, dtlapusta));
        
        // Compute strain rate, stress and stress change
        for (auto i : {ESP, EXY, SXY, EXX, SXX, EYY, SYY, EII, EIIVP, SII, DSII, DIS}) {
            i.setZero();
        }

        // Process internal basic nodes
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                // ESP = .5 *(dVy / dx - dVx / dy), EXY, SXY
                ESP(i, j) = .5 * ((vys(i, j + 1) - vys(i, j)) / dx - (vxs(i + 1, j) - vxs(i, j)) / dy);
                EXY(i, j) = .5 * ((vxs(i + 1, j) - vxs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dx);
                const double KXY = divplus(dt * GGG(i, j), ETA(i, j));
                SXY(i, j) = 2 * ETA(i, j) * EXY(i, j) * KXY + SXY0(i, j) * (1 - KXY);

                if (antiplane) {
                    EZX(i, j) = .5 *((vzs(i, j + 1) - vzs(i, j)) / dx + (vxs(i + 1, j) - vxs(i, j)) / dx);
                    EZY(i, j) = .5 *((vzs(i + 1, j) - vzs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dy);
                    SZX(i, j) = 2 * ETA(i, j) * EZX(i, j) * KXY + SZX0(i, j) * (1 - KXY);
                    SZY(i, j) = 2 * ETA(i, j) * EZY(i, j) * KXY + SZY0(i, j) * (1 - KXY);
                }
            }
        }
        
        // Process pressure cells
        // #pragma omp parallel for collapse(2) // about 2x faster with n = 4 but breaks the simulation
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx; j++) {
                if (antiplane) {
                    // EXY term is averaged from four surrounding basic nodes
                    EII(i, j) = sqrt(pow(EXX(i, j), 2) + (square_block(EXY.block(i - 1, j - 1, 2, 2)) + square_block(EZX.block(i - 1, j - 1, 2, 2)) + square_block(EZY.block(i - 1, j - 1, 2, 2))) / 4.);
                    
                    // Second strain rate invariant SII
                    // SXY term is averaged from four surrounding basic nodes
                    SII(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + pow(SYY(i, j), 2)) + (square_block(SXY.block(i - 1, j - 1, 2, 2)) + square_block(SZX.block(i - 1, j - 1, 2, 2)) + square_block(SZY.block(i - 1, j - 1, 2, 2))) / 4.);
                    
                    // Dissipation
                    const double DISXY = pow(SXY(i, j), 2) / ETA(i, j) + pow(SXY(i - 1, j), 2) / ETA(i - 1, j) + pow(SXY(i, j - 1), 2) / ETA(i, j - 1) + pow(SXY(i - 1, j - 1), 2) / ETA(i - 1, j - 1);
                    const double DISZX = pow(SZX(i, j), 2) / ETA(i, j) + pow(SZX(i - 1, j), 2) / ETA(i - 1, j) + pow(SZX(i, j - 1), 2) / ETA(i, j - 1) + pow(SZX(i - 1, j - 1), 2) / ETA(i - 1, j - 1);
                    const double DISZY = pow(SZY(i, j), 2) / ETA(i, j) + pow(SZY(i - 1, j), 2) / ETA(i - 1, j) + pow(SZY(i, j - 1), 2) / ETA(i, j - 1) + pow(SZY(i - 1, j - 1), 2) / ETA(i - 1, j - 1);
                    DIS(i, j) = (pow(SXX(i, j), 2) / ETAP(i, j) + pow(SYY(i, j), 2) / ETAP(i, j) + (DISXY + DISZX + DISZY) / 2.) / 2.;
                } else {
                    // EXX, SXX
                    EXX(i, j) = (2 * (vxs(i, j) - vxs(i, j - 1)) / dx - (vys(i, j) - vys(i - 1, j)) / dy) / 3.;
                    EYY(i, j) = (2 * (vys(i, j) - vys(i - 1, j)) / dy - (vxs(i, j) - vxs(i, j - 1)) / dx) / 3.;
                    const double KXX = divplus(dt * GGGP(i, j), ETAP(i, j));
                    SXX(i, j) = 2 * ETAP(i, j) * EXX(i, j) * KXX + SXX0(i, j) * (1 - KXX);
                    SYY(i, j) = 2 * ETAP(i, j) * EYY(i, j) * KXX + SYY0(i, j) * (1 - KXX);
            
                    // Compute stress and strain rate invariants and dissipation
                    Matrix2d temp;
                    temp << SXY(i - 1, j - 1) / ETA(i - 1, j - 1), SXY(i - 1, j) / ETA(i - 1, j), SXY(i, j - 1) / ETA(i, j - 1), SXY(i, j) / ETA(i, j);

                    // EXY term is averaged from four surrounding basic nodes
                    EII(i, j) = sqrt(pow(EXX(i, j), 2) + square_block(EXY.block(i - 1, j - 1, 2, 2)) / 4.);
                    EIIVP(i, j) = sqrt(.5 * (pow(SXX(i, j) / (2 * ETAP(i, j)), 2) + pow(SYY(i, j) / (2 * ETAP(i, j)), 2)) + square_block(temp) / 16.);
                    // Second strain rate invariant SII
                    // SXY term is averaged from four surrounding basic nodes
                    SII(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + pow(SYY(i, j), 2)) + square_block(SXY.block(i - 1, j - 1, 2, 2)) / 4.);
                    
                    // Dissipation
                    double DISXY = (pow(SXY(i, j), 2) / ETA(i, j) + pow(SXY(i - 1, j), 2) / ETA(i - 1, j) + pow(SXY(i, j - 1), 2) / ETA(i, j - 1) + pow(SXY(i - 1, j - 1), 2) / ETA(i - 1, j - 1)) / 2.;
                    DIS(i, j) = pow(SXX(i, j), 2) / (2 * ETAP(i, j)) + pow(SYY(i, j), 2) / (2 * ETAP(i, j)) + DISXY;
                }
                
                double PT0_ave, PF0_ave;
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
                const double BETADRAINED = (1 / GGGB(i, j) + BETASOLID) / (1 - POR(i, j));
                // Biott - Willis koefficient
                const double KBW = 1 - BETASOLID / BETADRAINED;
                EL_DECOM(i, j) = BETADRAINED * (pt_ave(i, j) - PT0_ave - KBW * pf_ave(i, j) + KBW * PF0_ave) / dt;
            }
        }

        // Move markers by nodal velocity field
        // #pragma omp parallel for // 3x slower for n = 4
        for (int m = 0; m < 0; m++) {
        // for (int m = 0; m < marknum; m++) {
            Vector4d vxm = Vector4d::Zero(), vym = Vector4d::Zero(), spm = Vector4d::Zero(); // Runge-Kutta velocity, spin array
            // Save marker position
            const double xold = xm(m);
            const double yold = ym(m);
            for (int rk = 0; rk < 4; rk++) {
                // vx - velocity interpolation
                // [i, j] -------- [i, j + 1]
                //   |                |
                //   |    o m         |
                //   |                |
                // [i + 1, j] ------- [i + 1, j + 1]
                // Indexes and distances                
                int j = check_bounds(fix(xm(m) / dx), Nx);
                int i = check_bounds(fix((ym(m) + dy / 2.) / dy), Ny);

                //Distances
                double dxm = (xm(m) - xvx(j)) / dx;
                double dym = (ym(m) - yvx(i)) / dy;
                // Weights
                // Interpolation
                vxm(rk) = vxs(i, j) * (1 - dxm) * (1 - dym) + vxs(i + 1, j) * (1 - dxm) * dym + vxs(i, j + 1) * dxm * (1 - dym) + vxs(i + 1, j + 1) * dxm * dym;
                
                // vy - velocity interpolation
                // [i, j] -------- [i, j + 1]
                //   |                |
                //   |    o m         |
                //   |                |
                // [i + 1, j] ------- [i + 1, j + 1]
                // Indexes and distances
                j = check_bounds(fix((xm(m) + dx / 2.) / dx), Nx);
                i = check_bounds(fix(ym(m) / dy), Ny);

                //Distances
                dxm = (xm(m) - xvy(j)) / dx;
                dym = (ym(m) - yvy(i)) / dy;
                // Weights
                // Interpolation
                vym(rk) = vys(i, j) * (1 - dxm) * (1 - dym) + vys(i + 1, j) * (1 - dxm) * dym + vys(i, j + 1) * dxm * (1 - dym) + vys(i + 1, j + 1) * dxm * dym;
                
                // ESP = .5 *(dVy / dx - dVx / dy) interpolation
                // [i, j] -------- [i, j + 1]
                //   |                |
                //   |    o m         |
                //   |                |
                // [i + 1, j] ------- [i + 1, j + 1]
                // Indexes and distances
                j = check_bounds(fix((xm(m)) / dx), Nx);
                i = check_bounds(fix((ym(m)) / dy), Ny);
                
                //Distances
                dxm = (xm(m) - x(j)) / dx;
                dym = (ym(m) - y(i)) / dy;
                // Weights
                // Interpolation ESP = .5 *(dVy / dx - dVx / dy) for the marker
                spm(rk) = ESP(i, j) * (1 - dxm) * (1 - dym) + ESP(i + 1, j) * (1 - dxm) * dym + ESP(i, j + 1) * dxm * (1 - dym) + ESP(i + 1, j + 1) * dxm * dym;
                
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
            // Rotate stress on marker according to its spin
            // Compute amount of rotation from spin rate:
            // Espin = .5 *(dvy / dx - dvx / dy) i.e. positive for clockwise rotation
            // (when x axis is directed rightward and y axis is directed downward)
            const double dspeff = spm(0) * dt;
            // Save old stresses
            const double msxxold = sxxm(m);
            const double msyyold = syym(m);
            const double msxyold = sxym(m);
            sxym(m) = .5 * (msxxold - msyyold) * sin(2 * dspeff) + msxyold * cos(2 * dspeff);
            sxxm(m) = msxxold * pow(cos(dspeff), 2) + msyyold * pow(sin(dspeff), 2) - msxyold * sin(2 * dspeff);
            syym(m) = msxxold * pow(sin(dspeff), 2) + msyyold * pow(cos(dspeff), 2) + msxyold * sin(2 * dspeff);
            
            // Move markers
            xm(m) = xold + (vxm(0) + 2 * vxm(1) + 2 * vxm(2) + vxm(3)) / 6. * dt;
            ym(m) = yold + (vym(0) + 2 * vym(1) + 2 * vym(2) + vym(3)) / 6. * dt;
            
            // Recycling
            if (xm(m) < 0) {
                xm(m) = xm(m) + xsize;
            }
            if (xm(m) > xsize) {
                xm(m) = xm(m) - xsize;
            }
        }

        const int temp = timestep - 1;
        
        // Update timesum
        timesum += dt;
        timesumcur(temp) = timesum;
        dtcur(temp) = dt;
        
        maxvxsmod(temp) = -1e30;
        minvxsmod(temp) = 1e30;
        maxvysmod(temp) = -1e30;
        minvysmod(temp) = 1e30;
        maxvzsmod(temp) = -1e30;
        minvzsmod(temp) = 1e30;

        VX0 = vxs;
        VY0 = vys;
        if (antiplane) {
            VZ0 = vzs;
        }
        
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                // Vx
                if (RHOX(i, j) > 2000 && i != 0) {
                    maxvxsmod(temp) = max(maxvxsmod(temp), vxs(i, j));
                    minvxsmod(temp) = min(minvxsmod(temp), vxs(i, j));
                }
                // Vy
                if (RHOY(i, j) > 2000 && j != 0) {
                    maxvysmod(temp) = max(maxvysmod(temp), vys(i, j));
                    minvysmod(temp) = min(minvysmod(temp), vys(i, j));
                }
                // Vz
                if (antiplane && RHOY(i, j) > 2000) {
                    maxvzsmod(temp) = max(maxvzsmod(temp), vzs(i, j));
                    minvzsmod(temp) = min(minvzsmod(temp), vzs(i, j));
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
        if (antiplane) {
            // Update SZX0
            SZX0 = SZX;
            // Update SZY0
            SZY0 = SZY;
        }
        // Update PTF0
        PTF0 = (pt - pf);
        PT0 = pt;
        PF0 = pf;

        ///////////////////////////////////////////////////////////////////////////////////////// 
        // output
        ///////////////////////////////////////////////////////////////////////////////////////// 
        
        cout << "====================================" << endl;
        cout << "total time:        " << timesum << " sec" << endl;
        cout << "time step:         " << dt << " sec" << endl;
        cout << "Vslip max:         " << Vmax << endl;
        cout << "iter - iterations:   " << iterstep + 1 << endl;
        cout << "global - iterations: " << ynlast + 1 << endl;
        
        if (timestep == 1) {
            ofstream out_fault;
            out_fault.open("x_fault.txt");
            out_fault << x << endl;
            out_fault.close();

            ofstream out_rsf;
            out_rsf.open("rsf_fault.txt", ios_base::app | ios_base::out);
            out_rsf << ARSF.row(line_fault) << "\n\n" << BRSF.row(line_fault) << "\n\n" << LRSF.row(line_fault) << endl;
            out_rsf.close();
        }

        // ========== save slip rate
        ofstream out_Vslip("EVO_Vslip.txt", ios_base::app | ios_base::out);
        out_Vslip << timesum << "    " << dt << "    " << VSLIPB.row(line_fault) << endl;
        out_Vslip.close();
        
        // ========== save viscosity
        ofstream out_viscosity("EVO_viscosity.txt", ios_base::app | ios_base::out);
        out_viscosity << timesum << "    " << dt << "    " << ETA.row(line_fault) << endl;
        out_viscosity.close();
        
        // ========== save fluid 
        ofstream out_press_flu("EVO_press_flu.txt", ios_base::app | ios_base::out);
        out_press_flu << timesum << "    " << dt << "    " << pf.row(line_fault) << endl;
        out_press_flu.close();
        
        // ========== save effective pressure
        ofstream out_press_eff("EVO_press_eff.txt", ios_base::app | ios_base::out);
        MatXd P_diff = pt - pf;
        out_press_eff << timesum << "    " << dt << "    " << P_diff.row(line_fault) << endl;
        out_press_eff.close();
        
        // ========== save SigmaY
        ofstream out_SigmaY("EVO_SigmaY.txt", ios_base::app | ios_base::out);
        out_SigmaY << timesum << "    " << dt << "    " << SigmaY.row(line_fault) << endl;
        out_SigmaY.close();
        
        // ========== save SII
        ofstream out_Sii("EVO_Sii.txt", ios_base::app | ios_base::out);
        out_Sii << timesum << "    " << dt << "    " << SII_fault.row(line_fault) << endl;
        out_Sii.close();
        
        // ========== save Theta
        ofstream out_Theta("EVO_Theta.txt", ios_base::app | ios_base::out);
        out_Theta << timesum << "    " << dt << "    " << OM.row(line_fault) << endl;
        out_Theta.close();
        
        // ========== save viscous compaction
        ofstream out_Visc_comp("EVO_Visc_comp.txt", ios_base::app | ios_base::out);
        out_Visc_comp << timesum << "    " << dt << "    " << VIS_COMP.row(line_fault) << endl;
        out_Visc_comp.close();
        
        // ========== save elastic compaction
        ofstream out_Elast_comp("EVO_Elast_comp.txt", ios_base::app | ios_base::out);
        out_Elast_comp << timesum << "    " << dt << "    " << EL_DECOM.row(line_fault) << endl;
        out_Elast_comp.close();

        // ========== save vx Darcy
        ofstream out_EVO_vxD("EVO_vxD.txt", ios_base::app | ios_base::out);
        out_EVO_vxD << timesum << "    " << dt << "    " << vxD.row(line_fault) << endl;
        out_EVO_vxD.close();
        
        // ========== save time, dt, vmax
        ofstream out_data("EVO_data.txt", ios_base::app | ios_base::out);
        out_data << setw(20) << timesum << setw(20) << dt << setw(20) << Vmax << setw(20) << ynlast + 1 << setw(20) << iterstep + 1 << endl;
        out_data.close();
        
        if (timestep % savestep == 0) {
            ofstream out_file;
            out_file.open("file.txt");
            out_file << timestep << endl;
            out_file.close();

            const string save_file_name = nname + to_string(timestep) + ".h5";
            const string group_matrix = "Matrix";
            const string group_vector = "Vector";
            const string group_values = "Value";

            create_file(save_file_name);
            for (auto i : {group_matrix,  group_vector, group_values}){
                add_group(save_file_name, i);
            }

            // save matrices and vectors to restart the simulation from a certain timestep.
            // !!! Always change both the string array and loop order !!!
            const hsize_t dims1[2] = {Ny, Nx};
            const string matrix_names[24] = {"SIIB", "OM0", "OM", "ARSF", "BRSF", "LRSF", "RHO", "ETA0", "ETA1", "ETA5", "ETA00", "IETAPLB", "SXY0", "YNY0", "KKK",
                                       "GGG", "COHC", "COHT", "FRIC", "FRIT", "DILC", "TTT", "EIIB", "VSLIPB"}; // {"names"} has to be the same as in *matrix
            int j = 0;
            for (auto i : {SIIB, OM0, OM, ARSF, BRSF, LRSF, RHO, ETA0, ETA1, ETA5, ETA00, IETAPLB, SXY0, YNY0, KKK, GGG, COHC, COHT, FRIC, FRIT, DILC, TTT, EIIB, VSLIPB}) { // {names} *matrix
                add_matrix(save_file_name, group_matrix, i, matrix_names[j], dims1);
                j++;
            }

            const hsize_t dims2[2] = {Ny1, Nx1};
            const string matrix_names_plus[30] = {"pt", "vxs", "vys", "pf", "vxD", "vyD", "ETAB", "ETAB0", "ETAP", "ETAP0", "POR", "GGGP", "GGGB", "PTF0", "PT0", "PF0",
                                            "SXX0", "SYY0", "RHOX", "RHOFX", "ETADX", "PORX", "VX0", "VXF0", "RHOY", "RHOFY", "ETADY", "PORY", "VY0", "VYF0"}; // {"names"} has to be the same as in *matrix_plus
            j = 0;
            for (auto i : {pt, vxs, vys, pf, vxD, vyD, ETAB, ETAB0, ETAP, ETAP0, POR, GGGP, GGGB, PTF0, PT0, PF0, SXX0, SYY0, RHOX, RHOFX, ETADX,
                           PORX, VX0, VXF0, RHOY, RHOFY, ETADY, PORY, VY0, VYF0}) { // {names} *matrix_plus
                add_matrix(save_file_name, group_matrix, i, matrix_names_plus[j], dims2);
                j++;
            }
            
            const hsize_t dim1[1] = {num_timesteps};
            const string vector_names[6] = {"timesumcur", "dtcur", "maxvxsmod", "minvxsmod", "maxvysmod", "minvysmod"}; // {"names"} has to be the same as in *vec
            j = 0;
            for (auto i : {timesumcur, dtcur, maxvxsmod, minvxsmod, maxvysmod, minvysmod}) { // {names} *vec
                add_vector(save_file_name, group_vector, i, vector_names[j], dim1);
                j++;
            }

            const hsize_t dim2[1] = {marknum};
            const string vector_names_marker[5] = {"xm", "ym", "sxxm", "syym", "sxym"}; // {"names"} has to be the same as in *vec2
            j = 0;
            for (auto i : {xm, ym, sxxm, syym, sxym}) { // {names} *vec2
                add_vector(save_file_name, group_vector, i, vector_names_marker[j], dim2);
                j++;
            }

            const hsize_t dim3[1] = {9};
            VecXd temp(9);
            temp << timesum, dt00, dtx, dty, dtlapusta, Vmax, maxvxy, dt, yndtdecrease;
            add_vector(save_file_name, group_values, temp, "values", dim3);

            if (antiplane) {
                const string group_antiplane = "Antiplane";
                add_group(save_file_name, group_antiplane);
                add_matrix(save_file_name, group_antiplane, SZX0, "SZX0", dims1);
                add_matrix(save_file_name, group_antiplane, SZY0, "SZY0", dims1);

                add_matrix(save_file_name, group_antiplane, VZ0, "VZ0", dims2);

                add_vector(save_file_name, group_antiplane, maxvzsmod, "maxvzsmod", dim1);
                add_vector(save_file_name, group_antiplane, minvzsmod, "minvzsmod", dim1);
            }
        }
    }
    return 0;
}