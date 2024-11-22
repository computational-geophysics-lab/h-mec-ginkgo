#ifndef GLOBALVAR_H
#define GLOBALVAR_H

#include <string>
#include <eigen3/Eigen/Eigen>

using namespace std;


// Basic nodes
extern Eigen::MatrixXd OM0;  // Old state parameter
extern Eigen::MatrixXd OM;   // State parameter
extern Eigen::MatrixXd OM5;
extern Eigen::MatrixXd ARSF; // a - parameter of RSF
extern Eigen::MatrixXd BRSF; // b - parameter of RSF
extern Eigen::MatrixXd LRSF; // L - parameter of RSF

// Unknown parameters
extern Eigen::MatrixXd pt;  // Total pressure
extern Eigen::MatrixXd vxs; // Solid vx - velocity
extern Eigen::MatrixXd vys; // Solid vy - velocity
extern Eigen::MatrixXd vzs; // Solid vz - velocity
extern Eigen::MatrixXd pf;  // Fluid pressure
extern Eigen::MatrixXd vxD; // Darsi vx - velocity
extern Eigen::MatrixXd vyD; // Darsi vy - velocity

// Nodal matrices
// Basic nodes
extern Eigen::MatrixXd RHO, ETA, ETA0, ETA1, ETA5, ETA00, IETAPLB, SXY, SXY0, SZX0, SZY0, YNY0, KKK,
      GGG, COHC, COHT, FRIC, FRIT, DILC, TTT, EIIB, VSLIPB, SZX, SZY;

// Pressure nodes
extern Eigen::MatrixXd ETAB, ETAB0, ETAP, ETAP0, POR, GGGP, GGGB, PTF0, PT0, PF0, pt_ave,
      pf_ave, SXX, SXX0, SYY, SYY0, DILP;
// Vx nodes
extern Eigen::MatrixXd RHOX, RHOFX, ETADX, PORX, VX0, VXF0;
// Vy nodes
extern Eigen::MatrixXd RHOY, RHOFY, ETADY, PORY, VY0, VYF0;
// Vz nodes
extern Eigen::MatrixXd VZ0;

extern Eigen::MatrixXd ESP, EXY, EZX, EZY, EXX, EYY, EII, EIIVP, SII, DSII, DIS;

extern Eigen::MatrixXd EL_DECOM;   // Elastic (de)compaction
extern Eigen::MatrixXd VIS_COMP;

// Lagrangian solid markers
extern Eigen::VectorXd t_marker;         // Marker rock type
extern Eigen::VectorXd rhom;          // Density of solid
extern Eigen::VectorXd etasm;         // Standard shear viscosity of bulk
extern Eigen::VectorXd etam;                                  // Shear viscosity of bulk
extern Eigen::VectorXd cohescm;      // Cohesion for confined fracture of solid
extern Eigen::VectorXd cohestm;      // Cohesion for tensile fracture of solid
extern Eigen::VectorXd frictcm;         // friction for confined fracture of solid
extern Eigen::VectorXd dilatcm; // dilatation for confined fracture of solid
extern Eigen::VectorXd fricttm;    // friction for tensile fracture of solid
extern Eigen::VectorXd porm;                                  // Porosity of solid
extern Eigen::VectorXd kkkm;         // Standard permeability of solid
extern Eigen::VectorXd xm;                                    // Horizontal coordinates of solid markers
extern Eigen::VectorXd ym;                                    // Vertical coordinates of solid markers
extern Eigen::VectorXd sxxm;                                  // Marker SIGMAxx', Pa
extern Eigen::VectorXd syym;                                  // Marker SIGMAyy', Pa
extern Eigen::VectorXd sxym;                                  // Marker SIGMAxy', Pa

// (3) Defining global matrixes
// according to the global number of unknowns
// Sparse Matrix L is not yet defined as it will be built from a set of Triplets each step
extern Eigen::VectorXd R; // Vector of the right parts of equations
extern Eigen::VectorXd S;

// variable type declaration
extern int ynlast, iterstep;

extern double timesum;
extern double dt;
extern bool yndtdecrease;

extern double dt00, dtx, dty, dtz, dtlapusta, Vmax, maxvxy;

extern Eigen::VectorXd DSYLSQ;

extern Eigen::MatrixXd DVX0, DVY0, DVZ0, DSY, YNY, SigmaY, SII_fault, SIIB; // DVX0, DVY0, DVZ0 unused

extern Eigen::VectorXd timesumcur, dtcur;
extern Eigen::VectorXd maxvxsmod, minvxsmod, maxvysmod, minvysmod, maxvzsmod, minvzsmod;

#endif