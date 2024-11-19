#include "run_simulation.h"

using namespace std;
using namespace H5;


void run_simulation(int &timestep){

    Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;

    // declaration of matrices that are set to zero each timestep
    Eigen::MatrixXd ETA0SUM(Ny, Nx), COHCSUM(Ny, Nx), FRICSUM(Ny, Nx), DILCSUM(Ny, Nx), COHTSUM(Ny, Nx), FRITSUM(Ny, Nx), WTSUM(Ny, Nx);

    // Interpolate ETA, RHO to nodal points
    // Basic nodes
    Eigen::MatrixXd RHOSUM(Ny, Nx), ETASUM(Ny, Nx), KKKSUM(Ny, Nx), TTTSUM(Ny, Nx), SXYSUM(Ny, Nx), GGGSUM(Ny, Nx);

    //LDZ
    Eigen::MatrixXd OM0SUM(Ny, Nx);  // Old state parameter
    Eigen::MatrixXd OMSUM(Ny, Nx);   // State parameter
    Eigen::MatrixXd ARSFSUM(Ny, Nx); // a - parameter of RSF
    Eigen::MatrixXd BRSFSUM(Ny, Nx); // b - parameter of RSF
    Eigen::MatrixXd LRSFSUM(Ny, Nx); // L - parameter of RSF

    // Pressure nodes
    Eigen::MatrixXd ETAPSUM(Ny1, Nx1), ETAP0SUM(Ny1, Nx1), ETAB0SUM(Ny1, Nx1), PORSUM(Ny1, Nx1), SXXSUM(Ny1, Nx1), SYYSUM(Ny1, Nx1), GGGPSUM(Ny1, Nx1), WTPSUM(Ny1, Nx1);
    // Vx nodes
    Eigen::MatrixXd RHOXSUM(Ny1, Nx1), RHOFXSUM(Ny1, Nx1), ETADXSUM(Ny1, Nx1), PORXSUM(Ny1, Nx1), WTXSUM(Ny1, Nx1);
    // Vy nodes
    Eigen::MatrixXd RHOYSUM(Ny1, Nx1), RHOFYSUM(Ny1, Nx1), ETADYSUM(Ny1, Nx1), PORYSUM(Ny1, Nx1), WTYSUM(Ny1, Nx1);

    Eigen::MatrixXd ETA50(Ny, Nx);

    // nucleation size
    const double hstar = pi / 2. * shearmod * lrsfm(1) * brsfm(1) / pow(brsfm(1) - arsfm(1), 2) / PTFDIFF;
    // cohesive zone size
    const double coh = 9. / 32. * pi * shearmod * lrsfm(1) / brsfm(1) / PTFDIFF;
    // Print information about discretization
    cout << ">> VW width = "                    << (((TS_4 + TS_3) / 2) - ((TS_2 + TS_1) / 2)) / 1e3    << " (km)" << endl;
    cout << ">> Critical nucleation size = "    << hstar / 1e3                                          << " (km)" << endl;
    cout << ">> Cohesive zone = "               << coh                                                  << " (m)"  << endl;


    // /////////////////////////////////////////////////////////////////////////////////////// 
    // actual computations start here
    // /////////////////////////////////////////////////////////////////////////////////////// 

    for (; timestep <= num_timesteps; timestep++) {
        for (auto i : {RHOSUM, ETASUM, KKKSUM, TTTSUM, SXYSUM, GGGSUM, ETA, ETA0SUM, COHCSUM, FRICSUM, DILCSUM, COHTSUM, FRITSUM, WTSUM, OM0SUM, OMSUM, ARSFSUM, BRSFSUM, LRSFSUM, ETAPSUM, ETAP0SUM, ETAB0SUM, PORSUM, SXXSUM, SYYSUM, GGGPSUM, WTPSUM, RHOXSUM, RHOFXSUM, ETADXSUM, PORXSUM, WTXSUM, RHOYSUM, RHOFYSUM, ETADYSUM, PORYSUM, WTYSUM}) {
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
                cohescmm    = cohescm(m) * (1 - porm(m)); // * exp(-alpha * porm(m));
                cohestmm    = cohestm(m) * (1 - porm(m)); // * exp(-alpha * porm(m));
                frictcmm    = frictcm(m);
                dilatcmm    = dilatcm(m);
                fricttmm    = fricttm(m);
                etasmm0     = etasm(m);
                etamm0      = etasm(m) * exp(-alpha * porm(m));
                etamm       = etam(m);
                // total density
                rhomm       = rhom(m) * (1 - porm(m)) + rhofm(m) * porm(m);
            } else {
                cohescmm    = cohescm(m);
                cohestmm    = cohestm(m);
                frictcmm    = frictcm(m);
                dilatcmm    = 0;
                fricttmm    = fricttm(m);
                etasmm0     = etasm(m);
                etamm0      = etasm(m);
                etamm       = etam(m);
                // total density
                rhomm       = rhom(m);
            }
            
            etamm0  = enforce_bounds(etamm0, etamin, etamax); // Matrix viscosity
            etamm   = enforce_bounds(etamm, etamin, etamax); // Effective viscosity
            etadm   = etafm(m) / kkkmm; // Darsi "viscosity"
            
            // Interpolate to basic nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            int j = check_bounds(fix_towards_zero(xm(m) / dx), Nx);
            int i = check_bounds(fix_towards_zero(ym(m) / dy), Ny);
            
            double dxm = (xm(m) - x(j)) / dx;
            double dym = (ym(m) - y(i)) / dy;
            
            // merged similar computations of matlab version
            // The computations on the 4 elements are now done on a 2x2-sub-block of the matrix
            // Similarly wtm is a 2x2 matrix that stores the corresponding value
            // Sub-block += variable * wtm

            Eigen::Matrix2d wtm;
            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;

            ETASUM    += etamm * wtm;
            RHOSUM    += rhomm * wtm;
            KKKSUM    += kkkmm * wtm;
            TTTSUM    += t_marker(m) * wtm;
            SXYSUM    += sxym(m) * wtm;
            GGGSUM    += 1. / gm(m) * wtm;
            ETA0SUM   += etamm0 * wtm;
            COHTSUM   += cohestmm * wtm;
            FRITSUM   += fricttmm * wtm;
            COHCSUM   += cohescmm * wtm;
            FRICSUM   += frictcmm * wtm;
            DILCSUM   += dilatcmm * wtm;

            WTSUM     += wtm;
            
            // Interpolate to pressure nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix_towards_zero((xm(m) + dx / 2.) / dx), Nx);
            i = check_bounds(fix_towards_zero((ym(m) + dy / 2.) / dy), Ny);
            
            dxm = (xm(m) - xp(j)) / dx;
            dym = (ym(m) - yp(i)) / dy;

            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;

            ETAPSUM.block(i, j, 2, 2)   += etamm * wtm;
            PORSUM.block(i, j, 2, 2)    += porm(m) * wtm;
            SXXSUM.block(i, j, 2, 2)    += sxxm(m) * wtm;
            SYYSUM.block(i, j, 2, 2)    += syym(m) * wtm;
            GGGPSUM.block(i, j, 2, 2)   += 1. / gm(m) * wtm;
            ETAP0SUM.block(i, j, 2, 2)  += etamm0 * wtm;
            ETAB0SUM.block(i, j, 2, 2)  += etasmm0 * wtm;

            WTPSUM.block(i, j, 2, 2)    += wtm;

            // Interpolate to Vx nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix_towards_zero((xm(m)) / dx), Nx);
            i = check_bounds(fix_towards_zero((ym(m) + dy / 2.) / dy), Ny);
            
            dxm = (xm(m) - xvx(j)) / dx;
            dym = (ym(m) - yvx(i)) / dy;
            
            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;

            RHOXSUM.block(i, j, 2, 2)   += rhomm * wtm;
            ETADXSUM.block(i, j, 2, 2)  += 1. / etadm * wtm;
            RHOFXSUM.block(i, j, 2, 2)  += rhofm(m) * wtm;
            PORXSUM.block(i, j, 2, 2)   += porm(m) * wtm;
            
            WTXSUM.block(i, j, 2, 2)    += wtm;
            
            // Interpolate to Vy nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix_towards_zero((xm(m) + dx / 2.) / dx), Nx);
            i = check_bounds(fix_towards_zero((ym(m)) / dy), Ny);
        
            dxm = (xm(m) - xvy(j)) / dx;
            dym = (ym(m) - yvy(i)) / dy;

            wtm << (1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym;
            
            RHOYSUM.block(i, j, 2, 2)   += rhomm * wtm;
            ETADYSUM.block(i, j, 2, 2)  += 1. / etadm * wtm;
            RHOFYSUM.block(i, j, 2, 2)  += rhofm(m) * wtm;
            PORYSUM.block(i, j, 2, 2)   += porm(m) * wtm;

            WTYSUM.block(i, j, 2, 2)    += wtm;
        } // ends loop through markers

        // Computing ETA and RHO
        for (int i = 0; i < Ny; i++) {
            // Basic nodes
            for (int j = 0; j < Nx; j++) {
                double wtsum = WTSUM(i, j);
                if (wtsum > 0) {
                    RHO(i, j)   = RHOSUM(i, j) / wtsum;
                    KKK(i, j)   = KKKSUM(i, j) / wtsum;
                    TTT(i, j)   = TTTSUM(i, j) / wtsum;
                    GGG(i, j)   = 1. / (GGGSUM(i, j) / wtsum);
                    ETA0(i, j)  = ETA0SUM(i, j) / wtsum;
                    COHT(i, j)  = COHTSUM(i, j) / wtsum;
                    FRIT(i, j)  = FRITSUM(i, j) / wtsum;
                    COHC(i, j)  = COHCSUM(i, j) / wtsum;
                    FRIC(i, j)  = FRICSUM(i, j) / wtsum;
                    DILC(i, j)  = DILCSUM(i, j) / wtsum;
                }
            }
            // Vy nodes
            for (int j = 0; j < Nx1; j++) {
                double wtysum = WTYSUM(i, j);
                if (wtysum > 0) {
                    RHOY(i, j)  = RHOYSUM(i, j) / wtysum;
                    ETADY(i, j) = 1. / (ETADYSUM(i, j) / wtysum);
                    RHOFY(i, j) = RHOFYSUM(i, j) / wtysum;
                    PORY(i, j)  = PORYSUM(i, j) / wtysum;
                }
            }
        }
        
        for (int i = 0; i < Ny1; i++) {
            // Vx nodes
            for (int j = 0; j < Nx; j++) {
                double wtxsum = WTXSUM(i, j);
                if (wtxsum > 0) {
                    RHOX(i, j)  = RHOXSUM(i, j) / wtxsum;
                    ETADX(i, j) = 1. / (ETADXSUM(i, j) / wtxsum);
                    RHOFX(i, j) = RHOFXSUM(i, j) / wtxsum;
                    PORX(i, j)  = PORXSUM(i, j) / wtxsum;
                }
            }
            //Pressure nodes
            for (int j = 0; j < Nx1; j++) {
                double wtpsum = WTPSUM(i, j);
                if (wtpsum > 0) {
                    ETAP(i, j)  = ETAPSUM(i, j) / wtpsum;
                    POR(i, j)   = PORSUM(i, j) / wtpsum;
                    ETAB0(i, j) = ETAB0SUM(i, j) / wtpsum / POR(i, j);
                    ETAP0(i, j) = ETAP0SUM(i, j) / wtpsum;
                    GGGP(i, j)  = 1. / (GGGPSUM(i, j) / wtpsum);
                }
            }
        }
        
        // Save viscosity
        if (timestep == 1) {
            ETA1    = ETA0;
            ETA     = ETA0;
            ETA50   = ETA0;
        } else {
            ETA1    = ETA00;
            ETA     = ETA00;
            ETA50   = ETA00;
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
                        const double SXX_temp   = SXX.block(i, j, 2, 2).sum();
                        const double SYY_temp   = SYY.block(i, j, 2, 2).sum();
                        SIIB(i, j)              = sqrt(pow(SXY(i, j), 2) + (pow(SXX_temp, 2) + pow(SYY_temp, 2) + 
                                                  pow(SXX_temp  + SYY_temp, 2)) / 32.);
                        IETAPLB(i, j)           = (1. / ETA(i, j) - 1. / ETA0(i, j));
                        EIIB(i, j)              = dy_faultw * SIIB(i, j) * IETAPLB(i, j) / 2.;
                    } else {
                        EIIB(i, j)      = 0;
                        IETAPLB(i, j)   = 0;
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

            // /////////////////////////////////////////////////////////////////////////////////////// 
            // Set up L and R matrices
            // /////////////////////////////////////////////////////////////////////////////////////// 
            Eigen::SparseMatrix<double> L(N, N); // Matrix of coefficients in the left part
            // LR_setup(L, pfscale, ptscale);
            vector<Eigen::Triplet<double>> Trip; // define triplet list to build sparse matrix
            Trip.reserve(N * 10); // reserving memory space for all the entries to be stored with some reserve

            R.setZero();

            for (int j = 0; j < Nx1; j++) {
                for (int i = 0; i < Ny1; i++) {
                    // Computing global indexes for vx, vy, p
                    const int kp    = (j * Ny1 + i) * Num_var;
                    const int kx    = kp + 1;
                    const int ky    = kp + 2;
                    const int kpf   = kp + 3;
                    const int kxf   = kp + 4;
                    const int kyf   = kp + 5;
                    const int kz    = kp + 6;
                    
                    // 5a) Composing equation for vxs
                    if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (j == Nx) {
                            Trip.push_back(Eigen::Triplet<double>(kx, kx, 1));
                        }
                        
                        // Upper boundary
                        // prescribed velocity
                        if (i == 0 && j < Nx) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kx, kx, 1), 
                                Eigen::Triplet<double>(kx, kx + Num_var, 1)
                            });
                            R(kx) = 2 * bcupper;
                        }
                        
                        // Lower boundary
                        // prescribed velocity
                        if (i == Ny && j < Nx) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kx, kx, 1), 
                                Eigen::Triplet<double>(kx, kx - Num_var, 1)
                            });
                            R(kx) = 2 * bclower;
                        }
                        
                        // Left boundary
                        if (j == 0 && i > 0 && i < Ny) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kx, kx, 1), 
                                Eigen::Triplet<double>(kx, kx + Num_var * Ny1, -1)
                            });
                        }
                        
                        // Right boundary
                        if (j == Nx - 1 && i > 0 && i < Ny) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kx, kx, 1), 
                                Eigen::Triplet<double>(kx, kx - Num_var * Ny1, -1)
                            });
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
                        double ETAXY1 = ETA(i - 1, j);
                        double ETAXY2 = ETA(i, j);
                        double ETAXX1 = ETAP(i, j);
                        double ETAXX2 = ETAP(i, j + 1);
                        // Shear modulus
                        const double GXY1 = GGG(i - 1, j);
                        const double GXY2 = GGG(i, j);
                        const double GXX1 = GGGP(i, j);
                        const double GXX2 = GGGP(i, j + 1);
                        // Viscoelasticity factor
                        const double KXY1 = divplus(dt * GXY1, ETAXY1);
                        const double KXY2 = divplus(dt * GXY2, ETAXY2);
                        const double KXX1 = divplus(dt * GXX1, ETAXX1);
                        const double KXX2 = divplus(dt * GXX2, ETAXX2);
                        // Numerical viscosity
                        ETAXY1 *= KXY1;
                        ETAXY2 *= KXY2;
                        ETAXX1 *= KXX1;
                        ETAXX2 *= KXX2;
                        // Numerical stresses
                        const double SXY1 = SXY0(i - 1, j) * (1 - KXY1);
                        const double SXY2 = SXY0(i, j) * (1 - KXY2);
                        const double SXX1 = SXX0(i, j) * (1 - KXX1);
                        const double SXX2 = SXX0(i, j + 1) * (1 - KXX2);
                        // Density derivatives
                        const double dRHOdx = (RHOX(i, j + 1) - RHOX(i, j - 1)) / (2. * dx);
                        const double dRHOdy = (RHO(i, j) - RHO(i - 1, j)) / dy;
                        // Left part
                        const double temp = gx * dt * dRHOdy / 4.;
                        Trip.insert(Trip.end(), {
                            // vxs3, vxs1, vxs5, vxs2, vxs4 
                            Eigen::Triplet<double>(kx, kx, -(ETAXX1 + ETAXX2) / dx2 - (ETAXY1 + ETAXY2) / dy2 - gx * dt * dRHOdx - inertia * RHOX(i, j) / dt), 
                            Eigen::Triplet<double>(kx, kx - Ny1 * Num_var, ETAXX1 / dx2), 
                            Eigen::Triplet<double>(kx, kx + Ny1 * Num_var, ETAXX2 / dx2), 
                            Eigen::Triplet<double>(kx, kx - Num_var, ETAXY1 / dy2), 
                            Eigen::Triplet<double>(kx, kx + Num_var, ETAXY2 / dy2), 
                            // vys1, vys2, vys3, vys4 
                            Eigen::Triplet<double>(kx, ky - Num_var, (ETAXY1 - ETAXX1) / dx_dy - temp), 
                            Eigen::Triplet<double>(kx, ky, (ETAXX1 - ETAXY2) / dx_dy - temp),
                            Eigen::Triplet<double>(kx, ky - Num_var + Ny1 * Num_var, (ETAXX2 - ETAXY1) / dx_dy - temp), 
                            Eigen::Triplet<double>(kx, ky + Ny1 * Num_var, (ETAXY2 - ETAXX2) / dx_dy - temp), 
                            // Pt1', Pt2' 
                            Eigen::Triplet<double>(kx, kp, ptscale / dx), 
                            Eigen::Triplet<double>(kx, kp + Ny1 * Num_var, -ptscale / dx)
                        }); 
                        // Right part
                        R(kx) = -RHOX(i, j) * (inertia * VX0(i, j) / dt + gx) - (SXX2 - SXX1) / dx - (SXY2 - SXY1) / dy;
                    }
                    
                    // 5b) Composing equation for vys
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vys = 0
                        if (i == Ny) {
                            Trip.push_back(Eigen::Triplet<double>(ky, ky, 1));
                        }
                        
                        // Left boundary
                        // Free Slip
                        if (j == 0) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(ky, ky, 1), 
                                Eigen::Triplet<double>(ky, ky + Ny1 * Num_var, 1)
                            });
                        }
                        
                        // Right boundary
                        // Free Slip
                        if (j == Nx) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(ky, ky, 1), 
                                Eigen::Triplet<double>(ky, ky - Ny1 * Num_var, 1)
                            });
                        }
                        
                        // Upper boundary: no penetration
                        if (i == 0 && j > 0 && j < Nx) {
                            Trip.push_back(Eigen::Triplet<double>(ky, ky, 1));
                        }
                        
                        // Lower boundary: no penetration
                        if (i == Ny - 1 && j > 0 && j < Nx) {
                            Trip.push_back(Eigen::Triplet<double>(ky, ky, 1));
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
                        double ETAXY1 = ETA(i, j - 1);
                        double ETAXY2 = ETA(i, j);
                        double ETAYY1 = ETAP(i, j);
                        double ETAYY2 = ETAP(i + 1, j);
                        // Shear modulus
                        const double GXY1 = GGG(i, j - 1);
                        const double GXY2 = GGG(i, j);
                        const double GYY1 = GGGP(i, j);
                        const double GYY2 = GGGP(i + 1, j);
                        // Viscoelasticity factor
                        const double KXY1 = divplus(dt * GXY1, ETAXY1);
                        const double KXY2 = divplus(dt * GXY2, ETAXY2);
                        const double KYY1 = divplus(dt * GYY1, ETAYY1);
                        const double KYY2 = divplus(dt * GYY2, ETAYY2);
                        // Numerical viscosity
                        ETAXY1 *= KXY1;
                        ETAXY2 *= KXY2;
                        ETAYY1 *= KYY1;
                        ETAYY2 *= KYY2;
                        // Numerical stresses
                        const double SXY1 = SXY0(i, j - 1) * (1 - KXY1);
                        const double SXY2 = SXY0(i, j) * (1 - KXY2);
                        const double SYY1 = SYY0(i, j) * (1 - KYY1);
                        const double SYY2 = SYY0(i + 1, j) * (1 - KYY2);
                        // Density derivatives
                        const double dRHOdy = (RHOY(i + 1, j) - RHOY(i - 1, j)) / 2. / dy;
                        const double dRHOdx = (RHO(i, j) - RHO(i, j - 1)) / dx;
                        // Left part
                        const double temp = gy * dt * dRHOdx / 4.;
                        Trip.insert(Trip.end(), {
                            // vys3, vys1, vys5, vys2, vys4 
                            Eigen::Triplet<double>(ky, ky, -(ETAYY1 + ETAYY2) / dy2 - (ETAXY1 + ETAXY2) / dx2 - gy * dt * dRHOdy - inertia * RHOY(i, j) / dt), 
                            Eigen::Triplet<double>(ky, ky - Ny1 * Num_var, ETAXY1 / dx2),
                            Eigen::Triplet<double>(ky, ky + Ny1 * Num_var, ETAXY2 / dx2), 
                            Eigen::Triplet<double>(ky, ky - Num_var, ETAYY1 / dy2), 
                            Eigen::Triplet<double>(ky, ky + Num_var, ETAYY2 / dy2), 
                            // vxs1, vxs2, vxs3, vxs4 
                            Eigen::Triplet<double>(ky, kx - Ny1 * Num_var, (ETAXY1 - ETAYY1) / dx_dy - temp), 
                            Eigen::Triplet<double>(ky, kx + Num_var - Ny1 * Num_var, (ETAYY2 - ETAXY1) / dx_dy - temp),
                            Eigen::Triplet<double>(ky, kx, (ETAYY1 - ETAXY2) / dx_dy - temp), 
                            Eigen::Triplet<double>(ky, kx + Num_var, (ETAXY2 - ETAYY2) / dx_dy - temp), 
                            // Pt1', Pt2' 
                            Eigen::Triplet<double>(ky, kp, ptscale / dy), 
                            Eigen::Triplet<double>(ky, kp + Num_var, -ptscale / dy)
                        }); 
                        // Right part
                        R(ky) = -RHOY(i, j) * (inertia * VY0(i, j) / dt + gy) - (SYY2 - SYY1) / dy - (SXY2 - SXY1) / dx;
                    }

                    // 5c) Composing equation for Pt
                    if (i == 0 || j == 0 || i == Ny || j == Nx) { // || (i == 2 && j == 2))
                        // BC equation: 1 * Pt = 0
                        Trip.push_back(Eigen::Triplet<double>(kp, kp, 1));
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
                        Trip.insert(Trip.end(), {
                            // vxs1, vxs2 
                            Eigen::Triplet<double>(kp, kx - Ny1 * Num_var, -1. / dx), 
                            Eigen::Triplet<double>(kp, kx, 1. / dx), 
                            // vys1, vys2 
                            Eigen::Triplet<double>(kp, ky - Num_var, -1. / dy), 
                            Eigen::Triplet<double>(kp, ky, 1. / dy), 
                            // Pt, Pf
                            Eigen::Triplet<double>(kp, kp, ptscale * (1. / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED / dt)), 
                            Eigen::Triplet<double>(kp, kpf, -pfscale * (1. / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / dt))}); 
                        // Right part
                        R(kp) = BETADRAINED * (PT0(i, j) - KBW * PF0(i, j)) / dt + DILP(i, j);
                    }
                    
                    // 5d) Composing equation for vxD
                    if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (j == Nx) {
                            Trip.push_back(Eigen::Triplet<double>(kxf, kxf, 1));
                        }
                        
                        // Upper boundary: symmetry
                        if (i == 0 && j < Nx) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kxf, kxf, 1), 
                                Eigen::Triplet<double>(kxf, kxf + Num_var, -1)
                            });
                        }
                        
                        // Lower boundary: symmetry
                        if (i == Ny && j < Nx) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kxf, kxf, 1), 
                                Eigen::Triplet<double>(kxf, kxf - Num_var, -1)
                            });
                        }
                        
                        // Left boundary
                        // no penetration
                        if (j == 0) {
                            Trip.push_back(Eigen::Triplet<double>(kxf, kxf, 1));
                        }
                        
                        // Right boundary
                        // no penetration
                        if (j == Nx - 1) {
                            Trip.push_back(Eigen::Triplet<double>(kxf, kxf, 1));
                        }
                        
                    } else {
                        // Fluid X - Darsi: - ETAfluid / K * VxD - dPf / dx = -RHOf * gx + RHOf * DVxs / Dt
                        //
                        //  Pf1 --- vxD, vxs --- Pf2
                        //
                        // Left part
                        Trip.insert(Trip.end(), {
                            // vxD, vxs 
                            Eigen::Triplet<double>(kxf, kxf, -ETADX(i, j) - RHOFX(i, j) / PORX(i, j) * inertia / dt), 
                            Eigen::Triplet<double>(kxf, kx, -RHOFX(i, j) * inertia / dt), 
                            // Pf1', Pf2' 
                            Eigen::Triplet<double>(kxf, kpf, pfscale / dx), 
                            Eigen::Triplet<double>(kxf, kpf + Ny1 * Num_var, -pfscale / dx)
                        }); 
                        // Right part
                        R(kxf) = -RHOFX(i, j) * (inertia * VXF0(i, j) / dt + gx);
                    }
                    
                    // 5e) Composing equation for vyD
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (i == Ny) {
                            Trip.push_back(Eigen::Triplet<double>(kyf, kyf, 1));
                        }
                        
                        // Left boundary
                        // symmetry
                        if (j == 0 && i > 0 && i < Ny - 1) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kyf, kyf, 1), 
                                Eigen::Triplet<double>(kyf, kyf + Ny1 * Num_var, -1)
                            });
                        }
                        
                        // Right boundary
                        // symmetry
                        if (j == Nx && i > 0 && i < Ny - 1) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kyf, kyf, 1), 
                                Eigen::Triplet<double>(kyf, kyf - Ny1 * Num_var, -1)
                            });
                        }
                        
                        // Upper boundary: no penetration
                        if (i == 0) {
                            Trip.push_back(Eigen::Triplet<double>(kyf, kyf, 1));
                            R(kyf) = bcvyflower;
                        }
                        
                        // Lower boundary: no penetration
                        if (i == Ny - 1) {
                            Trip.push_back(Eigen::Triplet<double>(kyf, kyf, 1));
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
                        Trip.insert(Trip.end(), {
                            // vyD, vys 
                            Eigen::Triplet<double>(kyf, kyf, -ETADY(i, j) - RHOFY(i, j) / PORY(i, j) * inertia / dt), 
                            Eigen::Triplet<double>(kyf, ky, -RHOFY(i, j) * inertia / dt), 
                            // Pf1', Pf2'
                            Eigen::Triplet<double>(kyf, kpf, pfscale / dy), 
                            Eigen::Triplet<double>(kyf, kpf + Num_var, -pfscale / dy)
                        }); 
                        // Right part
                        R(kyf) = -RHOFY(i, j) * (inertia * VYF0(i, j) / dt + gy);
                    }
                                        
                    // 5f) Composing equation for Pf
                    if (j == 0 || j == Nx || i <= 1 || i >= Ny - 1) {
                        // BC equation: 1 * Pf = 0
                        // Real BC
                        if (i == 1 || i == Ny - 1) {
                            Trip.insert(Trip.end(), {
                                Eigen::Triplet<double>(kpf, kpf, pfscale), 
                                Eigen::Triplet<double>(kpf, kp, -ptscale)
                            });
                            R(kpf) = -PTFDIFF;
                        } else {
                            Trip.push_back(Eigen::Triplet<double>(kpf, kpf, 1));
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
                        Trip.insert(Trip.end(), {
                            // vxs1, vxs2  
                            Eigen::Triplet<double>(kpf, kxf - Ny1 * Num_var, -1. / dx), 
                            Eigen::Triplet<double>(kpf, kxf, 1. / dx), 
                            // vys1, vys2 
                            Eigen::Triplet<double>(kpf, kyf - Num_var, -1. / dy), 
                            Eigen::Triplet<double>(kpf, kyf, 1. / dy), 
                            // Pt, Pf
                            Eigen::Triplet<double>(kpf, kp, -ptscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / dt)), 
                            Eigen::Triplet<double>(kpf, kpf, pfscale * (1 / ETAB(i, j) / (1 - POR(i, j)) + BETADRAINED * KBW / KSK / dt))
                        }); 
                        // Right part
                        R(kpf) = -BETADRAINED * KBW * (PT0(i, j) - 1. / KSK * PF0(i, j)) / dt - DILP(i, j);
                    }
                    if (antiplane) {
                        // 5o) Composing equation for vzs (out-of-plane component)
                        if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                            // Ghost nodes: 1*vzs=0
                            if (j == Nx1) {
                                Trip.push_back(Eigen::Triplet<double>(kz, kz, 1));
                            }
                            
                            // Upper boundary
                            if (i == 1 && j < Nx1) {
                                Trip.insert(Trip.end(), {
                                    Eigen::Triplet<double>(kz, kz, 1), 
                                    Eigen::Triplet<double>(kz, kz + Num_var, 1)
                                });
                                R(kz) = 2 * bclower;
                            }
                            
                            // Lower boundary
                            if (i == Ny1 && j < Nx1) {
                                Trip.insert(Trip.end(), {
                                    Eigen::Triplet<double>(kz, kz, 1), 
                                    Eigen::Triplet<double>(kz, kz - Num_var, 1)
                                });
                                R(kz) = -2 * bclower;
                            }
                            
                            // Left boundary
                            if (j == 1 && i > 1 && i < Ny1) {
                                Trip.insert(Trip.end(), {
                                    Eigen::Triplet<double>(kz, kz, 1), 
                                    Eigen::Triplet<double>(kz, kz + Num_var * Ny1, -1)
                                });
                            }
                            
                            // Right boundary
                            if (j == Nx && i > 1 && i < Ny1) {
                                Trip.insert(Trip.end(), {
                                    Eigen::Triplet<double>(kz, kz, 1), 
                                    Eigen::Triplet<double>(kz, kz - Num_var * Ny1, -1)
                                });
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
                            Trip.insert(Trip.end(), {
                                // vzs3, vzs3, vzs3, vzs2, vzs4
                                Eigen::Triplet<double>(kz, kz, -2 * (ETAXY / dx2 + ETAXY / dy2) - RHOX(i, j) / dt), 
                                Eigen::Triplet<double>(kz, kz - Ny1 * Num_var, ETAXY / dx2), 
                                Eigen::Triplet<double>(kz, kz + Ny1 * Num_var, ETAXY / dx2),
                                Eigen::Triplet<double>(kz, kz - Num_var, ETAXY / dy2), 
                                Eigen::Triplet<double>(kz, kz + Num_var, ETAXY / dy2)
                            });
                            // Right part
                            R(kz) = -RHOX(i, j) * (VZ0(i, j) / dt) - (SZX2 - SZX1) / dx - (SZY2 - SZY1) / dy;
                        }
                    }
                }
            }

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
                    const int kp    = (j * Ny1 + i) * Num_var;
                    // Reload solution
                    pt(i, j)        = S(kp) * ptscale;
                    vxs(i, j)       = S(kp + 1);
                    vys(i, j)       = S(kp + 2);
                    pf(i, j)        = S(kp + 3) * pfscale;
                    vxD(i, j)       = S(kp + 4);
                    vyD(i, j)       = S(kp + 5);
                    if (antiplane) {
                        vzs(i, j)   = S(kp + 6);
                    }
                }
            }

            Vmax = VSLIPB.maxCoeff();
            
            /*
            if (dt > 1e4 && Vmax < 1e-7) {
                double avgpt = pt.sum() / (double)(pt.rows() * pt.cols()); //calculate average total pressure
                double diffpt = (PCONF + PTFDIFF) - avgpt;
                pt += Eigen::MatrixXd::Constant(Ny1, Nx1, diffpt);
            }
            */
            
            // Velocity change
            if (antiplane) {
                vzs.col(Nx) = vzs.col(Nx - 1);
                DVZ0        = vzs - VZ0;
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
                    EXY(i, j)           = .5 * ((vxs(i + 1, j) - vxs(i, j)) / dy + (vys(i, j + 1) - vys(i, j)) / dx);
                    const double KXY    = divplus(dt * GGG(i, j), ETA(i, j));
                    SXY(i, j)           = 2 * ETA(i, j) * EXY(i, j) * KXY + SXY0(i, j) * (1 - KXY);

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
                    EXX(i, j)           = (2 * (vxs(i, j) - vxs(i, j - 1)) / dx - (vys(i, j) - vys(i - 1, j)) / dy) / 3.;
                    EYY(i, j)           = (2 * (vys(i, j) - vys(i - 1, j)) / dy - (vxs(i, j) - vxs(i, j - 1)) / dx) / 3.;
                    const double KXX    = divplus(dt * GGGP(i, j), ETAP(i, j));
                    SXX(i, j)           = 2 * ETAP(i, j) * EXX(i, j) * KXX + SXX0(i, j) * (1 - KXX);
                    SYY(i, j)           = 2 * ETAP(i, j) * EYY(i, j) * KXX + SYY0(i, j) * (1 - KXX);
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
                        EII(i, j) = sqrt(.5 * (pow(EXX(i, j), 2) + pow(EYY(i, j), 2)) + 
                                        (square_block(EXY.block(i - 1, j - 1, 2, 2)) + 
                                        square_block(EZX.block(i - 1, j - 1, 2, 2)) + 
                                        square_block(EZY.block(i - 1, j - 1, 2, 2))) / 4.
                                    );
                        
                        // Second strain rate invariant SII
                        // SXY term is averaged from four surrounding basic nodes
                        SII(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + pow(SYY(i, j), 2)) + 
                                        (square_block(SXY.block(i - 1, j - 1, 2, 2)) + 
                                        square_block(SZX.block(i - 1, j - 1, 2, 2)) + 
                                        square_block(SZY.block(i - 1, j - 1, 2, 2))) / 4.
                                    );
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

            int ynpl    = 0;
            double ddd  = 0;
            dtlapusta   = 1e7;
            OM5         = OM;

            // Power law plasticity model Yi et al., 2018
            // Journal of Offshore Mechanics and Arctic Engineering
            double dtslip = 1e30;
            // power_law(timestep, ynpl, ddd);

            if (timestep > tyield) {
            for (int i = 0; i < Ny; i++) {
                // if (y(i)  >= upper_block && y(i)  <= lower_block) {
                if (i == line_fault) {
                    for (int j = 0; j < Nx; j++) {
                        // reducing matrix calls
                        const double arsf_temp = ARSF(i, j);
                        const double brsf_temp = BRSF(i, j);
                        const double lrsf_temp = LRSF(i, j);
                        const double fric_temp = FRIC(i, j);
                        const double eta0_temp = ETA0(i, j);

                        const double sxx_temp = SXX.block(i, j, 2, 2).sum();
                        const double syy_temp = SYY.block(i, j, 2, 2).sum();

                        // SXX, pt are averaged from four surrounding pressure nodes
                        if (antiplane) {
                            SIIB(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + 
                                                    pow(SYY(i, j), 2)) + 
                                                    pow(SXY(i, j), 2) + 
                                                    pow(SZX(i, j), 2) + 
                                                    pow(SZY(i, j), 2)
                                                );                        
                        } else {
                            SIIB(i, j)    = sqrt(pow(SXY(i, j), 2) + 
                                                (pow(sxx_temp, 2) + 
                                                pow(syy_temp, 2) + 
                                                pow(sxx_temp + syy_temp, 2))/ 32.
                                            );
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
                        double ETAPL        = SIIB1 / 2. / EIISLIP;
                        double ETAVP        = 1. / (1. / eta0_temp + 1. / ETAPL);
                        // Compute new stress invariant
                        double SIIB2        = siiel * ETAVP / (GGG(i, j) * dt + ETAVP);
                        const double DSIIB1 = SIIB2 - SIIB1;
                        
                        // Compute slip velocity for current stress invariant and state
                        double V            = 2 * V0 * sinh(max(SIIB2, 0.) / arsf_temp / prB) * 
                                            exp(-(brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                        
                        EIISLIP = V / dx / 2.;
                        
                        // Compute new ETAVP
                        ETAPL               = SIIB2 / 2. / EIISLIP;
                        ETAVP               = 1. / (1. / eta0_temp + 1. / ETAPL);
                        // Compute new stress invariant
                        const double DSIIB2 = siiel * ETAVP / (GGG(i, j) * dt + ETAVP) - SIIB2;
                        double SIIB4        = 0.;
                        
                        if ((DSIIB1 >= 0 && DSIIB2 <= 0) || (DSIIB1 <= 0 && DSIIB2 >= 0)) {
                            double DSIIB = 1e9;
                            
                            while(abs(DSIIB) > 1e-3) {
                                SIIB4   = (SIIB1 + SIIB2) / 2.;
                                
                                //Compute slip velocity for current stress invariant and state
                                V       = 2 * V0 * sinh(max((SIIB4), 0.) / arsf_temp / prB) * 
                                        exp( - (brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                                
                                EIISLIP = V / dx / 2.;
                                
                                // Compute new ETAVP
                                ETAPL   = SIIB4 / 2. / EIISLIP;
                                ETAVP   = 1. / (1. / eta0_temp + 1. / ETAPL);
                                // Compute new stress invariant
                                DSIIB   = siiel * ETAVP / (GGG(i, j) * dt + ETAVP) - SIIB4;
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
                        SigmaY(i, j)        = syield;
                        VSLIPB(i, j)        = V;
                        SII_fault(i, j)     = SIIB4;

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
                            ynn         = true;
                            DSY(i, j)   = SIIB(i, j) - syield;
                            ddd         += pow(DSY(i, j), 2);
                            ynpl++;
                        }

                        // Update viscosity
                        const double A = syield / siiel;
                        if (A < 1) {
                            // New viscosity for the basic node
                            etapl = dt * GGG(i, j) * A / (1 - A);
                            if (etapl < eta0_temp) {
                                // Update plastic nodes
                                ETA5(i, j)      = pow(etapl, 1 - etawt) * pow(ETA(i, j), etawt);
                                YNY(i, j)       = 1;
                                // Count yelding nodes
                                if (!ynn) {
                                    DSY(i, j)   = SIIB(i, j) - syield;
                                    ddd         += pow(DSY(i, j), 2);
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
                dtpl    = dt / dtkoef;
                yn      = true;
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
                dtx     = dx / dtkoefv * stpmaxcur / maxvxs;
                yn      = true;
            }
            dty = dt;
            if (dt > dy * stpmaxcur / maxvys) {
                dty     = dy / dtkoefv * stpmaxcur / maxvys;
                yn      = true;
            }
            if (antiplane) {
                dtz = dt;
                if(dt > dy * stpmaxcur / maxvzs) {
                    dtz     = dx / dtkoefv * stpmaxcur / maxvzs;
                    yn      = true;
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
                dtx     = dx / dtkoefv * stpmaxcur / maxvxs;
                yn      = true;
            }
            if (dt > dy * stpmaxcur / maxvys) {
                dty     = dy / dtkoefv * stpmaxcur / maxvys;
                yn      = true;
            }
            if (antiplane && dt > dy * stpmaxcur / maxvzs) {
                dtz     = dy / dtkoefv * stpmaxcur / maxvzs;
                yn      = true;
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
                yn      = true;
                dtslip  = dtslip / dtkoefv;
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
        process_p_cells();

        // /////////////////////////////////////////////////////////////////////////////////////// 
        // Move markers by nodal velocity field
        // /////////////////////////////////////////////////////////////////////////////////////// 
        move_markers();

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
        
        write_output(timestep);
    }
}
