#include "LR_setup.h"

using namespace std;
using namespace H5;


/* ====================== Ginkgo ===============================
Todo:
 - No changes made yet, switch to ginkgo
*/


void LR_setup(Eigen::SparseMatrix<double> &L, double pfscale, double ptscale){

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
                    R[kx] = 2 * bcupper;
                }
                
                // Lower boundary
                // prescribed velocity
                if (i == Ny && j < Nx) {
                    Trip.insert(Trip.end(), {
                        Eigen::Triplet<double>(kx, kx, 1), 
                        Eigen::Triplet<double>(kx, kx - Num_var, 1)
                    });
                    R[kx] = 2 * bclower;
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
                double ETAXY2 = ETA[i * Nx + j];
                double ETAXX1 = ETAP[i * Nx1 + j];
                double ETAXX2 = ETAP(i, j + 1);
                // Shear modulus
                const double GXY1 = GGG(i - 1, j);
                const double GXY2 = GGG[i * Nx + j];
                const double GXX1 = GGGP[i * Nx1 + j];
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
                const double SXY2 = SXY0[i * Nx + j] * (1 - KXY2);
                const double SXX1 = SXX0[i * Nx1 + j] * (1 - KXX1);
                const double SXX2 = SXX0(i, j + 1) * (1 - KXX2);
                // Density derivatives
                const double dRHOdx = (RHOX(i, j + 1) - RHOX(i, j - 1)) / (2. * dx);
                const double dRHOdy = (RHO[i * Nx + j] - RHO(i - 1, j)) / dy;
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
                R[kx] = -RHOX(i, j) * (inertia * VX0(i, j) / dt + gx) - (SXX2 - SXX1) / dx - (SXY2 - SXY1) / dy;
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
                R[ky] = -RHOY(i, j) * (inertia * VY0(i, j) / dt + gy) - (SYY2 - SYY1) / dy - (SXY2 - SXY1) / dx;
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
                        R[kz] = 2 * bclower;
                    }
                    
                    // Lower boundary
                    if (i == Ny1 && j < Nx1) {
                        Trip.insert(Trip.end(), {
                            Eigen::Triplet<double>(kz, kz, 1), 
                            Eigen::Triplet<double>(kz, kz - Num_var, 1)
                        });
                        R[kz] = -2 * bclower;
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
                        R[kz]=0;
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
                    R[kz] = -RHOX(i, j) * (VZ0(i, j) / dt) - (SZX2 - SZX1) / dx - (SZY2 - SZY1) / dy;
                }
            }
        }
    }

    L.setFromTriplets(Trip.begin(), Trip.end()); // Build Sparse Matrix
}