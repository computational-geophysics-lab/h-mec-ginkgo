#include "process_p_cells.h"


using namespace std;
using namespace H5;

/* ======================== Ginkgo =====================
Todo:
Convert everything to Ginkgo. Example Sxx.block(i,j,2,2).sum(); needs to be implemented with ginkgo
*/

void process_p_cells(){
    // #pragma omp parallel for collapse(2) // about 2x faster with n = 4 but breaks the simulation
    for (int i = 1; i < Ny1-1; i++) {
        for (int j = 1; j < Nx1-1; j++) {
            if (antiplane) {
                // EXY term is averaged from four surrounding basic nodes
                EII[i * Nx1 + j] = sqrt(pow(EXX[i * Nx1 + j], 2) + (square_block(EXY, i - 1, j - 1, Nx)) + square_block(EZX, i - 1, j - 1, Nx) + square_block(EZY, i - 1, j - 1, Nx)) / 4.;
                
                // Second strain rate invariant SII
                // SXY term is averaged from four surrounding basic nodes
                SII[i * Nx + j] = sqrt(.5 * (pow(SXX[i * Nx + j], 2) + pow(SYY[i * Nx + j], 2)) + (square_block(SXY, i - 1, j - 1, Nx)) + square_block(SZX, i - 1, j - 1, Nx)) + square_block(SZY, i - 1, j - 1, Nx) / 4.;
                
                // Dissipation
                const double DISXY = pow(SXY[i * Nx + j], 2) / ETA[i * Nx + j] + pow(SXY[(i - 1) * Nx+ j], 2) / ETA[(i - 1) * Nx+ j] + pow(SXY[i * Nx + j - 1], 2) / ETA[i * Nx + j - 1] + pow(SXY[(i - 1) * Nx + j - 1], 2) / ETA[(i - 1) * Nx + j - 1];
                const double DISZX = pow(SZX[i * Nx + j], 2) / ETA[i * Nx + j] + pow(SZX[(i - 1) * Nx+ j], 2) / ETA[(i - 1) * Nx+ j] + pow(SZX[i * Nx + j - 1], 2) / ETA[i * Nx + j - 1] + pow(SZX[(i - 1) * Nx + j - 1], 2) / ETA[(i - 1) * Nx + j - 1];
                const double DISZY = pow(SZY[i * Nx + j], 2) / ETA[i * Nx + j] + pow(SZY[(i - 1) * Nx+ j], 2) / ETA[(i - 1) * Nx+ j] + pow(SZY[i * Nx + j - 1], 2) / ETA[i * Nx + j - 1] + pow(SZY[(i - 1) * Nx + j - 1], 2) / ETA[(i - 1) * Nx + j - 1];
                DIS[i * Nx + j] = (pow(SXX[i * Nx + j], 2) / ETAP[i * Nx + j] + pow(SYY[i * Nx + j], 2) / ETAP[i * Nx + j] + (DISXY + DISZX + DISZY) / 2.) / 2.;
            } else {
                // EXX, SXX
                /*EXX[i * Nx1 + j] = (2 * (vxs[i * Nx1 + j] - vxs[i * Nx1 + j - 1]) / dx - (vys[i * Nx1 + j] - vys[(i - 1) * Nx1 + j]) / dy) / 3.;
                EYY[i * Nx1 + j] = (2 * (vys[i * Nx1 + j] - vys[(i - 1) * Nx1 + j]) / dy - (vxs[i * Nx1 + j] - vxs[i * Nx1 + j - 1]) / dx) / 3.;*/
                double dvxs = safe_difference(vxs[i * Nx1 + j], vxs[i * Nx1 + j - 1], significant_digits);
                double dvys = safe_difference(vys[i * Nx1 + j], vys[(i-1) * Nx1 + j], significant_digits);
                EXX[i * Nx1 + j] = (2 * dvxs / dx - dvys / dy) / 3.;
                EYY[i * Nx1 + j] = (2 * dvys / dy - dvxs / dx) / 3.;
                const double KXX = divplus(dt * GGGP[i * Nx1 + j], ETAP[i * Nx1 + j]);
                SXX[i * Nx1 + j] = 2 * ETAP[i * Nx1 + j] * EXX[i * Nx1 + j] * KXX + SXX0[i * Nx1 + j] * (1 - KXX);
                SYY[i * Nx1 + j] = 2 * ETAP[i * Nx1 + j] * EYY[i * Nx1 + j] * KXX + SYY0[i * Nx1 + j] * (1 - KXX);

                // Compute stress and strain rate invariants and dissipation
                auto temp_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(2,2));
                double* temp = temp_gko->get_values();
                temp[0] = SXY[(i - 1) * Nx + j - 1] / ETA[(i - 1) * Nx + j - 1];
                temp[1] = SXY[(i - 1) * Nx + j] / ETA[(i - 1) * Nx+ j];
                temp[2] = SXY[i * Nx + j - 1] / ETA[i * Nx + j - 1];
                temp[3] = SXY[i * Nx + j] / ETA[i * Nx + j];
                //temp << SXY[(i - 1) * Nx + j - 1] / ETA[(i - 1) * Nx + j - 1], SXY[(i - 1) * Nx+ j] / ETA[(i - 1) * Nx+ j], SXY[i * Nx + j - 1] / ETA[i * Nx + j - 1], SXY[i * Nx + j] / ETA[i * Nx + j];

                // EXY term is averaged from four surrounding basic nodes
                EII[i * Nx1 + j] = sqrt(pow(EXX[i * Nx1 + j], 2) + square_block(EXY, i - 1, j - 1, Nx) / 4.);
                EIIVP[i * Nx1 + j] = sqrt(.5 * (pow(SXX[i * Nx1 + j] / (2 * ETAP[i * Nx1 + j]), 2) + pow(SYY[i * Nx1 + j] / (2 * ETAP[i * Nx1 + j]), 2)) + square_block(temp, 0, 0, 2) / 16.);
                // Second strain rate invariant SII
                // SXY term is averaged from four surrounding basic nodes
                SII[i * Nx1 + j] = sqrt(.5 * (pow(SXX[i * Nx1 + j], 2) + pow(SYY[i * Nx1 + j], 2)) + square_block(SXY, i - 1, j - 1, Nx) / 4.);

                // Dissipation
                double DISXY = (pow(SXY[i * Nx + j], 2) / ETA[i * Nx + j] + pow(SXY[(i - 1) * Nx+ j], 2) / ETA[(i - 1) * Nx+ j] + pow(SXY[i * Nx + j - 1], 2) / ETA[i * Nx + j - 1] + pow(SXY[(i - 1) * Nx + j - 1], 2) / ETA[(i - 1) * Nx + j - 1]) / 2.;
                DIS[i * Nx1 + j] = pow(SXX[i * Nx1 + j], 2) / (2 * ETAP[i * Nx1 + j]) + pow(SYY[i * Nx1 + j], 2) / (2 * ETAP[i * Nx1 + j]) + DISXY;
            }
            
            double PT0_ave, PF0_ave;
            if (i < Ny - 1) {
                pt_ave[i * Nx1 + j] = (pt[i * Nx1 + j] + pt[(i + 1) * Nx1 + j]) / 2.;
                pf_ave[i * Nx1 + j] = (pf[i * Nx1 + j] + pf[(i + 1) * Nx1 + j]) / 2.;
                PT0_ave = (PT0[i * Nx1 + j] + PT0[(i + 1) * Nx1 + j]) / 2.;
                PF0_ave = (PF0[i * Nx1 + j] + PF0[(i + 1) * Nx1 + j]) / 2.;
            } else {
                pt_ave[i * Nx1 + j] = pt[i * Nx1 + j];
                pf_ave[i * Nx1 + j] = pf[i * Nx1 + j];
                PT0_ave = PT0[i * Nx1 + j];
                PF0_ave = PF0[i * Nx1 + j];
            }
            
            // Compute elastic and viscous compaction
            VIS_COMP[i * Nx1 + j] = (pt_ave[i * Nx1 + j] - pf_ave[i * Nx1 + j]) / (ETAB[i * Nx1 + j] * (1 - POR[i * Nx1 + j]));
            // Drained compressibility
            const double BETADRAINED = (1 / GGGB[i * Nx1 + j] + BETASOLID) / (1 - POR[i * Nx1 + j]);
            // Biott - Willis koefficient
            const double KBW = 1 - BETASOLID / BETADRAINED;
            EL_DECOM[i * Nx1 + j] = BETADRAINED * (pt_ave[i * Nx1 + j] - PT0_ave - KBW * pf_ave[i * Nx1 + j] + KBW * PF0_ave) / dt;
        }
    }
}