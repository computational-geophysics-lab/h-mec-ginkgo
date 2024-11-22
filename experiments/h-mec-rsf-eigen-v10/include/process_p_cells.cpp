#include "process_p_cells.h"


using namespace std;
using namespace H5;

void process_p_cells(){
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
                Eigen::Matrix2d temp;
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
}