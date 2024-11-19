#include "move_markers.h"

using namespace std;
using namespace H5;

//This part of the program is currently unused. This is done by setting the for loop condition (int m=0; m<0; m++)
/*=======================Ginkgo===============================
Todo:
 - Ask Luca if he intends to reenable this part of the code
*/


void move_markers(){

    // #pragma omp parallel for // 3x slower for n = 4
    for (int m = 0; m < 0; m++) {
    // for (int m = 0; m < marknum; m++) {
        Eigen::Vector4d vxm = Eigen::Vector4d::Zero(), vym = Eigen::Vector4d::Zero(), spm = Eigen::Vector4d::Zero(); // Runge-Kutta velocity, spin array


        auto vxm = gko::matrix::Dense<double>::create(exec, gko::dim<2>(4,1)); //Runge-Kutta velocity
        auto vym = gko::matrix::Dense<double>::create(exec, gko::dim<2>(4,1)); //Runge-Kutta velocity
        auto spm = gko::matrix::Dense<double>::create(exec, gko::dim<2>(4,1)); // spin array

        // Save marker position
        const double xold = xm[m];
        const double yold = ym[m];
        for (int rk = 0; rk < 4; rk++) {
            // vx - velocity interpolation
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances                
            int j = check_bounds(fix_towards_zero(xm[m] / dx), Nx);
            int i = check_bounds(fix_towards_zero((ym[m] + dy / 2.) / dy), Ny);

            //Distances
            double dxm = (xm[m] - xvx(j)) / dx;
            double dym = (ym[m] - yvx(i)) / dy;
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
            j = check_bounds(fix_towards_zero((xm[m] + dx / 2.) / dx), Nx);
            i = check_bounds(fix_towards_zero(ym[m] / dy), Ny);

            //Distances
            dxm = (xm[m] - xvy(j)) / dx;
            dym = (ym[m] - yvy(i)) / dy;
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
            j = check_bounds(fix_towards_zero((xm[m]) / dx), Nx);
            i = check_bounds(fix_towards_zero((ym[m]) / dy), Ny);
            
            //Distances
            dxm = (xm[m] - x(j)) / dx;
            dym = (ym[m] - y(i)) / dy;
            // Weights
            // Interpolation ESP = .5 *(dVy / dx - dVx / dy) for the marker
            spm(rk) = ESP(i, j) * (1 - dxm) * (1 - dym) + ESP(i + 1, j) * (1 - dxm) * dym + ESP(i, j + 1) * dxm * (1 - dym) + ESP(i + 1, j + 1) * dxm * dym;
            
            // Moving between A, B, C, D points
            if (rk < 2) {
                // Moving A -> B and A -> C
                xm[m] = xold + vxm(rk) * dt / 2.;
                ym[m] = yold + vym(rk) * dt / 2.;
            } else if (rk == 3) {
                // Moving A -> D
                xm[m] = xold + vxm(rk) * dt;
                ym[m] = yold + vym(rk) * dt;
            }
        }         
        // Rotate stress on marker according to its spin
        // Compute amount of rotation from spin rate:
        // Espin = .5 *(dvy / dx - dvx / dy) i.e. positive for clockwise rotation
        // (when x axis is directed rightward and y axis is directed downward)
        const double dspeff = spm(0) * dt;
        // Save old stresses
        const double msxxold = sxxm[m];
        const double msyyold = syym[m];
        const double msxyold = sxym[m];
        sxym[m] = .5 * (msxxold - msyyold) * sin(2 * dspeff) + msxyold * cos(2 * dspeff);
        sxxm[m] = msxxold * pow(cos(dspeff), 2) + msyyold * pow(sin(dspeff), 2) - msxyold * sin(2 * dspeff);
        syym[m] = msxxold * pow(sin(dspeff), 2) + msyyold * pow(cos(dspeff), 2) + msxyold * sin(2 * dspeff);
        
        // Move markers
        xm[m] = xold + (vxm(0) + 2 * vxm(1) + 2 * vxm(2) + vxm(3)) / 6. * dt;
        ym[m] = yold + (vym(0) + 2 * vym(1) + 2 * vym(2) + vym(3)) / 6. * dt;
        
        // Recycling
        if (xm[m] < 0) {
            xm[m] = xm[m] + xsize;
        }
        if (xm[m] > xsize) {
            xm[m] = xm[m] - xsize;
        }
    }
}
