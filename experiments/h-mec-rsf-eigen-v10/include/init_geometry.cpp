#include "init_geometry.h"



// void init_geometry(VecXd &xm, VecXd &ym, VecXd &porm, VecXd &etam, VecXd &etasm, VecXd &t_marker, VecXd &rhom, VecXd &kkkm, VecXd &cohescm, VecXd &cohestm, VecXd &frictcm, VecXd &dilatcm, VecXd &fricttm, int marknum){
void init_geometry(){
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
            etam(m)     = 1e23;
            etasm(m)    = 1e23;
            rhom(m)     = 2800;
            kkkm(m)     = 2e-16; // * (dy_faultw)^2;
            cohescm(m)  = cohes * 1e3;
            cohestm(m)  = cohes * 1e3;
            frictcm(m)  = .8;
            dilatcm(m)  = dilatation;
            fricttm(m)  = tensile;
        }
    }
}
