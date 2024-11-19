#include "init_geometry.h"


void init_geometry(){
    bool randomize_markers = true;
    std::cout << "\nMarker positions are randomized = " << randomize_markers << ". -> init_geometry.cpp, line 5.\n";
    for (int m = 0; m < marknum; m++) {
        // Define randomized regular coordinates.
        xm[m] = xbeg + floor(m / Ny_markers) * dxms + randd() * dxms * randomize_markers;
        ym[m] = ybeg + m % Ny_markers * dyms + randd() * dyms * randomize_markers;
        // Matrix
        porm[m] = .01 * (1 + .0 * (randomize_markers * randd() - .5));
        etam[m] = etasm[m] * exp(-alpha * porm[m]);

        // Air, wedge, slab
        if (ym[m] < upper_block || ym[m] > lower_block) {
            t_marker[m] = -1;
            etam[m]     = 1e23;
            etasm[m]    = 1e23;
            rhom[m]     = 2800;
            kkkm[m]     = 2e-16; // * (dy_faultw)^2;
            cohescm[m]  = cohes * 1e3;
            cohestm[m]  = cohes * 1e3;
            frictcm[m]  = .8;
            dilatcm[m]  = dilatation;
            fricttm[m]  = tensile;
        }
    }
}
