#include "run_simulation.h"
using namespace std;
using namespace H5;


/*
Todo:  - Check all ->at statements. Are i, j in the right or wrong order?
 - Use different method to assemble the stencil. This looks promising: https://ginkgo-project.github.io/ginkgo-generated-documentation/doc/develop/classgko_1_1matrix__assembly__data.html#a0927a396f936167366bce13745156cf5
 The add_value function should do exactly what I need, and what emplace_back couldn't do
*/

void run_simulation(int &timestep){
    bool save_first_matrix = true;
    std::cout << "Check if both the input vector is identical to eigen version and also if the L matrix is identical. \nMaybe I have to go line by line and compare to make sure I calculate it as in the eigen version.\nUse the clean matrix assembly method, it will also make the code more readable hopefully!" << std::endl;
    // Declaration of matrices that are set to zero each timestep
    std::unique_ptr<gko::matrix::Dense<>> ETA0SUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* ETA0SUM = ETA0SUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> COHCSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* COHCSUM = COHCSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> FRICSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* FRICSUM = FRICSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> DILCSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* DILCSUM = DILCSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> COHTSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* COHTSUM = COHTSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> FRITSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* FRITSUM = FRITSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> WTSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* WTSUM = WTSUM_gko->get_values(); // Ny x Nx

    // interpolate ETA and RHO to nodal points
    // basic nodes
    std::unique_ptr<gko::matrix::Dense<>> RHOSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* RHOSUM = RHOSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> ETASUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* ETASUM = ETASUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> KKKSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* KKKSUM = KKKSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> TTTSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* TTTSUM = TTTSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> SXYSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* SXYSUM = SXYSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> GGGSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* GGGSUM = GGGSUM_gko->get_values(); // Ny x Nx

    // LDZ
    std::unique_ptr<gko::matrix::Dense<>> OM0SUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* OM0SUM = OM0SUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> OMSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* OMSUM = OMSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> ARSFSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* ARSFSUM = ARSFSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> BRSFSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* BRSFSUM = BRSFSUM_gko->get_values(); // Ny x Nx
    std::unique_ptr<gko::matrix::Dense<>> LRSFSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny, Nx)); // Ny x Nx
    double* LRSFSUM = LRSFSUM_gko->get_values(); // Ny x Nx
 
    // Pressure nodes
    std::unique_ptr<gko::matrix::Dense<>> ETAPSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* ETAPSUM = ETAPSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> ETAP0SUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* ETAP0SUM = ETAP0SUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> ETAB0SUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* ETAB0SUM = ETAB0SUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> PORSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* PORSUM = PORSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> SXXSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* SXXSUM = SXXSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> SYYSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* SYYSUM = SYYSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> GGGPSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* GGGPSUM = GGGPSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> WTPSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* WTPSUM = WTPSUM_gko->get_values(); // Ny1 x Nx1

    // Vx nodes
    std::unique_ptr<gko::matrix::Dense<>> RHOXSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* RHOXSUM = RHOXSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> RHOFXSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* RHOFXSUM = RHOFXSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> ETADXSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* ETADXSUM = ETADXSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> PORXSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* PORXSUM = PORXSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> WTXSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* WTXSUM = WTXSUM_gko->get_values(); // Ny1 x Nx1

    // Vy nodes
    std::unique_ptr<gko::matrix::Dense<>> RHOYSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* RHOYSUM = RHOYSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> RHOFYSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* RHOFYSUM = RHOFYSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> ETADYSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* ETADYSUM = ETADYSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> PORYSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* PORYSUM = PORYSUM_gko->get_values(); // Ny1 x Nx1
    std::unique_ptr<gko::matrix::Dense<>> WTYSUM_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* WTYSUM = WTYSUM_gko->get_values(); // Ny1 x Nx1
 
    std::unique_ptr<gko::matrix::Dense<>> ETA50_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1)); // Ny1 x Nx1
    double* ETA50 = ETA50_gko->get_values(); // Ny1 x Nx1

    auto L = gko::share(gko::matrix::Csr<double>::create(exec)); // Matrix of coefficients in the left part
    gko::matrix_data<> empty_csr_matrix{gko::dim<2>(1,1)}; // Contains nothing, used to reset the L matrix after each timestep
    std::unique_ptr<gko::matrix::Dense<>> L_dense_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(N,N));
    double* L_dense = L_dense_gko->get_values();

    // To subtract gko matrices I only found the sub_scaled function, that also takes a factor alpha, where alpha has to be a gko::matrix, I just set it to 1.0 I assume there is a better solution, but I didn't find it
    std::unique_ptr<gko::matrix::Dense<>> sub_alpha = gko::matrix::Dense<double>::create(exec, gko::dim<2>(1,1)); // Need a alpha parameter when subtracting matrices -> set to 1.0
    sub_alpha->at(0,0)=1.0; // Need a alpha parameter when subtracting matrices -> set to 1.0

    // nucleation size
    const double hstar = pi / 2. * shearmod * lrsfm[1] * brsfm[1] / pow(brsfm[1] - arsfm[1], 2) / PTFDIFF;
    // cohesive zone size
    const double coh = 9. / 32. * pi * shearmod * lrsfm[1] / brsfm[1] / PTFDIFF;
    // Print information about discretization
    cout << ">> VW width = "                    << (((TS_4 + TS_3) / 2) - ((TS_2 + TS_1) / 2)) / 1e3    << " (km)" << endl;
    cout << ">> Critical nucleation size = "    << hstar / 1e3                                          << " (km)" << endl;
    cout << ">> Cohesive zone = "               << coh                                                  << " [m]"  << endl;



    // /////////////////////////////////////////////////////////////////////////////////////// 
    // actual computations start here
    // /////////////////////////////////////////////////////////////////////////////////////// 


    for (; timestep <= num_timesteps; timestep++) {
        
        RHOSUM_gko->fill(0);
        ETASUM_gko->fill(0);
        KKKSUM_gko->fill(0);
        TTTSUM_gko->fill(0);
        SXYSUM_gko->fill(0);
        GGGSUM_gko->fill(0);
        ETA_gko->fill(0);
        ETA0SUM_gko->fill(0);
        COHCSUM_gko->fill(0);
        FRICSUM_gko->fill(0);
        DILCSUM_gko->fill(0);
        COHTSUM_gko->fill(0);
        FRITSUM_gko->fill(0);
        WTSUM_gko->fill(0);
        OM0SUM_gko->fill(0);
        OMSUM_gko->fill(0);
        ARSFSUM_gko->fill(0);
        BRSFSUM_gko->fill(0);
        LRSFSUM_gko->fill(0);
        ETAPSUM_gko->fill(0);
        ETAP0SUM_gko->fill(0);
        ETAB0SUM_gko->fill(0);
        PORSUM_gko->fill(0);
        SXXSUM_gko->fill(0);
        SYYSUM_gko->fill(0);
        GGGPSUM_gko->fill(0);
        WTPSUM_gko->fill(0);
        RHOXSUM_gko->fill(0);
        RHOFXSUM_gko->fill(0);
        ETADXSUM_gko->fill(0);
        PORXSUM_gko->fill(0);
        WTXSUM_gko->fill(0);
        RHOYSUM_gko->fill(0);
        RHOFYSUM_gko->fill(0);
        ETADYSUM_gko->fill(0);
        PORYSUM_gko->fill(0);
        WTYSUM_gko->fill(0);


        //Debugging:
        int i_max = 0;
        int j_max = 0;

        // Cycle on markers
        #pragma omp parallel for // about 3-4x faster with n = 4
        for (int m = 0; m < marknum; m++) {
            double cohescmm, cohestmm, frictcmm, dilatcmm, fricttmm, etasmm0, etamm0, etamm, rhomm, etadm;

            // Marker properties
            double kkkmm = kkkm[m] * pow(porm[m] / POR0, 3);
            // Checking permeability limits
            kkkmm = enforce_bounds(kkkmm, kkkmin, kkkmax);


            // Viscosity of porous matrix
            if (t_marker[m] != 0) {
                cohescmm    = cohescm[m] * (1 - porm[m]); // * exp(-alpha * porm[m]);
                cohestmm    = cohestm[m] * (1 - porm[m]); // * exp(-alpha * porm[m]);
                frictcmm    = frictcm[m];
                dilatcmm    = dilatcm[m];
                fricttmm    = fricttm[m];
                etasmm0     = etasm[m];
                etamm0      = etasm[m] * exp(-alpha * porm[m]);
                etamm       = etam[m];
                // total density
                rhomm       = rhom[m] * (1 - porm[m]) + rhofm[m] * porm[m];
            } else {
                cohescmm    = cohescm[m];
                cohestmm    = cohestm[m];
                frictcmm    = frictcm[m];
                dilatcmm    = 0;
                fricttmm    = fricttm[m];
                etasmm0     = etasm[m];
                etamm0      = etasm[m];
                etamm       = etam[m];
                // total density
                rhomm       = rhom[m];
            }
            
            etamm0  = enforce_bounds(etamm0, etamin, etamax); // Matrix viscosity
            etamm   = enforce_bounds(etamm, etamin, etamax); // Effective viscosity
            etadm   = etafm[m] / kkkmm; // Darsi "viscosity"

            // Interpolate to basic nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            int j = check_bounds(fix_towards_zero(xm[m] / dx), Nx);
            int i = check_bounds(fix_towards_zero(ym[m] / dy), Ny);

            double dxm = (xm[m] - x[j]) / dx;
            double dym = (ym[m] - y[i]) / dy;
            
            // merged similar computations of matlab version
            // The computations on the 4 elements are now done on a 2x2-sub-block of the matrix
            // Similarly wtm is a 2x2 matrix that stores the corresponding value
            // Sub-block += variable * wtm

            double wtm[4] = {(1 - dxm) * (1 - dym), dxm * (1 - dym), (1 - dxm) * dym, dxm * dym};


            add_block(ETASUM, wtm, etamm, i, j, Nx, Ny);
            add_block(RHOSUM, wtm, rhomm, i, j, Nx, Ny);
            add_block(KKKSUM, wtm, kkkmm, i, j, Nx, Ny);
            add_block(TTTSUM, wtm, t_marker[m], i, j, Nx, Ny);
            add_block(SXYSUM, wtm, sxym[m], i, j, Nx, Ny);
            add_block(GGGSUM, wtm, 1.0/gm[m], i, j, Nx, Ny);
            add_block(ETA0SUM, wtm, etamm0, i, j, Nx, Ny);
            add_block(COHTSUM, wtm, cohestmm, i, j, Nx, Ny);
            add_block(FRITSUM, wtm, fricttmm, i, j, Nx, Ny);
            add_block(COHCSUM, wtm, cohescmm, i, j, Nx, Ny);
            add_block(FRICSUM, wtm, frictcmm, i, j, Nx, Ny);
            add_block(DILCSUM, wtm, dilatcmm, i, j, Nx, Ny);
            add_block(WTSUM, wtm, 1.0, i, j, Nx, Ny);

            // Interpolate to pressure nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix_towards_zero((xm[m] + dx / 2.) / dx), Nx);
            i = check_bounds(fix_towards_zero((ym[m] + dy / 2.) / dy), Ny);
            if (i>i_max) {
                i_max = i;
            }
            if (j>j_max) {
                j_max = j;
            }
            
            dxm = (xm[m] - xp[j]) / dx;
            dym = (ym[m] - yp[i]) / dy;
            if (m<200 && false) {
                std::cout << "xm[" << m << "]=" << xm[m] << " xp[" << j << "]=" << xp[j] << "  ym[" << m << "]=" << ym[m] << " yp[" << i << "]=" << yp[i] << std::endl;
            }
            wtm[0] = (1 - dxm) * (1 - dym);
            wtm[1] = dxm * (1 - dym);
            wtm[2] = (1 - dxm) * dym;
            wtm[3] = dxm * dym;

            if (m<0 && false) {
                std::cout << wtm[0] << " "<< wtm[1] << " "<< wtm[2] << " " << wtm[3] << std::endl;
            }

            add_block(ETAPSUM, wtm, etamm, i, j, Nx1, Ny1);
            add_block(PORSUM, wtm, porm[m], i, j, Nx1, Ny1);
            add_block(SXXSUM, wtm, sxxm[m], i, j, Nx1, Ny1);
            add_block(SYYSUM, wtm, syym[m], i, j, Nx1, Ny1);
            add_block(GGGPSUM, wtm, 1.0/gm[m], i, j, Nx1, Ny1);
            add_block(ETAP0SUM, wtm, etamm0, i, j, Nx1, Ny1);
            add_block(ETAB0SUM, wtm, etasmm0, i, j, Nx1, Ny1);

            add_block(WTPSUM, wtm, 1.0, i, j, Nx1, Ny1);



            // Interpolate to Vx nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix_towards_zero((xm[m]) / dx), Nx);
            i = check_bounds(fix_towards_zero((ym[m] + dy / 2.) / dy), Ny);
            
            dxm = (xm[m] - xvx[j]) / dx;
            dym = (ym[m] - yvx[i]) / dy;
            
            wtm[0] = (1 - dxm) * (1 - dym);
            wtm[1] = dxm * (1 - dym);
            wtm[2] = (1 - dxm) * dym;
            wtm[3] = dxm * dym;

            add_block(RHOXSUM, wtm, rhomm, i, j, Nx1, Ny1);
            add_block(ETADXSUM, wtm, 1.0/etadm, i, j, Nx1, Ny1);
            add_block(RHOFXSUM, wtm, rhofm[m], i, j, Nx1, Ny1);
            add_block(PORXSUM, wtm, porm[m], i, j, Nx1, Ny1);
            
            add_block(WTXSUM, wtm, 1.0, i, j, Nx1, Ny1);
            
            // Interpolate to Vy nodes
            // [i, j] -------- [i, j + 1]
            //   |                |
            //   |    o m         |
            //   |                |
            // [i + 1, j] ------- [i + 1, j + 1]
            // Indexes and distances
            j = check_bounds(fix_towards_zero((xm[m] + dx / 2.) / dx), Nx);
            i = check_bounds(fix_towards_zero((ym[m]) / dy), Ny);
        
            dxm = (xm[m] - xvy[j]) / dx;
            dym = (ym[m] - yvy[i]) / dy;

            wtm[0] = (1 - dxm) * (1 - dym);
            wtm[1] = dxm * (1 - dym);
            wtm[2] = (1 - dxm) * dym;
            wtm[3] = dxm * dym;

            add_block(RHOYSUM, wtm, rhomm, i, j, Nx1, Ny1);
            add_block(ETADYSUM, wtm, 1.0/etadm, i, j, Nx1, Ny1);
            add_block(RHOFYSUM, wtm, rhofm[m], i, j, Nx1, Ny1);
            add_block(PORYSUM, wtm, porm[m], i, j, Nx1, Ny1);

            add_block(WTYSUM, wtm, 1.0, i, j, Nx1, Ny1);
        } // ends loop through markers


        // Computing ETA and RHO
        for (int i = 0; i < Ny; i++) {
            // Basic nodes
            for (int j = 0; j < Nx; j++) {
                double wtsum = WTSUM[i*Nx + j];
                if (wtsum > 0) {
                    RHO[i*Nx + j]   = RHOSUM[i*Nx + j] / wtsum;
                    KKK[i*Nx + j]   = KKKSUM[i*Nx + j] / wtsum;
                    TTT[i*Nx + j]   = TTTSUM[i*Nx + j] / wtsum;
                    GGG[i*Nx + j]   = 1. / (GGGSUM[i*Nx + j] / wtsum);
                    ETA0[i*Nx + j]  = ETA0SUM[i*Nx + j] / wtsum;
                    COHT[i*Nx + j]  = COHTSUM[i*Nx + j] / wtsum;
                    FRIT[i*Nx + j]  = FRITSUM[i*Nx + j] / wtsum;
                    COHC[i*Nx + j]  = COHCSUM[i*Nx + j] / wtsum;
                    FRIC[i*Nx + j]  = FRICSUM[i*Nx + j] / wtsum;
                    DILC[i*Nx + j]  = DILCSUM[i*Nx + j] / wtsum;
                }
            }
            // Vy nodes
            for (int j = 0; j < Nx1; j++) {
                double wtysum = WTYSUM[i*Nx1 + j];
                if (wtysum > 0) {
                    RHOY[i*Nx1 + j]  = RHOYSUM[i*Nx1 + j] / wtysum;
                    ETADY[i*Nx1 + j] = 1. / (ETADYSUM[i*Nx1 + j] / wtysum);
                    RHOFY[i*Nx1 + j] = RHOFYSUM[i*Nx1 + j] / wtysum;
                    PORY[i*Nx1 + j]  = PORYSUM[i*Nx1 + j] / wtysum;
                }
            }
        }


        for (int i = 0; i < Ny1; i++) {
            // Vx nodes
            for (int j = 0; j < Nx; j++) {
                double wtxsum = WTXSUM[i*Nx1 + j];
                if (wtxsum > 0) {
                    RHOX[i*Nx1 + j]  = RHOXSUM[i*Nx1 + j] / wtxsum;
                    ETADX[i*Nx1 + j] = 1. / (ETADXSUM[i*Nx1 + j] / wtxsum);
                    RHOFX[i*Nx1 + j] = RHOFXSUM[i*Nx1 + j] / wtxsum;
                    PORX[i*Nx1 + j]  = PORXSUM[i*Nx1 + j] / wtxsum;
                }
            }
            //Pressure nodes
            for (int j = 0; j < Nx1; j++) {
                double wtpsum = WTPSUM[i*Nx1 + j];
                if (wtpsum > 0) {
                    ETAP[i*Nx1 + j]  = ETAPSUM[i*Nx1 + j] / wtpsum;
                    POR[i*Nx1 + j]   = PORSUM[i*Nx1 + j] / wtpsum;
                    ETAB0[i*Nx1 + j] = ETAB0SUM[i*Nx1 + j] / wtpsum / POR[i*Nx1 + j];
                    if (POR[i*Nx1+j]==0) {
                        std::cout << "\nline 390: This should produce a nan: ETAB0 = " << ETAB0SUM[i*Nx1 + j] / wtpsum / POR[i*Nx1 + j];
                    }
                    ETAP0[i*Nx1 + j] = ETAP0SUM[i*Nx1 + j] / wtpsum;
                    GGGP[i*Nx1 + j]  = 1. / (GGGPSUM[i*Nx1 + j] / wtpsum);
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
        GkoSetZero(DSYLSQ_gko);


        for (iterstep = 0; iterstep < niterglobal; iterstep++) {
            // Limiting viscosity
            double etamincur = dt * shearmod * 1e-4;
            double ptscale, pfscale;

            // External P - nodes: symmetry
            copy_bounds(pt_gko, Nx1, Ny1);
            copy_bounds(pf_gko, Nx1, Ny1);

            // Basic nodes
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    if (ETA[i*Nx + j] < etamincur) {
                        ETA[i*Nx + j] = etamincur;
                    }
                    // Compute plastic strain rate
                    if (ETA[i*Nx + j] < ETA0[i*Nx + j]) {
                        double SXX_temp = BlockSum(SXX, i, j, 2, 2, Ny1, Nx1);
                        double SYY_temp = BlockSum(SYY, i, j, 2, 2, Ny1, Nx1);
                        SIIB[i*Nx + j]              = sqrt(pow(SXY[i*Nx + j], 2) + (pow(SXX_temp, 2) + pow(SYY_temp, 2) +
                                                  pow(SXX_temp  + SYY_temp, 2)) / 32.);
                        IETAPLB[i*Nx + j]           = (1. / ETA[i*Nx + j] - 1. / ETA0[i*Nx + j]);
                        EIIB[i*Nx + j]              = dy_faultw * SIIB[i*Nx + j] * IETAPLB[i*Nx + j] / 2.;
                    } else {
                        EIIB[i*Nx + j]      = 0;
                        IETAPLB[i*Nx + j]   = 0;
                    }
                }
            }


            std::cout << "run_sim 1" << std::endl;
            // Computing viscosity and dilatation in pressure nodes
            for (int i = 1; i < Ny; i++) {
                for (int j = 1; j < Nx; j++) {
                    // Compute viscoplastic viscosity
                    double IETAPL = BlockSum(IETAPLB, i - 1, j - 1, 2, 2, Ny, Nx)/4.0;
                    if (YNY0[(i - 1)*Nx + j - 1] > 0 || YNY0[i*Nx + j - 1] > 0 || YNY0[(i - 1)*Nx + j] > 0 || YNY0[i*Nx + j] > 0) {
                        ETAP[i*Nx1 + j] = 1. / (1. / ETAP0[i*Nx1 + j] + IETAPL);
                        ETAB[i*Nx1 + j] = 1. / (1. / ETAB0[i*Nx1 + j] + dy_faultw * IETAPL * POR[i*Nx1 + j]);
                        if (isnan(POR[i*Nx1 + j]==0.0) && false) {
                            std::cout << "POR[" << i*Nx1+j << "] is " << POR[i*Nx1+j];
                        }
                    } else {
                        ETAP[i*Nx1 + j] = ETAP0[i*Nx1 + j];
                        ETAB[i*Nx1 + j] = ETAB0[i*Nx1 + j];
                    }
                    // Check viscosity
                    if (ETAP[i*Nx1 + j] < etamincur) {
                        ETAP[i*Nx1 + j] = etamincur;
                    }
                    if (ETAB[i*Nx1 + j] * POR[i*Nx1 + j] < etamincur) {
                        ETAB[i*Nx1 + j] = etamincur / POR[i*Nx1 + j];
                    }
                    // Pores compressibility
                    GGGB[i*Nx1 + j] = GGGP[i*Nx1 + j] / POR[i*Nx1 + j];
                    if (isnan(GGGB[i*Nx1+j]) && true) {
                        std::cout << "POR[" << i*Nx1+j << "] is zero and therefore GGGB[" << i*Nx1+j << "] is " << GGGB[i*Nx1+j]<< "!" << std::endl;
                    }
                    if (isnan(GGGB[i * Nx1 + j]) && true) {
                        std::cout << "GGGB is nan: Cause might be that POR is 0: POR["<< i*Nx1+j << "] = " << POR[i * Nx1 + j] << std::endl;
                    }
                    // Dilation
                    // Zhao and Cai, International Journal of Rock Mechanics & Mining Sciences 47 (2010) 368–384
                    // Weak sandstone parameters
                    // double ss3 = min(max((pt[i*Nx + j] - pf[i*Nx + j]) * 1e-6, 0.), 100.); // SIGMA3, MPa
                    // double aa = aa1 + aa2 * exp(-ss3 / aa3);
                    // double bb = bb1 + bb2 * exp(-ss3 / bb3);
                    // double cc = cc1 + cc2 / 100. * pow(ss3, cc3);
                    // double dil = sin(aa * bb * (exp(-bb * gammap) - exp(-cc * gammap)) / (cc - bb) / 180. * pi);
                    DILP[i*Nx1 + j] = 0; // 2 * (dil * EIIB[(i - 1)*Nx + j - 1] + dil * EIIB[i*Nx + j - 1] + dil * EIIB[(i - 1)*Nx + j] + dil * EIIB[i*Nx + j]) / 4;
                }
            }
            std::cout << "Works for Nx-1" << std::endl;

            if (dt > 2e5) {
                double pfscale = BlockMin(ETADX, 1, 0, Ny - 1, Nx, Ny1, Nx1) * dx * 1e17 / pow(dt, 2);
            } else {
                double pfscale = BlockMin(ETADX, 1, 0, Ny - 1, Nx, Ny1, Nx1) * dx;            }

            ptscale = pfscale;

            // ///////////////////////////////////////////////////////////////////////////////////////
            // Set up L and R matrices
            // ///////////////////////////////////////////////////////////////////////////////////////


            R_gko->fill(0.0);
            L_dense_gko->fill(0.0);

            //Set the CSR matrix to zero before building it:
            L->read(empty_csr_matrix);



            std::cout << "Starting to assemble Stencil matrix" << std::endl;
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
                            L_dense_gko->at(kx,kx) += 1;
                        }

                        // Upper boundary
                        // prescribed velocity
                        if (i == 0 && j < Nx) {
                            L_dense_gko->at(kx,kx) += 1;
                            L_dense_gko->at(kx, kx+Num_var)+=1;
                            R[kx] = 2 * bcupper;
                        }

                        // Lower boundary
                        // prescribed velocity
                        if (i == Ny && j < Nx) {
                            L_dense_gko->at(kx, kx)+=1;
                            L_dense_gko->at(kx, kx-Num_var)+=1;
                            R[kx] = 2 * bclower;
                        }
                        
                        // Left boundary
                        if (j == 0 && i > 0 && i < Ny) {
                            L_dense_gko->at(kx, kx)+=1;
                            L_dense_gko->at(kx, kx+Num_var*Ny1)+=-1;
                        }
                        
                        // Right boundary
                        if (j == Nx - 1 && i > 0 && i < Ny) {
                            L_dense_gko->at(kx, kx)+=1;
                            L_dense_gko->at(kx, kx-Num_var*Ny1)+=-1;
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
                        double ETAXY1 = ETA[(i - 1)*Nx + j];
                        double ETAXY2 = ETA[i * Nx + j];
                        double ETAXX1 = ETAP[i * Nx1 + j];
                        double ETAXX2 = ETAP[i * Nx1 + j + 1];
                        // Shear modulus
                        const double GXY1 = GGG[(i - 1)*Nx + j];
                        const double GXY2 = GGG[i * Nx + j];
                        const double GXX1 = GGGP[i * Nx1 + j];
                        const double GXX2 = GGGP[i * Nx1 + j + 1];
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
                        const double SXY1 = SXY0[(i - 1)*Nx + j] * (1 - KXY1);
                        const double SXY2 = SXY0[i * Nx + j] * (1 - KXY2);
                        const double SXX1 = SXX0[i * Nx1 + j] * (1 - KXX1);
                        const double SXX2 = SXX0[i * Nx1 + j + 1] * (1 - KXX2);
                        // Density derivatives
                        const double dRHOdx = (RHOX[i * Nx1 + j + 1] - RHOX[i * Nx1 + j - 1]) / (2. * dx);
                        const double dRHOdy = (RHO[i * Nx + j] - RHO[(i - 1)*Nx + j]) / dy;
                        // Left part
                        const double temp = gx * dt * dRHOdy / 4.;
                        if (dx2==0 || dy2==0 || dt==0) {
                            std::cout << "Line 608 produces nans";
                        }
                        L_dense_gko->at(kx, kx) += -(ETAXX1 + ETAXX2) / dx2 - (ETAXY1 + ETAXY2) / dy2 - gx * dt * dRHOdx - inertia * RHOX[i * Nx1 + j]/dt;
                        L_dense_gko->at(kx, kx -Ny1 * Num_var) +=  ETAXX1 / dx2;
                        L_dense_gko->at(kx, kx -Ny1 * Num_var) += ETAXX1 / dx2;
                        L_dense_gko->at(kx, kx + Ny1 * Num_var) += ETAXX2 / dx2;
                        L_dense_gko->at(kx, kx - Num_var) += ETAXY1 / dy2;
                        L_dense_gko->at(kx, kx + Num_var) += ETAXY2 / dy2;
                        // vys1, vys2, vys3, vys4
                        L_dense_gko->at(kx, ky - Num_var) += (ETAXY1 - ETAXX1) / dx_dy - temp;
                        L_dense_gko->at(kx, ky) += (ETAXX1 - ETAXY2) / dx_dy - temp;
                        L_dense_gko->at(kx, ky - Num_var + Ny1 * Num_var) += (ETAXX2 - ETAXY1) / dx_dy - temp;
                        L_dense_gko->at(kx, ky + Ny1 * Num_var) += (ETAXY2 - ETAXX2) / dx_dy - temp;
                        // Pt1', Pt2'
                        L_dense_gko->at(kx, kp) += ptscale / dx;
                        L_dense_gko->at(kx, kp + Ny1 * Num_var) += -ptscale / dx;
                        // Right part
                        R[kx] = -RHOX[i * Nx1 + j] * (inertia * VX0[i * Nx1 + j] / dt + gx) - (SXX2 - SXX1) / dx - (SXY2 - SXY1) / dy;
                    }
                    
                    // 5b) Composing equation for vys
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vys = 0
                        if (i == Ny) {
                            L_dense_gko->at(ky, ky) += 1;
                        }
                            
                        // Left boundary
                        // Free Slip
                        if (j == 0) {
                            L_dense_gko->at(ky, ky) += 1;
                            L_dense_gko->at(ky, ky + Ny1 * Num_var) += 1;
                        }

                        // Right boundary
                        // Free Slip
                        if (j == Nx) {
                            L_dense_gko->at(ky, ky) += 1;
                            L_dense_gko->at(ky, ky - Ny1 * Num_var) += 1;
                        }

                        // Upper boundary: no penetration
                        if (i == 0 && j > 0 && j < Nx) {
                            L_dense_gko->at(ky, ky) += 1;
                        }

                        // Lower boundary: no penetration
                        if (i == Ny - 1 && j > 0 && j < Nx) {
                            L_dense_gko->at(ky, ky) +=  1;
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
                        double ETAXY1 = ETA[i * Nx + j - 1];
                        double ETAXY2 = ETA[i * Nx + j];
                        double ETAYY1 = ETAP[i * Nx1 + j];
                        double ETAYY2 = ETAP[(i + 1) * Nx1 + j];
                        // Shear modulus
                        const double GXY1 = GGG[i * Nx + j - 1];
                        const double GXY2 = GGG[i * Nx + j];
                        const double GYY1 = GGGP[i * Nx1 + j];
                        const double GYY2 = GGGP[(i + 1) * Nx1 + j];
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
                        const double SXY1 = SXY0[i * Nx + j - 1] * (1 - KXY1);
                        const double SXY2 = SXY0[i * Nx + j] * (1 - KXY2);
                        const double SYY1 = SYY0[i * Nx1 + j] * (1 - KYY1);
                        const double SYY2 = SYY0[(i + 1) * Nx1 + j] * (1 - KYY2);
                        // Density derivatives
                        const double dRHOdy = (RHOY[(i + 1) * Nx1 + j] - RHOY[(i - 1) * Nx1 + j]) / 2. / dy;
                        const double dRHOdx = (RHO[i * Nx + j] - RHO[i * Nx + j - 1]) / dx;
                        // Left part
                        const double temp = gy * dt * dRHOdx / 4.;
                        // vys3, vys1, vys5, vys2, vys4
                        L_dense_gko->at(ky, ky) += -(ETAYY1 + ETAYY2) / dy2 - (ETAXY1 + ETAXY2) / dx2 - gy * dt * dRHOdy - inertia * RHOY[i * Nx1 + j] / dt;
                        L_dense_gko->at(ky, ky - Ny1 * Num_var) += ETAXY1 / dx2;
                        L_dense_gko->at(ky, ky + Ny1 * Num_var) += ETAXY2 / dx2;
                        L_dense_gko->at(ky, ky - Num_var) += ETAYY1 / dy2;
                        L_dense_gko->at(ky, ky + Num_var) += ETAYY2 / dy2;
                        // vxs1, vxs2, vxs3, vxs4
                        L_dense_gko->at(ky, kx - Ny1 * Num_var) += (ETAXY1 - ETAYY1) / dx_dy - temp;
                        L_dense_gko->at(ky, kx + Num_var - Ny1 * Num_var) += (ETAYY2 - ETAXY1) / dx_dy - temp;
                        L_dense_gko->at(ky, kx) += (ETAYY1 - ETAXY2) / dx_dy - temp;
                        L_dense_gko->at(ky, kx + Num_var) += (ETAXY2 - ETAYY2) / dx_dy - temp;
                        // Pt1', Pt2'
                        L_dense_gko->at(ky, kp) += ptscale / dy;
                        L_dense_gko->at(ky, kp + Num_var) += -ptscale / dy;
                        // Right part
                        R[ky] = -RHOY[i * Nx1 + j] * (inertia * VY0[i * Nx1 + j] / dt + gy) - (SYY2 - SYY1) / dy - (SXY2 - SXY1) / dx;
                    }

                    // 5c) Composing equation for Pt
                    if (i == 0 || j == 0 || i == Ny || j == Nx) { // || (i == 2 && j == 2))
                        // BC equation: 1 * Pt = 0
                        L_dense_gko->at(kp, kp) += 1;
                    } else {
                        // Solid Continuity: dVxs / dx + dVys / dy + (Pt - Pf) / ETAbulk = 0
                        //              vys1
                        //               |
                        //        vxs1 -- Pt, Pf -- vxs2
                        //               |
                        //              vys2
                        // Drained compressibility
                        const double BETADRAINED = (1. / GGGB[i * Nx1 + j] + BETASOLID) / (1 - POR[i * Nx1 + j]);
                        // Biott - Willis koefficient
                        const double KBW = 1 - BETASOLID / BETADRAINED;
                        // Left part
                        // vxs1, vxs2
                        L_dense_gko->at(kp, kx - Ny1 * Num_var) += -1. / dx;
                        L_dense_gko->at(kp, kx) += 1. / dx;
                        // vys1, vys2
                        L_dense_gko->at(kp, ky - Num_var) += -1. / dy;
                        L_dense_gko->at(kp, ky) += 1. / dy;
                        // Pt, Pf
                        L_dense_gko->at(kp, kp) += ptscale * (1. / ETAB[i * Nx1 + j] / (1 - POR[i * Nx1 + j]) + BETADRAINED / dt);
                        L_dense_gko->at(kp, kpf) += -pfscale * (1. / ETAB[i * Nx1 + j] / (1 - POR[i * Nx1 + j]) + BETADRAINED * KBW / dt);
                        // Right part
                        R[kp] = BETADRAINED * (PT0[i * Nx1 + j] - KBW * PF0[i * Nx1 + j]) / dt + DILP[i * Nx1 + j];
                    }


                    // 5d) Composing equation for vxD
                    if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (j == Nx) {
                            L_dense_gko->at(kxf, kxf) += 1;
                        }

                        // Upper boundary: symmetry
                        if (i == 0 && j < Nx) {
                            L_dense_gko->at(kxf, kxf) += 1;
                            L_dense_gko->at(kxf, kxf + Num_var) += -1;
                        }

                        // Lower boundary: symmetry
                        if (i == Ny && j < Nx) {
                            L_dense_gko->at(kxf, kxf) += 1;
                            L_dense_gko->at(kxf, kxf - Num_var) += -1;
                        }

                        // Left boundary
                        // no penetration
                        if (j == 0) {
                            L_dense_gko->at(kxf, kxf) += 1;
                        }

                        // Right boundary
                        // no penetration
                        if (j == Nx - 1) {
                            L_dense_gko->at(kxf, kxf) += 1;
                        }

                    } else {
                        // Fluid X - Darsi: - ETAfluid / K * VxD - dPf / dx = -RHOf * gx + RHOf * DVxs / Dt
                        //
                        //  Pf1 --- vxD, vxs --- Pf2
                        //
                        // Left part
                        // vxD, vxs
                        L_dense_gko->at(kxf, kxf) += -ETADX[i * Nx1 + j] - RHOFX[i * Nx1 + j] / PORX[i * Nx1 + j] * inertia / dt;
                        L_dense_gko->at(kxf, kx) += -RHOFX[i * Nx1 + j] * inertia / dt;
                        // Pf1', Pf2'
                        L_dense_gko->at(kxf, kpf) += pfscale / dx;
                        L_dense_gko->at(kxf, kpf + Ny1 * Num_var) += -pfscale / dx;
                        // Right part
                        R[kxf] = -RHOFX[i * Nx1 + j] * (inertia * VXF0[i * Nx1 + j] / dt + gx);
                    }


                    // 5e) Composing equation for vyD
                    if (j == 0 || j == Nx || i == 0 || i >= Ny - 1) {
                        // Ghost nodes: 1 * vxs = 0
                        if (i == Ny) {
                            L_dense_gko->at(kyf, kyf) += 1;
                        }

                        // Left boundary
                        // symmetry
                        if (j == 0 && i > 0 && i < Ny - 1) {
                            L_dense_gko->at(kyf, kyf) += 1;
                            L_dense_gko->at(kyf, kyf + Ny1 * Num_var) += -1;
                        }

                        // Right boundary
                        // symmetry
                        if (j == Nx && i > 0 && i < Ny - 1) {
                            L_dense_gko->at(kyf, kyf) += 1;
                            L_dense_gko->at(kyf, kyf - Ny1 * Num_var) += -1;
                        }

                        // Upper boundary: no penetration
                        if (i == 0) {
                            L_dense_gko->at(kyf, kyf) += 1;
                            R[kyf] = bcvyflower;
                        }

                        // Lower boundary: no penetration
                        if (i == Ny - 1) {
                            L_dense_gko->at(kyf, kyf) += 1;
                            R[kyf] = bcvyflower;
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
                        // vyD, vys
                        L_dense_gko->at(kyf, kyf) += -ETADY[i * Nx1 + j] - RHOFY[i * Nx1 + j] / PORY[i * Nx1 + j] * inertia / dt;
                        L_dense_gko->at(kyf, ky) += -RHOFY[i * Nx1 + j] * inertia / dt;
                        // Pf1', Pf2'
                        L_dense_gko->at(kyf, kpf) += pfscale / dy;
                        L_dense_gko->at(kyf, kpf + Num_var) += -pfscale / dy;
                        // Right part
                        R[kyf] = -RHOFY[i * Nx1 + j] * (inertia * VYF0[i * Nx1 + j] / dt + gy);
                    }



                    // 5f) Composing equation for Pf
                    if (j == 0 || j == Nx || i <= 1 || i >= Ny - 1) {
                        // BC equation: 1 * Pf = 0
                        // Real BC
                        if (i == 1 || i == Ny - 1) {
                            L_dense_gko->at(kpf, kpf) += pfscale;
                            L_dense_gko->at(kpf, kp) += -ptscale;
                            R[kpf] = -PTFDIFF;
                        } else {
                            L_dense_gko->at(kpf, kpf) += 1;
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
                        const double BETADRAINED = (1 / GGGB[i * Nx1 + j] + BETASOLID) / (1 - POR[i * Nx1 + j]);

                        // Biott - Willis koefficient
                        const double KBW = 1 - BETASOLID / BETADRAINED;
                        // Skempton koefficient
                        const double KSK = (BETADRAINED - BETASOLID) / (BETADRAINED - BETASOLID + POR[i * Nx1 + j] * (BETAFLUID - BETASOLID));
                        // Left part
                        // vxs1, vxs2
                        L_dense_gko->at(kpf, kxf - Ny1 * Num_var) += -1. / dx;
                        L_dense_gko->at(kpf, kxf) += 1. / dx;
                        // vys1, vys2
                        L_dense_gko->at(kpf, kyf - Num_var) += -1. / dy;
                        L_dense_gko->at(kpf, kyf) +=  1. / dy;
                        // Pt, Pf
                        L_dense_gko->at(kpf, kp) += -ptscale * (1 / ETAB[i * Nx1 + j] / (1 - POR[i * Nx1 + j]) + BETADRAINED * KBW / dt);
                        L_dense_gko->at(kpf, kpf) +=  pfscale * (1 / ETAB[i * Nx1 + j] / (1 - POR[i * Nx1 + j]) + BETADRAINED * KBW / KSK / dt);
                        // Right part
                        R[kpf] = -BETADRAINED * KBW * (PT0[i * Nx1 + j] - 1. / KSK * PF0[i * Nx1 + j]) / dt - DILP[i * Nx1 + j];
                    }

                    if (antiplane) {
                        // 5o) Composing equation for vzs (out-of-plane component)
                        if (i == 0 || i == Ny || j == 0 || j >= Nx - 1) {
                            // Ghost nodes: 1*vzs=0
                            if (j == Nx1) {
                                L_dense_gko->at(kz, kz) += 1;
                            }

                            // Upper boundary
                            if (i == 1 && j < Nx1) {
                                L_dense_gko->at(kz, kz) += 1;
                                L_dense_gko->at(kz, kz + Num_var) += 1;
                                R[kz] = 2 * bclower;
                            }

                            // Lower boundary
                            if (i == Ny1 && j < Nx1) {
                                L_dense_gko->at(kz, kz) += 1;
                                L_dense_gko->at(kz, kz - Num_var) += 1;
                                R[kz] = -2 * bclower;
                            }

                            // Left boundary
                            if (j == 1 && i > 1 && i < Ny1) {
                                L_dense_gko->at(kz, kz) += 1;
                                L_dense_gko->at(kz, kz + Num_var * Ny1) += -1;
                            }

                            // Right boundary
                            if (j == Nx && i > 1 && i < Ny1) {
                                L_dense_gko->at(kz, kz) += 1;
                                L_dense_gko->at(kz, kz - Num_var * Ny1) += -1;
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

                            double ETAXY = ETA[i * Nx + j];
                            // Shear modulus
                            const double GXY = GGG[i * Nx + j];
                            // Viscoelasticity factor
                            const double KXY = dt * GXY / (dt * GXY + ETAXY);
                            // Numerical viscosity
                            ETAXY = ETAXY * KXY;
                            // Numerical stresses
                            const double SZY1 = SZY0[(i - 1) * Nx + j] * (1 - KXY);
                            const double SZY2 = SZY0[i * Nx + j] * (1 - KXY);
                            const double SZX1 = SZX0[i * Nx + j -1] * (1 - KXY);
                            const double SZX2 = SZX0[i * Nx + j] * (1 - KXY);
                            // Left part
                            // vzs3, vzs3, vzs3, vzs2, vzs4
                            L_dense_gko->at(kz, kz) += -2 * (ETAXY / dx2 + ETAXY / dy2) - RHOX[i * Nx1 + j] / dt;
                            L_dense_gko->at(kz, kz - Ny1 * Num_var) += ETAXY / dx2;
                            L_dense_gko->at(kz, kz + Ny1 * Num_var) += ETAXY / dx2;
                            L_dense_gko->at(kz, kz - Num_var) += ETAXY / dy2;
                            L_dense_gko->at(kz, kz + Num_var) += ETAXY / dy2;
                            // Right part
                            R[kz] = -RHOX[i * Nx1 + j] * (VZ0[i * Nx1 + j] / dt) - (SZX2 - SZX1) / dx - (SZY2 - SZY1) / dy;
                        }
                    }
                }
            }
            std::cout << "Matrix Assembled" << std::endl;
            std::cout << "The L matrix looks like this:\n" << N;
            int counter = 0;
            int nancounter=0;
            for (int l = 0; l<N*N; l++) {
                if (L_dense[l]!=0) {
                    if (counter < 0) {
                        std::cout << L_dense[l] << "\n";
                    }
                    counter++;
                }
                if (isnan(L_dense[l])) {
                    nancounter++;
                }
            }
            std::cout << "All in all there are " <<counter<< " nonzero entries!\nAlso there are " << nancounter << " nans :(" << std::endl;

            // Copy the data from the Dense matrix to the CSR matrix
            L->copy_from(L_dense_gko.get());

            /*
                        double* L_print = L->get_values();
                        std::cout << "Number of stored elements: " << L->get_num_stored_elements() << std::endl;
                        std::cout << "The matrix is " << N << "x" << N << "shaped" << std::endl;
                        int counter = 0;
                        for (int i=0; i<L->get_num_stored_elements(); i++) {
                            if (abs(L_print[i])<1.0e-100) {
                                counter++;
                            }
                            else {
                                std::cout << L_print[i] << " ";
                            }
                        }
                        std::cout << "/n The number of elements<e-100 is: " << counter << std::endl;
                        std::cout << "This means a total of " << L->get_num_stored_elements()-counter << " nonzero elements are in the matrix!" << std::endl;
            */
            double reduction_factor = 1e-10;
            std::cout << "A reduction factor of " << reduction_factor << " is used for the solver!" << std::endl;
            // Build the solver from the CSR matrix
            /*
            auto solver = gko::solver::Cg<>::build()
                //.with_preconditioner(gko::preconditioner::Ic<>::build().on(exec))
                .with_criteria(gko::stop::ResidualNorm<>::build()
                    .with_baseline(gko::stop::mode::rhs_norm)
                    .with_reduction_factor(reduction_factor)
                    .on(exec))
                .on(exec)
                ->generate(L);*/

            //Alternative solver: Careful, it is experimental and I will have to check out what that means:
            auto solver = gko::experimental::solver::Direct<double, int>::build()
                .with_factorization(gko::experimental::factorization::Lu<double, int>::build())
                .on(exec)
                ->generate(L);

            std::cout << "Solver built" << std::endl;

            /*            std::cout << "Printing the R vector" << std::endl;

                        for (int i = 0; i<N; i++) {
                            std::cout << R_gko->at(i,0) << " ";
                        }*/

            //Saving the dense matrix for easier viewing:
            if (save_first_matrix) {
                save_first_matrix = false;
                std::string L_matrix_data = "./output_data/L_matrix.txt";
                std::ofstream reference_output;
                reference_output.open(L_matrix_data, std::ios_base::app | std::ios_base::out);
                for (std::size_t l = 0; l < N*N; ++l) {
                    if (l%N==0) {
                        reference_output << "\n";
                    }
                    reference_output << L_dense[l] << " ";
                }
                reference_output << std::endl;
                reference_output.close();
            }
            //Saving the R vector and the L matrix to a file:
            if (false){
                save_first_matrix = false;
                std::string L_matrix_data = "./output_data/L_matrix.txt";
                std::ofstream reference_output;
                double* L_ptr = L->get_values();
                reference_output.open(L_matrix_data, std::ios_base::app | std::ios_base::out);
                for (std::size_t i = 0; i < L->get_num_stored_elements(); ++i) {
                    reference_output << L_ptr[i] << " ";
                }
                reference_output << std::endl;
                reference_output.close();


                //Saving the R vector and the L matrix to a file:
                std::string R_matrix_data = "./output_data/R_vector.txt";
                std::ofstream reference_output_R;
                reference_output_R.open(R_matrix_data, std::ios_base::app | std::ios_base::out);
                for (std::size_t i = 0; i < N; ++i) {
                    reference_output_R << R[i] << " ";
                }
                reference_output_R << std::endl;
                reference_output_R.close();

                std::cout << "Files generated!" << std::endl;
            }

            std::cout << "Apply uses initial_guess is " << solver->apply_uses_initial_guess() << std::endl;
            solver->apply(R_gko, S_gko);
            std::cout << "Solved successfully for timestep=" << timestep << std::endl;

            // slightly slower in parallel
            // #pragma omp parallel for collapse(2)
            for (int j = 0; j < Nx1; j++) {
                for (int i = 0; i < Ny1; i++) {
                    // Global indexes for vx, vy, P
                    const int kp    = (j * Ny1 + i) * Num_var;
                    // Reload solution
                    pt[i * Nx1 + j]        = S[kp] * ptscale;
                    vxs[i * Nx1 + j]       = S[kp + 1];
                    vys[i * Nx1 + j]       = S[kp + 2];
                    pf[i * Nx1 + j]        = S[kp + 3] * pfscale;
                    vxD[i * Nx1 + j]       = S[kp + 4];
                    vyD[i * Nx1 + j]       = S[kp + 5];
                    if (antiplane) {
                        vzs[i * Nx1 + j]   = S[kp + 6];
                    }
                }
            }

            Vmax = findmax(VSLIPB, Nx, Ny);

            /*
                        if (dt > 1e4 && Vmax < 1e-7) {
                            double avgpt = pt.sum() / (double)(pt.rows() * pt.cols()); //calculate average total pressure
                            double diffpt = (PCONF + PTFDIFF) - avgpt;
                            pt += Eigen::MatrixXd::Constant(Ny1, Nx1, diffpt);
                        }
                        */

            // Velocity change
            if (antiplane) {
                for (int i=0; i<Ny; i++) {
                    vzs_gko->at(i,Nx) = vzs_gko->at(i, Nx-1);
                }
                DVZ0_gko->copy_from(vzs_gko.get());
                DVZ0_gko->sub_scaled(sub_alpha, VZ0_gko);
            }

            DVX0_gko->copy_from(vxs_gko.get());
            DVX0_gko->sub_scaled(sub_alpha, VX0_gko);
            DVY0_gko->copy_from(vys_gko.get());
            DVY0_gko->sub_scaled(sub_alpha, VY0_gko);

            // Define timestep
            bool yn = false;


            // Plastic iterations
            // Compute strain rate, stress and stress change
            // Process internal basic nodes
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    EXY[i * Nx + j]           = .5 * ((vxs[(i + 1) * Nx1 + j] - vxs[i * Nx1 + j]) / dy + (vys[i * Nx1 + j + 1] - vys[i * Nx1 + j]) / dx);
                    const double KXY    = divplus(dt * GGG[i * Nx + j], ETA[i * Nx + j]);
                    SXY[i * Nx + j]           = 2 * ETA[i * Nx + j] * EXY[i * Nx + j] * KXY + SXY0[i * Nx + j] * (1 - KXY);

                    if (antiplane) {
                        EZX[i * Nx + j] = .5 * ((vzs[i * Nx1 + j + 1] - vzs[i * Nx1 + j]) / dx + (vxs[(i + 1) * Nx1 + j] - vxs[i * Nx1 + j]) / dx);
                        EZY[i * Nx + j] = .5 * ((vzs[(i + 1) * Nx1 + j] - vzs[i * Nx1 + j]) / dy + (vys[i * Nx1 + j + 1] - vys[i * Nx1 + j]) / dy);
                        SZX[i * Nx + j] = 2 * ETA[i * Nx + j] * EZX[i * Nx + j] * KXY + SZX0[i * Nx + j] * (1 - KXY);
                        SZY[i * Nx + j] = 2 * ETA[i * Nx + j] * EZY[i * Nx + j] * KXY + SZY0[i * Nx + j] * (1 - KXY);
                    }
                }
            }

            // Process pressure cells
            for (int i = 1; i < Ny; i++) {
                for (int j = 1; j < Nx; j++) {
                    // EXX, SXX
                    EXX[i * Nx1 + j]           = (2 * (vxs[i * Nx1 + j] - vxs[i * Nx1 + j - 1]) / dx - (vys[i * Nx1 + j] - vys[(i - 1) * Nx1 + j]) / dy) / 3.;
                    EYY[i * Nx1 + j]           = (2 * (vys[i * Nx1 + j] - vys[(i - 1) * Nx1 + j]) / dy - (vxs[i * Nx1 + j] - vxs[i * Nx1 + j - 1]) / dx) / 3.;
                    const double KXX    = divplus(dt * GGGP[i * Nx1 + j], ETAP[i * Nx1 + j]);
                    SXX[i * Nx1 + j]           = 2 * ETAP[i * Nx1 + j] * EXX[i * Nx1 + j] * KXX + SXX0[i * Nx1 + j] * (1 - KXX);
                    SYY[i * Nx1 + j]           = 2 * ETAP[i * Nx1 + j] * EYY[i * Nx1 + j] * KXX + SYY0[i * Nx1 + j] * (1 - KXX);
                }
            }




            // This can be implemented as for loop, change it later on:
            // External P - nodes: symmetry
            copy_bounds(pt_gko, Nx, Ny);
            copy_bounds(pf_gko, Nx, Ny);
            copy_bounds(EXX_gko, Nx, Ny);
            copy_bounds(SXX_gko, Nx, Ny);
            copy_bounds(SXX0_gko, Nx, Ny);
            copy_bounds(EYY_gko, Nx, Ny);
            copy_bounds(SYY_gko, Nx, Ny);
            copy_bounds(SYY0_gko, Nx, Ny);
            copy_bounds(ETAP_gko, Nx, Ny);
            copy_bounds(ETAB_gko, Nx, Ny);
            copy_bounds(GGGP_gko, Nx, Ny);
            copy_bounds(GGGB_gko, Nx, Ny);

            if (antiplane) {
                for (int i=0; i<Ny; i++) {
                    SZX_gko->at(i, Nx - 1) = SZX_gko->at(i, Nx - 2);
                    EZX_gko->at(i, Nx - 1) = EZX_gko->at(i, Nx - 2);
                }



                // Compute stress and strain rate invariants and dissipation
                // Process pressure cells
                for (int i = 1; i < Ny; i++) {
                    for (int j = 1; j < Nx; j++) {
                        // EXY term is averaged from four surrounding basic nodes
                        EII[i * Nx1 + j] = sqrt(.5 * (pow(EXX[i * Nx1 + j], 2) + pow(EYY[i * Nx1 + j], 2)) +
                                        (square_block(EXY, i - 1, j - 1, Nx) +
                                        square_block(EZX, i - 1, j - 1, Nx) +
                                        square_block(EZY, i - 1, j - 1, Nx)) / 4.
                                    );

                        // Second strain rate invariant SII
                        // SXY term is averaged from four surrounding basic nodes
                        SII[i * Nx1 + j] = sqrt(.5 * (pow(SXX[i * Nx1 + j], 2) + pow(SYY[i * Nx1 + j], 2)) +
                                        square_block(SXY, i - 1, j - 1, Nx) +
                                        square_block(SZX, i - 1, j - 1, Nx) +
                                        square_block(SZY, i - 1, j - 1, Nx)) / 4.;
                    }
                }
            }


            // Update viscosity for yielding
            // dt0 = dt; dt = dt * 1.1;
            ETA5_gko->copy_from(ETA0_gko.get());
            // Basic nodes


            GkoSetZero(DSY_gko);
            GkoSetZero(YNY_gko);
            GkoSetZero(SigmaY_gko);
            GkoSetZero(SII_fault_gko);


            int ynpl    = 0;
            double ddd  = 0;
            dtlapusta   = 1e7;
            OM5_gko->copy_from(OM_gko.get());

            // Power law plasticity model Yi et al., 2018
            // Journal of Offshore Mechanics and Arctic Engineering
            double dtslip = 1e30;
            // power_law(timestep, ynpl, ddd);




            if (timestep > tyield) {
            for (int i = 0; i < Ny; i++) {
                // if (y[i]  >= upper_block && y[i]  <= lower_block) {
                if (i == line_fault) {
                    for (int j = 0; j < Nx; j++) {
                        // reducing matrix calls
                        const double arsf_temp = ARSF[i * Nx + j];
                        const double brsf_temp = BRSF[i * Nx + j];
                        const double lrsf_temp = LRSF[i * Nx + j];
                        const double fric_temp = FRIC[i * Nx + j];
                        const double eta0_temp = ETA0[i * Nx + j];


                        const double sxx_temp = BlockSum(SXX, i, j, 2, 2, Nx1, Ny1);
                        const double syy_temp = BlockSum(SYY, i, j, 2, 2, Nx1, Ny1);

                        // SXX, pt are averaged from four surrounding pressure nodes
                        if (antiplane) {
                            SIIB[i * Nx + j] = sqrt(.5 * (pow(SXX[i * Nx1 + j], 2) +
                                                    pow(SYY[i * Nx1 + j], 2)) +
                                                    pow(SXY[i * Nx + j], 2) +
                                                    pow(SZX[i * Nx + j], 2) +
                                                    pow(SZY[i * Nx + j], 2)
                                                );
                        } else {
                            SIIB[i * Nx + j]    = sqrt(pow(SXY[i * Nx + j], 2) +
                                                (pow(sxx_temp, 2) +
                                                pow(syy_temp, 2) +
                                                pow(sxx_temp + syy_temp, 2))/ 32.
                                            );
                        }
                        // Computing "elastic" stress invariant
                        const double siiel = SIIB[i * Nx + j] / (ETA[i * Nx + j] / (GGG[i * Nx + j] * dt + ETA[i * Nx + j]));

                        // Compute old viscoplastic slip rate
                        // Compute PEFF
                        const double prB_temp = (BlockSum(pt, i, j, 2, 2, Nx1, Ny1) - BlockSum(pt, i, j, 2, 2, Nx1, Ny1))/4.0;
                        double prB = prB_temp;

                        if (prB < 1e3) {
                            prB = 1e3;
                        }
                        // Compute old power law strain rate
                        double SIIB1 = SIIB[i * Nx + j];

                        // Compute slip velocity for current stress invariant and state
                        double EIISLIP = V0 * sinh(max(SIIB1, 0.) / arsf_temp / prB) * exp(-(brsf_temp * OM[i * Nx + j] + fric_temp) / arsf_temp) / dx;

                        // Compute new ETAVP
                        double ETAPL        = SIIB1 / 2. / EIISLIP;
                        double ETAVP        = 1. / (1. / eta0_temp + 1. / ETAPL);
                        // Compute new stress invariant
                        double SIIB2        = siiel * ETAVP / (GGG[i * Nx + j] * dt + ETAVP);
                        const double DSIIB1 = SIIB2 - SIIB1;

                        // Compute slip velocity for current stress invariant and state
                        double V            = 2 * V0 * sinh(max(SIIB2, 0.) / arsf_temp / prB) *
                                            exp(-(brsf_temp * OM[i * Nx + j] + fric_temp) / arsf_temp);

                        EIISLIP = V / dx / 2.;

                        // Compute new ETAVP
                        ETAPL               = SIIB2 / 2. / EIISLIP;
                        ETAVP               = 1. / (1. / eta0_temp + 1. / ETAPL);
                        // Compute new stress invariant
                        const double DSIIB2 = siiel * ETAVP / (GGG[i * Nx + j] * dt + ETAVP) - SIIB2;
                        double SIIB4        = 0.;

                        if ((DSIIB1 >= 0 && DSIIB2 <= 0) || (DSIIB1 <= 0 && DSIIB2 >= 0)) {
                            double DSIIB = 1e9;

                            while(abs(DSIIB) > 1e-3) {
                                SIIB4   = (SIIB1 + SIIB2) / 2.;

                                //Compute slip velocity for current stress invariant and state
                                V       = 2 * V0 * sinh(max((SIIB4), 0.) / arsf_temp / prB) *
                                        exp( - (brsf_temp * OM[i * Nx + j] + fric_temp) / arsf_temp);

                                EIISLIP = V / dx / 2.;

                                // Compute new ETAVP
                                ETAPL   = SIIB4 / 2. / EIISLIP;
                                ETAVP   = 1. / (1. / eta0_temp + 1. / ETAPL);
                                // Compute new stress invariant
                                DSIIB   = siiel * ETAVP / (GGG[i * Nx + j] * dt + ETAVP) - SIIB4;
                                if ((DSIIB >= 0 && DSIIB1 >= 0) || (DSIIB <= 0 && DSIIB1 <= 0)) {
                                    SIIB1 = SIIB4;
                                } else {
                                    SIIB2 = SIIB4;
                                }
                            }
                        }

                        if (V * dt / lrsf_temp > 1e-6) {
                            OM5[i * Nx + j] = log(V0 / V + (exp(OM0[i * Nx + j]) - V0 / V) * exp( - V * dt / lrsf_temp));
                        } else {
                            OM5[i * Nx + j] = log(exp(OM0[i * Nx + j]) * (1 - V * dt / lrsf_temp) + V0 * dt / lrsf_temp);
                        }

                        // Compute yielding stress
                        const double syield = max(syieldmin, prB_temp * arsf_temp * asinh(V / 2. / V0 * exp((brsf_temp * OM5[i * Nx + j] + fric_temp) / arsf_temp)));

                        // Compute visco - plastic viscosity
                        double etapl = eta0_temp * syield / (eta0_temp * V + syield);

                        // Save syield
                        SigmaY[i * Nx + j]        = syield;
                        VSLIPB[i * Nx + j]        = V;
                        SII_fault[i * Nx + j]     = SIIB4;

                        // reduces calls on matrix
                        const double g_temp = 2 * GGG[i * Nx + j]; // reduces calls on matrix

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
                        if (YNY0[i * Nx + j] > 0) {
                            ynn         = true;
                            DSY[i * Nx + j]   = SIIB[i * Nx + j] - syield;
                            ddd         += pow(DSY[i * Nx + j], 2);
                            ynpl++;
                        }

                        // Update viscosity
                        const double A = syield / siiel;
                        if (A < 1) {
                            // New viscosity for the basic node
                            etapl = dt * GGG[i * Nx + j] * A / (1 - A);
                            if (etapl < eta0_temp) {
                                // Update plastic nodes
                                ETA5[i * Nx + j]      = pow(etapl, 1 - etawt) * pow(ETA[i * Nx + j], etawt);
                                YNY[i * Nx + j]       = 1;
                                // Count yelding nodes
                                if (!ynn) {
                                    DSY[i * Nx + j]   = SIIB[i * Nx + j] - syield;
                                    ddd         += pow(DSY[i * Nx + j], 2);
                                    ynpl++;
                                }
                            } else {
                                ETA5[i * Nx + j] = eta0_temp;
                            }
                        } else {
                            ETA5[i * Nx + j] = eta0_temp;
                        }
                    }
                }
            }
        }



        // Compute Error
        DSYLSQ[iterstep] = 0;
        if (ynpl > 0) {
            DSYLSQ[iterstep] = sqrt(ddd / ynpl);
        }
        if (ynpl == 0) {
            ETA = ETA0;
        }

        // connot calculate DSYLSQ(iterstep - 1) if iterstep = 0
        // need to know what to do when iterstep = 0
        const double D_iter = DSYLSQ[iterstep];
        double D_iter_quot = 0;
        if (iterstep != 0) {
            D_iter_quot = D_iter / DSYLSQ[iterstep - 1];
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
            double maxvxs = max(findmax(vxs, Nx1, Ny1), abs(findmin(vxs, Nx1, Ny1)));
            double maxvys = max(findmax(vys, Nx1, Ny1), abs(findmin(vys, Nx1, Ny1)));
            double maxvzs;
            if (antiplane) {
                double maxvzs = max(findmax(vzs, Nx1, Ny1), abs(findmin(vzs, Nx1, Ny1)));
                maxvxy = sqrt(pow(findmax(vxs, Nx1, Ny1) - findmin(vxs, Nx1, Ny1), 2) + pow(findmax(vys, Nx1, Ny1) - findmin(vys, Nx1, Ny1), 2) + pow(findmax(vzs, Nx1, Ny1) - findmin(vzs, Nx1, Ny1), 2));
            } else {
                maxvxy = sqrt(pow(findmax(vxs, Nx1, Ny1) - findmin(vxs, Nx1, Ny1), 2) + pow(findmax(vys, Nx1, Ny1) - findmin(vys, Nx1, Ny1), 2));
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
                    if (yvx[i] >= upper_block && yvx[i] <= lower_block) {
                        maxvxs = max(maxvxs, abs(vxs[i * Nx1 + j]));
                    }
                    if (yvy[i] >= upper_block && yvy[i] <= lower_block) {
                        maxvys = max(maxvys, abs(vys[i * Nx1 + j]));
                    }
                    if (antiplane && yvy[i] >= upper_block && yvy[i] <= lower_block) {
                        maxvzs = max(maxvzs, abs(vzs[i * Nx1 + j]));
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
                    if (VSLIPB[i * Nx + j] > 0) {
                        dtslip = min(dtslip, dx * stpmax / VSLIPB[i * Nx + j]);
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
            if (!yn && (ynpl == 0 || (DSYLSQ[iterstep] < errmin && iterstep > 0 && abs(vratio) < vratiomax))) {
                ynstop = true;
            } else {
                // Recomputing ETA
                for (int i = 0; i < Ny; i++) {
                    for (int j = 0; j < Nx; j++) {
                        ETA[i * Nx + j] = max(min(ETA5[i * Nx + j], ETA0[i * Nx + j]), etamin);
                    }
                }
                // Save current viscosity

                ETA50_gko->copy_from(ETA_gko.get());
                YNY0_gko->copy_from(YNY_gko.get());
                OM_gko->copy_from(OM5_gko.get());

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
        ETA00_gko->copy_from(ETA50_gko.get());
        OM0_gko->copy_from(OM_gko.get());

        ETA00 = ETA50;
        OM0 = OM;
        // Recheck displacement timestep
        dt = min(min(dtx, dty), min(dt, dtlapusta));

        // This can be done in a loop instead
        // Compute strain rate, stress and stress change
        GkoSetZero(ESP_gko);
        GkoSetZero(EXY_gko);
        GkoSetZero(SXY_gko);
        GkoSetZero(EXX_gko);
        GkoSetZero(SXX_gko);
        GkoSetZero(EYY_gko);
        GkoSetZero(SYY_gko);
        GkoSetZero(EII_gko);
        GkoSetZero(EIIVP_gko);
        GkoSetZero(SII_gko);
        GkoSetZero(DSII_gko);
        GkoSetZero(DIS_gko);
                
                
            // Process internal basic nodes
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    // ESP = .5 *(dVy / dx - dVx / dy), EXY, SXY
                    ESP[i * Nx + j] = .5 * ((vys[i * Nx1 + j + 1] - vys[i * Nx1 + j]) / dx - (vxs[(i + 1) * Nx1 + j] - vxs[i * Nx1 + j]) / dy);
                    EXY[i * Nx + j] = .5 * ((vxs[(i + 1) * Nx1 + j] - vxs[i * Nx1 + j]) / dy + (vys[i * Nx1 + j + 1] - vys[i * Nx1 + j]) / dx);
                    const double KXY = divplus(dt * GGG[i * Nx + j], ETA[i * Nx + j]);
                    SXY[i * Nx + j] = 2 * ETA[i * Nx + j] * EXY[i * Nx + j] * KXY + SXY0[i * Nx + j] * (1 - KXY);

                    if (antiplane) {
                        EZX[i * Nx + j] = .5 *((vzs[i * Nx1 + j + 1] - vzs[i * Nx1 + j]) / dx + (vxs[(i + 1) * Nx1 + j] - vxs[i * Nx1 + j]) / dx);
                        EZY[i * Nx + j] = .5 *((vzs[(i + 1) * Nx1 + j] - vzs[i * Nx1 + j]) / dy + (vys[i * Nx1 + j + 1] - vys[i * Nx1 + j]) / dy);
                        SZX[i * Nx + j] = 2 * ETA[i * Nx + j] * EZX[i * Nx + j] * KXY + SZX0[i * Nx + j] * (1 - KXY);
                        SZY[i * Nx + j] = 2 * ETA[i * Nx + j] * EZY[i * Nx + j] * KXY + SZY0[i * Nx + j] * (1 - KXY);
                    }
                }
            }
            
            
            // Process pressure cells
            process_p_cells();

            // /////////////////////////////////////////////////////////////////////////////////////// 
            // Move markers by nodal velocity field
            // /////////////////////////////////////////////////////////////////////////////////////// 
            // This line is disabled also in the move_markers() function by setting the for loop to   for (int m=0;m<0;m++)
            //move_markers();

            const int temp = timestep - 1;
        
            // Update timesum
            timesum += dt;
            timesumcur[temp] = timesum;
            dtcur[temp] = dt;
        
            maxvxsmod[temp] = -1e30;
            minvxsmod[temp] = 1e30;
            maxvysmod[temp] = -1e30;
            minvysmod[temp] = 1e30;
            maxvzsmod[temp] = -1e30;
            minvzsmod[temp] = 1e30;

            VX0 = vxs;
            VY0 = vys;
            if (antiplane) {
                VZ0 = vzs;
            }
                
                
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    // Vx
                    if (RHOX[i * Nx1 + j] > 2000 && i != 0) {
                        maxvxsmod[temp] = max(maxvxsmod[temp], vxs[i * Nx1 + j]);
                        minvxsmod[temp] = min(minvxsmod[temp], vxs[i * Nx1 + j]);
                    }
                    // Vy
                    if (RHOY[i * Nx1 + j] > 2000 && j != 0) {
                        maxvysmod[temp] = max(maxvysmod[temp], vys[i * Nx1 + j]);
                        minvysmod[temp] = min(minvysmod[temp], vys[i * Nx1 + j]);
                    }
                    // Vz
                    if (antiplane && RHOY[i * Nx1 + j] > 2000) {
                        maxvzsmod[temp] = max(maxvzsmod[temp], vzs[i * Nx1 + j]);
                        minvzsmod[temp] = min(minvzsmod[temp], vzs[i * Nx1 + j]);
                    }
                }
            }



            for (int i = 0; i < Ny1; i++) {
                for (int j = 0; j < Nx1; j++) {
                    // Update VX0
                    if (PORX[i * Nx1 + j] > 0) {
                        VXF0[i * Nx1 + j] = vxs[i * Nx1 + j] + vxD[i * Nx1 + j] / PORX[i * Nx1 + j];
                    }
                    // Update VY0
                    if (PORY[i * Nx1 + j] > 0) {
                        VYF0[i * Nx1 + j] = vys[i * Nx1 + j] + vyD[i * Nx1 + j] / PORY[i * Nx1 + j];
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
            PTF0_gko->copy_from(pt_gko.get());
            PTF0_gko->sub_scaled(sub_alpha, pf_gko);

            PT0 = pt;
            PF0 = pf;

            /////////////////////////////////////////////////////////////////////////////////////////
            // output
            /////////////////////////////////////////////////////////////////////////////////////////

            write_output(timestep);
        std::cout << "\n\n\n\n\n\n\n\n\n\n\nFinished one run through run_simulation successfully!" << std::endl;
        std::cout << "\n\n\n\n\n\n\n\n\n======================Push to git===================" << std::endl;
        }


}