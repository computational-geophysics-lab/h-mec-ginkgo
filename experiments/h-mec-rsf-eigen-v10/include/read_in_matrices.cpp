#include "read_in_matrices.h"
#include <eigen3/Eigen/Eigen>



void read_in_matrices(int &timestep){
    if (timestep > 0) {
        const string filename = nname + to_string(timestep) + ".h5";
        const string group_matrix = "Matrix";
        const string group_vector = "Vector";
        const string group_values = "Value";

        const hsize_t dims1[2] = {Ny, Nx};
        const string matrix_names[24] = {"SIIB", "OM0", "OM", "ARSF", "BRSF", "LRSF", "RHO", "ETA0", "ETA1", "ETA5", "ETA00", "IETAPLB", "SXY0", "YNY0", "KKK", "GGG", "COHC", "COHT", "FRIC", "FRIT", "DILC", "TTT", "EIIB", "VSLIPB"}; // {"names"} has to be the same as in *matrix
        int j = 0;
        for (auto i : {&SIIB, &OM0, &OM, &ARSF, &BRSF, &LRSF, &RHO, &ETA0, &ETA1, &ETA5, &ETA00, &IETAPLB, &SXY0, &YNY0, &KKK, &GGG, &COHC, &COHT, &FRIC, &FRIT, &DILC, &TTT, &EIIB, &VSLIPB}) { // {names} *matrix
            *i = read_matrix(filename, group_matrix, matrix_names[j], dims1);
            j++;
        }

        const hsize_t dims2[2] = {Ny1, Nx1};
        const string matrix_names_plus[30] = {"pt", "vxs", "vys", "pf", "vxD", "vyD", "ETAB", "ETAB0", "ETAP", "ETAP0", "POR", "GGGP", "GGGB", "PTF0", "PT0", "PF0", "SXX0", "SYY0", "RHOX", "RHOFX", "ETADX", "PORX", "VX0", "VXF0", "RHOY", "RHOFY", "ETADY", "PORY", "VY0", "VYF0"}; // {"names"} has to be the same as in *matrix_plus
        j = 0;
        for (auto i : {&pt, &vxs, &vys, &pf, &vxD, &vyD, &ETAB, &ETAB0, &ETAP, &ETAP0, &POR, &GGGP, &GGGB, &PTF0, &PT0, &PF0, &SXX0, &SYY0, &RHOX, &RHOFX, &ETADX, &PORX, &VX0, &VXF0, &RHOY, &RHOFY, &ETADY, &PORY, &VY0, &VYF0}) { // {names} *matrix_plus
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
        const Eigen::VectorXd temp = read_vector(filename, group_values, "values", dim2);
        timesum         = temp(0); 
        dt00            = temp(1); 
        dtx             = temp(2); 
        dty             = temp(3); 
        dtlapusta       = temp(4); 
        Vmax            = temp(5); 
        maxvxy          = temp(6); 
        dt              = temp(7), 
        yndtdecrease    = temp(8);

        if (antiplane) {
            const string group_antiplane = "Antiplane";
            SZX0    = read_matrix(filename, group_antiplane, "SZX0", dims1);
            SZY0    = read_matrix(filename, group_antiplane, "SZY0", dims1);
            VZ0     = read_matrix(filename, group_antiplane, "VZ0", dims2);

            maxvzsmod = read_vector(filename, group_antiplane, "maxvzsmod", dim1);
            minvzsmod = read_vector(filename, group_antiplane, "minvzsmod", dim1);
        }

        timestep++;

        int num_lines = 0;
        string line;
        ifstream myfile("./output_data/EVO_data.txt");
        while (getline(myfile, line)) {++num_lines;}
        myfile.close();
        
        if (num_lines >= timestep) {
            const string files[11] = {"./output_data/EVO_Vslip.txt", "./output_data/EVO_viscosity.txt", "./output_data/EVO_press_flu.txt", "./output_data/EVO_press_eff.txt", "./output_data/EVO_SigmaY.txt", "./output_data/EVO_Sii.txt", "./output_data/EVO_Theta.txt", "./output_data/EVO_Visc_comp.txt", "./output_data/EVO_Elast_comp.txt", "./output_data/EVO_VxD.txt", "./output_data/EVO_data.txt"};
            for (int i = 0; i < 11; i++) {
                remove_double_output(files[i], timestep - 1);
            }
        }
    } else {
        timestep = 1;
        // Basic nodes
        OM0     = Eigen::MatrixXd::Constant(Ny, Nx, omm(0));    // Old state parameter
        ARSF    = Eigen::MatrixXd::Constant(Ny, Nx, arsfm(1)); // a - parameter of RSF
        BRSF    = Eigen::MatrixXd::Constant(Ny, Nx, brsfm(1)); // b - parameter of RSF
        LRSF    = Eigen::MatrixXd::Constant(Ny, Nx, lrsfm(0)); // L - parameter of RSF
        
        // set matrices to 0
        for (auto i : {pt, vxs, vys, pf, vxD, vyD, RHO, ETA0, IETAPLB, SXY, SXY0, SZX0, SZY0, YNY0, KKK, GGG, COHC, COHT, FRIC, FRIT, DILC, TTT, EIIB, ETAB, ETAB0, ETAP, ETAP0, POR, GGGP, GGGB, PTF0, pt_ave, pf_ave, SXX, SXX0, SYY, SYY0, RHOX, RHOFX, ETADX, PORX, VX0, VXF0, RHOY, RHOFY, ETADY, PORY, VY0, VYF0, VSLIPB, VZ0}) {
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

        PT0 = Eigen::MatrixXd::Constant(Ny1, Nx1, PCONF + PTFDIFF);
        PF0 = Eigen::MatrixXd::Constant(Ny1, Nx1, PCONF);
    }
}