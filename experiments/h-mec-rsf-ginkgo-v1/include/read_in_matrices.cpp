#include "read_in_matrices.h"
#include <ginkgo/ginkgo.hpp>

/*
remove " && false" from the if condition in line 9. This causes a segfault on gcc compiled machines
 */

void read_in_matrices(int &timestep){
    if (timestep > 0) {
        const string filename = nname + to_string(timestep) + ".h5";
        const string group_matrix = "Matrix";
        const string group_vector = "Vector";
        const string group_values = "Value";

        const hsize_t dims1[2] = {Ny, Nx};
        const string matrix_names[24] = {"SIIB", "OM0", "OM", "ARSF", "BRSF", "LRSF", "RHO", "ETA0", "ETA1", "ETA5", "ETA00", "IETAPLB", "SXY0", "YNY0", "KKK", "GGG", "COHC", "COHT", "FRIC", "FRIT", "DILC", "TTT", "EIIB", "VSLIPB"}; // {"names"} has to be the same as in *matrix
        int j = 0;



        for (auto i : {&SIIB_gko, &OM0_gko, &OM_gko, &ARSF_gko, &BRSF_gko, &LRSF_gko, &RHO_gko, &ETA0_gko, &ETA1_gko, &ETA5_gko, &ETA00_gko, &IETAPLB_gko, &SXY0_gko, &YNY0_gko, &KKK_gko, &GGG_gko, &COHC_gko, &COHT_gko, &FRIC_gko, &FRIT_gko, &DILC_gko, &TTT_gko, &EIIB_gko, &VSLIPB_gko}) { // {names} *matrix
            *i = read_matrix(filename, group_matrix, matrix_names[j], dims1);
            j++;
        }

        const hsize_t dims2[2] = {Ny1, Nx1};
        const string matrix_names_plus[30] = {"pt", "vxs", "vys", "pf", "vxD", "vyD", "ETAB", "ETAB0", "ETAP", "ETAP0", "POR", "GGGP", "GGGB", "PTF0", "PT0", "PF0", "SXX0", "SYY0", "RHOX", "RHOFX", "ETADX", "PORX", "VX0", "VXF0", "RHOY", "RHOFY", "ETADY", "PORY", "VY0", "VYF0"}; // {"names"} has to be the same as in *matrix_plus
        j = 0;
        for (auto i : {&pt_gko, &vxs_gko, &vys_gko, &pf_gko, &vxD_gko, &vyD_gko, &ETAB_gko, &ETAB0_gko, &ETAP_gko, &ETAP0_gko, &POR_gko, &GGGP_gko, &GGGB_gko, &PTF0_gko, &PT0_gko, &PF0_gko, &SXX0_gko, &SYY0_gko, &RHOX_gko, &RHOFX_gko, &ETADX_gko, &PORX_gko, &VX0_gko, &VXF0_gko, &RHOY_gko, &RHOFY_gko, &ETADY_gko, &PORY_gko, &VY0_gko, &VYF0_gko}) { // {names} *matrix_plus
            *i = read_matrix(filename, group_matrix, matrix_names_plus[j], dims2);
            j++;
        }
        
        const hsize_t dim1[1] = {num_timesteps};
        const string vector_names[6] = {"timesumcur", "dtcur", "maxvxsmod", "minvxsmod", "maxvysmod", "minvysmod"}; // {"names"} has to be the same as in *vec
        j = 0;
        for (auto i : {&timesumcur_gko, &dtcur_gko, &maxvxsmod_gko, &minvxsmod_gko, &maxvysmod_gko, &minvysmod_gko}) { // {names} *vec
            *i = read_vector(filename, group_vector, vector_names[j], dim1);
            j++;
        }

        const hsize_t dim3[1] = {marknum};
        const string vector_names_marker[5] = {"xm", "ym", "sxxm", "syym", "sxym"}; // {"names"} has to be the same as in *vec2
        j = 0;
        for (auto i : {&xm_gko, &ym_gko, &sxxm_gko, &syym_gko, &sxym_gko}) { // {names} *vec2
            *i = read_vector(filename, group_vector, vector_names_marker[j], dim3);
            j++;
        }

        const hsize_t dim2[1] = {9};
        auto temp_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(dim2[0], 1));
        double* temp = temp_gko->get_values();
        temp_gko = read_vector(filename, group_values, "values", dim2);
        timesum         = temp[0];
        dt00            = temp[1]; 
        dtx             = temp[2];
        dty             = temp[3];
        dtlapusta       = temp[4];
        Vmax            = temp[5];
        maxvxy          = temp[6];
        dt              = temp[7],
        yndtdecrease    = temp[8];

        if (antiplane) {
            const string group_antiplane = "Antiplane";
            SZX0_gko    = read_matrix(filename, group_antiplane, "SZX0", dims1);
            SZY0_gko    = read_matrix(filename, group_antiplane, "SZY0", dims1);
            VZ0_gko     = read_matrix(filename, group_antiplane, "VZ0", dims2);

            maxvzsmod_gko = read_vector(filename, group_antiplane, "maxvzsmod", dim1);
            minvzsmod_gko = read_vector(filename, group_antiplane, "minvzsmod", dim1);
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
        OM0_gko->fill(omm[0]);// Old state parameter
        ARSF_gko->fill(arsfm[1]); // a - parameter of RSF
        BRSF_gko->fill(brsfm[1]); // b - parameter of RSF
        LRSF_gko->fill(lrsfm[0]); // L - parameter of RSF


        // This list is stupid, however I don't know how to implement it as a loop
        pt_gko->fill(0);
        pt_gko->fill(0.0);
        vxs_gko->fill(0.0);
        vys_gko->fill(0.0);
        pf_gko->fill(0.0);
        vxD_gko->fill(0.0);
        vyD_gko->fill(0.0);
        RHO_gko->fill(0.0);
        ETA0_gko->fill(0.0);
        IETAPLB_gko->fill(0.0);
        SXY_gko->fill(0.0);
        SXY0_gko->fill(0.0);
        SZX0_gko->fill(0.0);
        SZY0_gko->fill(0.0);
        YNY0_gko->fill(0.0);
        KKK_gko->fill(0.0);
        GGG_gko->fill(0.0);
        COHC_gko->fill(0.0);
        COHT_gko->fill(0.0);
        FRIC_gko->fill(0.0);
        FRIT_gko->fill(0.0);
        DILC_gko->fill(0.0);
        TTT_gko->fill(0.0);
        EIIB_gko->fill(0.0);
        ETAB_gko->fill(0.0);
        ETAB0_gko->fill(0.0);
        ETAP_gko->fill(0.0);
        ETAP0_gko->fill(0.0);
        POR_gko->fill(0.0);
        GGGP_gko->fill(0.0);
        GGGB_gko->fill(0.0);
        PTF0_gko->fill(0.0);
        pt_ave_gko->fill(0.0);
        pf_ave_gko->fill(0.0);
        SXX_gko->fill(0.0);
        SXX0_gko->fill(0.0);
        SYY_gko->fill(0.0);
        SYY0_gko->fill(0.0);
        RHOX_gko->fill(0.0);
        RHOFX_gko->fill(0.0);
        ETADX_gko->fill(0.0);
        PORX_gko->fill(0.0);
        VX0_gko->fill(0.0);
        VXF0_gko->fill(0.0);
        RHOY_gko->fill(0.0);
        RHOFY_gko->fill(0.0);
        ETADY_gko->fill(0.0);
        PORY_gko->fill(0.0);
        VY0_gko->fill(0.0);
        VYF0_gko->fill(0.0);
        VSLIPB_gko->fill(0.0);
        VZ0_gko->fill(0.0);


        // Define Fault
        // #pragma omp parallel for (slower with n = 4)
        for (int i = 0; i < Ny; i++) {
            if (y[i] > upper_block && y[i] < lower_block) {
                for (int j = 0; j < Nx; j++) {
                    OM0[i * Nx + j] = omm[1];
                
                    if (x[j] < TS_1) {
                        BRSF[i * Nx + j] = brsfm[0];
                        ARSF[i * Nx + j] = arsfm[0];
                    }
                    if (x[j] >= TS_1 && x[j] < TS_2) {
                        BRSF[i * Nx + j] = brsfm[0] - (brsfm[0] - brsfm[1]) * ((x[j] - TS_1) / (TS_2 - TS_1));
                        ARSF[i * Nx + j] = arsfm[0] - (arsfm[0] - arsfm[1]) * ((x[j] - TS_1) / (TS_2 - TS_1));
                        LRSF[i * Nx + j] = lrsfm[0] - (lrsfm[0] - lrsfm[1]) * ((x[j] - TS_1) / (TS_2 - TS_1));
                    }
                    if (x[j] >= TS_2 && x[j] <= TS_3) {
                        LRSF[i * Nx + j] = lrsfm[1];
                    }
                    if (x[j] > TS_3 && x[j] <= TS_4) {
                        BRSF[i * Nx + j] = brsfm[1] - (brsfm[1] - brsfm[0]) * ((x[j] - TS_3) / (TS_4 - TS_3));
                        ARSF[i * Nx + j] = arsfm[1] - (arsfm[1] - arsfm[0]) * ((x[j] - TS_3) / (TS_4 - TS_3));
                        LRSF[i * Nx + j] = lrsfm[1] - (lrsfm[1] - lrsfm[0]) * ((x[j] - TS_3) / (TS_4 - TS_3));
                    }
                    if (x[j] > TS_4) {
                        BRSF[i * Nx + j] = brsfm[0];
                        ARSF[i * Nx + j] = arsfm[0];
                    }
                }
            }
        }
        OM = OM0;

        // set vectors to 0
        sxxm_gko->fill(0.0);
        syym_gko->fill(0.0);
        sxym_gko->fill(0.0);


        // The following ginkgo code replaces the two lines of Eigen::MatrixXd::Constant
        // PT0 = Eigen::MatrixXd::Constant(Ny1, Nx1, PCONF + PTFDIFF);
        // PF0 = Eigen::MatrixXd::Constant(Ny1, Nx1, PCONF);

        auto PT0_vec = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1));
        double* PT0 = PT0_vec->get_values();
        auto PF0_vec = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1));
        double* PF0 = PF0_vec->get_values();

        double PCOMBINED = PCONF + PTFDIFF; // only used for this for loop
        for (int i = 0; i<Ny1; i++){
            for (int j=0; j<Nx1; j++){
                PT0[i*Nx1 + j] = PCOMBINED;
                PF0[i*Nx1 + j] = PCONF;
            }
        }

        PT0_gko->fill(PCONF+PTFDIFF);
        PF0_gko->fill(PCONF);

    }
}