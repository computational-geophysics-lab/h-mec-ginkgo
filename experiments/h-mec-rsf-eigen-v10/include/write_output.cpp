#include "write_output.h"

using namespace std;
using namespace H5;

void write_output(int &timestep){
    cout << "====================================" << endl;
    cout << "total time:        " << timesum << " sec" << endl;
    cout << "time step:         " << dt << " sec" << endl;
    cout << "Vslip max:         " << Vmax << endl;
    cout << "iter - iterations:   " << iterstep + 1 << endl;
    cout << "global - iterations: " << ynlast + 1 << endl;
    
    if (timestep == 1) {
        ofstream out_fault;
        out_fault.open("./output_data/x_fault.txt");
        out_fault << x << endl;
        out_fault.close();

        ofstream out_rsf;
        out_rsf.open("./output_data/rsf_fault.txt", ios_base::app | ios_base::out);
        out_rsf << ARSF.row(line_fault) << "\n\n" << BRSF.row(line_fault) << "\n\n" << LRSF.row(line_fault) << endl;
        out_rsf.close();
    }

    // ========== save slip rate
    ofstream out_Vslip("./output_data/EVO_Vslip.txt", ios_base::app | ios_base::out);
    out_Vslip << timesum << "    " << dt << "    " << VSLIPB.row(line_fault) << endl;
    out_Vslip.close();
    
    // ========== save viscosity
    ofstream out_viscosity("./output_data/EVO_viscosity.txt", ios_base::app | ios_base::out);
    out_viscosity << timesum << "    " << dt << "    " << ETA.row(line_fault) << endl;
    out_viscosity.close();
    
    // ========== save fluid 
    ofstream out_press_flu("./output_data/EVO_press_flu.txt", ios_base::app | ios_base::out);
    out_press_flu << timesum << "    " << dt << "    " << pf.row(line_fault) << endl;
    out_press_flu.close();
    
    // ========== save effective pressure
    ofstream out_press_eff("./output_data/EVO_press_eff.txt", ios_base::app | ios_base::out);
    Eigen::MatrixXd P_diff = pt - pf;
    out_press_eff << timesum << "    " << dt << "    " << P_diff.row(line_fault) << endl;
    out_press_eff.close();
    
    // ========== save SigmaY
    ofstream out_SigmaY("./output_data/EVO_SigmaY.txt", ios_base::app | ios_base::out);
    out_SigmaY << timesum << "    " << dt << "    " << SigmaY.row(line_fault) << endl;
    out_SigmaY.close();
    
    // ========== save SII
    ofstream out_Sii("./output_data/EVO_Sii.txt", ios_base::app | ios_base::out);
    out_Sii << timesum << "    " << dt << "    " << SII_fault.row(line_fault) << endl;
    out_Sii.close();
    
    // ========== save Theta
    ofstream out_Theta("./output_data/EVO_Theta.txt", ios_base::app | ios_base::out);
    out_Theta << timesum << "    " << dt << "    " << OM.row(line_fault) << endl;
    out_Theta.close();
    
    // ========== save viscous compaction
    ofstream out_Visc_comp("./output_data/EVO_Visc_comp.txt", ios_base::app | ios_base::out);
    out_Visc_comp << timesum << "    " << dt << "    " << VIS_COMP.row(line_fault) << endl;
    out_Visc_comp.close();
    
    // ========== save elastic compaction
    ofstream out_Elast_comp("./output_data/EVO_Elast_comp.txt", ios_base::app | ios_base::out);
    out_Elast_comp << timesum << "    " << dt << "    " << EL_DECOM.row(line_fault) << endl;
    out_Elast_comp.close();

    // ========== save vx Darcy
    ofstream out_EVO_vxD("./output_data/EVO_vxD.txt", ios_base::app | ios_base::out);
    out_EVO_vxD << timesum << "    " << dt << "    " << vxD.row(line_fault) << endl;
    out_EVO_vxD.close();
    
    // ========== save time, dt, vmax
    ofstream out_data("./output_data/EVO_data.txt", ios_base::app | ios_base::out);
    out_data << setw(20) << timesum << setw(20) << dt << setw(20) << Vmax << setw(20) << ynlast + 1 << setw(20) << iterstep + 1 << endl;
    out_data.close();
    
    if (timestep % savestep == 0) {
        ofstream out_file;
        out_file.open("./input_data/StartingTimestep.txt");
        out_file << timestep << endl;
        out_file.close();

        const string save_file_name = "./output_data/" + nname + to_string(timestep) + ".h5";
        const string group_matrix = "Matrix";
        const string group_vector = "Vector";
        const string group_values = "Value";

        create_file(save_file_name);
        for (auto i : {group_matrix,  group_vector, group_values}){
            add_group(save_file_name, i);
        }

        // save matrices and vectors to restart the simulation from a certain timestep.
        // !!! Always change both the string array and loop order !!!
        const hsize_t dims1[2] = {Ny, Nx};
        const string matrix_names[24] = {"SIIB", "OM0", "OM", "ARSF", "BRSF", "LRSF", "RHO", "ETA0", "ETA1", "ETA5", "ETA00", "IETAPLB", "SXY0", "YNY0", "KKK",
                                    "GGG", "COHC", "COHT", "FRIC", "FRIT", "DILC", "TTT", "EIIB", "VSLIPB"}; // {"names"} has to be the same as in *matrix
        int j = 0;
        for (auto i : {SIIB, OM0, OM, ARSF, BRSF, LRSF, RHO, ETA0, ETA1, ETA5, ETA00, IETAPLB, SXY0, YNY0, KKK, GGG, COHC, COHT, FRIC, FRIT, DILC, TTT, EIIB, VSLIPB}) { // {names} *matrix
            add_matrix(save_file_name, group_matrix, i, matrix_names[j], dims1);
            j++;
        }

        const hsize_t dims2[2] = {Ny1, Nx1};
        const string matrix_names_plus[30] = {"pt", "vxs", "vys", "pf", "vxD", "vyD", "ETAB", "ETAB0", "ETAP", "ETAP0", "POR", "GGGP", "GGGB", "PTF0", "PT0", "PF0",
                                        "SXX0", "SYY0", "RHOX", "RHOFX", "ETADX", "PORX", "VX0", "VXF0", "RHOY", "RHOFY", "ETADY", "PORY", "VY0", "VYF0"}; // {"names"} has to be the same as in *matrix_plus
        j = 0;
        for (auto i : {pt, vxs, vys, pf, vxD, vyD, ETAB, ETAB0, ETAP, ETAP0, POR, GGGP, GGGB, PTF0, PT0, PF0, SXX0, SYY0, RHOX, RHOFX, ETADX,
                        PORX, VX0, VXF0, RHOY, RHOFY, ETADY, PORY, VY0, VYF0}) { // {names} *matrix_plus
            add_matrix(save_file_name, group_matrix, i, matrix_names_plus[j], dims2);
            j++;
        }
        
        const hsize_t dim1[1] = {num_timesteps};
        const string vector_names[6] = {"timesumcur", "dtcur", "maxvxsmod", "minvxsmod", "maxvysmod", "minvysmod"}; // {"names"} has to be the same as in *vec
        j = 0;
        for (auto i : {timesumcur, dtcur, maxvxsmod, minvxsmod, maxvysmod, minvysmod}) { // {names} *vec
            add_vector(save_file_name, group_vector, i, vector_names[j], dim1);
            j++;
        }

        const hsize_t dim2[1] = {marknum};
        const string vector_names_marker[5] = {"xm", "ym", "sxxm", "syym", "sxym"}; // {"names"} has to be the same as in *vec2
        j = 0;
        for (auto i : {xm, ym, sxxm, syym, sxym}) { // {names} *vec2
            add_vector(save_file_name, group_vector, i, vector_names_marker[j], dim2);
            j++;
        }

        const hsize_t dim3[1] = {9};
        Eigen::VectorXd temp(9);
        temp << timesum, dt00, dtx, dty, dtlapusta, Vmax, maxvxy, dt, yndtdecrease;
        add_vector(save_file_name, group_values, temp, "values", dim3);

        if (antiplane) {
            const string group_antiplane = "Antiplane";
            add_group(save_file_name, group_antiplane);
            add_matrix(save_file_name, group_antiplane, SZX0, "SZX0", dims1);
            add_matrix(save_file_name, group_antiplane, SZY0, "SZY0", dims1);

            add_matrix(save_file_name, group_antiplane, VZ0, "VZ0", dims2);

            add_vector(save_file_name, group_antiplane, maxvzsmod, "maxvzsmod", dim1);
            add_vector(save_file_name, group_antiplane, minvzsmod, "minvzsmod", dim1);
        }
    }
}