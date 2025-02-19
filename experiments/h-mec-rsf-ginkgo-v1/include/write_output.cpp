#include "write_output.h"

using namespace std;
using namespace H5;



// prints row of gko matrix
void write_row(std::string savefile_name, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& matrix, std::size_t row_index) {
    if (row_index >= matrix->get_size()[0]) {
        throw std::out_of_range("Row index is out of bounds");
    }

    // Get the data pointer to the beginning of the desired row
    auto row_data = matrix->get_const_values() + row_index * matrix->get_stride();

    // Get the number of columns
    auto num_cols = matrix->get_size()[1];

    // Print the row values
    std::ofstream output;
    output.open(savefile_name, std::ios_base::app | std::ios_base::out);
    for (std::size_t i = 0; i < num_cols; ++i) {
        output << row_data[i] << " ";
    }
    output << std::endl;
    output.close();
}

void write_vector(std::string savefile_name, const std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& vector) {

    // Get the data pointer to the beginning of the desired row
    auto vec_data = vector->get_const_values();

    // Get the number of columns
    auto length_vec = vector->get_size()[0];

    // Print the row values
    std::ofstream output;
    output.open(savefile_name, std::ios_base::app | std::ios_base::out);
    for (std::size_t i = 0; i < length_vec; ++i) {
        output << vec_data[i] << " ";
    }
    output << std::endl;
    output.close();
}

void write_EVO(std::string filename, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& matrix, std::size_t row_index, double timesum, double dt) {
    if (row_index >= matrix->get_size()[0]) {
        throw std::out_of_range("Row index is out of bounds");
    }

    // Get the data pointer to the beginning of the desired row
    auto row_data = matrix->get_const_values() + row_index * matrix->get_stride();

    // Get the number of columns
    auto num_cols = matrix->get_size()[1];

    // Print the row values
    std::ofstream output;
    output.open(filename, std::ios_base::app | std::ios_base::out);


    output << timesum << "    " << dt << "    ";
    for (std::size_t i = 0; i < num_cols; ++i) {
        output << row_data[i] << " ";
    }
    output << std::endl;
    output.close();
}

void write_RSF_param(std::string filename, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& A_param, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& B_param, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& L_param, std::size_t row_index) {
    if (row_index >= A_param->get_size()[0]) {
        throw std::out_of_range("Row index is out of bounds");
    }

    // Get the data pointer to the beginning of the desired row
    auto A_row_data = A_param->get_const_values() + row_index * A_param->get_stride();
    auto B_row_data = B_param->get_const_values() + row_index * A_param->get_stride();
    auto L_row_data = L_param->get_const_values() + row_index * A_param->get_stride();

    // Get the number of columns
    auto num_cols = A_param->get_size()[1];

    // Print the row values
    std::ofstream output;
    output.open(filename, std::ios_base::app | std::ios_base::out);

    for (std::size_t i = 0; i < num_cols; ++i) {
        output << A_row_data[i] << " ";
    }
    output << "\n\n";
    for (std::size_t i = 0; i < num_cols; ++i) {
        output << B_row_data[i] << " ";
    }
    output << "\n\n";
    for (std::size_t i = 0; i < num_cols; ++i) {
        output << L_row_data[i] << " ";
    }
    output << std::endl;
    output.close();
}

void write_output(int &timestep){
    // To subtract gko matrices I only found the sub_scaled function, that also takes a factor alpha, where alpha has to be a gko::matrix. I just set it to 1.0 I assume there is a better solution, but I didn't find it
    auto sub_alpha = gko::matrix::Dense<double>::create(exec, gko::dim<2>(1,1)); // Need a alpha parameter when subtracting matrices, set to 1.0
    sub_alpha->at(0,0)=1.0; // Need a alpha parameter when subtracting matrices -> set to 1.0

    cout << "====================================" << endl;
    cout << "total time:        " << timesum << " sec" << endl;
    cout << "time step:         " << dt << " sec" << endl;
    cout << "Vslip max:         " << Vmax << endl;
    cout << "iter - iterations:   " << iterstep + 1 << endl;
    cout << "global - iterations: " << ynlast + 1 << endl;
    
    if (timestep == 1) {

        write_vector("./output_data/x_fault.txt", x_gko);
        /*
        ofstream out_fault;
        out_fault.open("./output_data/x_fault.txt");
        out_fault << x << endl;
        out_fault.close();
        */

       //write_row(".output_data/rsf_fault.txt", ARSF_gko, line_fault);
       write_RSF_param("./output_data/rsf_fault.txt", ARSF_gko, BRSF_gko, LRSF_gko, line_fault);
        /*ofstream out_rsf;
        out_rsf.open("./output_data/rsf_fault.txt", ios_base::app | ios_base::out);
        out_rsf << ARSF.row(line_fault) << "\n\n" << BRSF.row(line_fault) << "\n\n" << LRSF.row(line_fault) << endl;
        out_rsf.close();*/
        
    }



    // ========== save slip rate
    write_EVO("./output_data/EVO_Vslip.txt", VSLIPB_gko, line_fault, timesum, dt);
    //ofstream out_Vslip("./output_data/EVO_Vslip.txt", ios_base::app | ios_base::out);
    //out_Vslip << timesum << "    " << dt << "    " << VSLIPB.row(line_fault) << endl;
    //out_Vslip.close();
    
    // ========== save viscosity
    write_EVO("./output_data/EVO_viscosity.txt", ETA_gko, line_fault, timesum, dt);
    //ofstream out_viscosity("./output_data/EVO_viscosity.txt", ios_base::app | ios_base::out);
    //out_viscosity << timesum << "    " << dt << "    " << ETA.row(line_fault) << endl;
    //out_viscosity.close();
    
    // ========== save fluid
    write_EVO("./output_data/EVO_press_flu.txt", pf_gko, line_fault, timesum, dt);
    //ofstream out_press_flu("./output_data/EVO_press_flu.txt", ios_base::app | ios_base::out);
    //out_press_flu << timesum << "    " << dt << "    " << pf.row(line_fault) << endl;
    //out_press_flu.close();
    
    // ========== save effective pressure
    //ofstream out_press_eff("./output_data/EVO_press_eff.txt", ios_base::app | ios_base::out);
    auto P_diff_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(Ny1, Nx1));
    P_diff_gko->copy_from(pt_gko.get());
    P_diff_gko->sub_scaled(sub_alpha, pf_gko);
    //out_press_eff << timesum << "    " << dt << "    " << P_diff_gko.row(line_fault) << endl;
    //out_press_eff.close();
    write_EVO("./output_data/EVO_press_eff.txt", P_diff_gko, line_fault, timesum, dt);

    // ========== save SigmaY
    write_EVO("./output_data/EVO_SigmaY.txt", SigmaY_gko, line_fault, timesum, dt);
    //ofstream out_SigmaY("./output_data/EVO_SigmaY.txt", ios_base::app | ios_base::out);
    //out_SigmaY << timesum << "    " << dt << "    " << SigmaY.row(line_fault) << endl;
    //out_SigmaY.close();
    
    // ========== save SII
    write_EVO("./output_data/EVO_Sii.txt", SII_fault_gko, line_fault, timesum, dt);
    //ofstream out_Sii("./output_data/EVO_Sii.txt", ios_base::app | ios_base::out);
    //out_Sii << timesum << "    " << dt << "    " << SII_fault.row(line_fault) << endl;
    //out_Sii.close();
    
    // ========== save Theta
    write_EVO("./output_data/EVO_Theta.txt", OM_gko, line_fault, timesum, dt);
    //ofstream out_Theta("./output_data/EVO_Theta.txt", ios_base::app | ios_base::out);
    //out_Theta << timesum << "    " << dt << "    " << OM.row(line_fault) << endl;
    //out_Theta.close();
    
    // ========== save viscous compaction
    write_EVO("./output_data/EVO_Visc.txt", VIS_COMP_gko, line_fault, timesum, dt);
    //ofstream out_Visc_comp("./output_data/EVO_Visc_comp.txt", ios_base::app | ios_base::out);
    //out_Visc_comp << timesum << "    " << dt << "    " << VIS_COMP.row(line_fault) << endl;
    //out_Visc_comp.close();
    
    // ========== save elastic compaction
    write_EVO("./output_data/EVO_Elast_comp.txt", EL_DECOM_gko, line_fault, timesum, dt);
    //ofstream out_Elast_comp("./output_data/EVO_Elast_comp.txt", ios_base::app | ios_base::out);
    //out_Elast_comp << timesum << "    " << dt << "    " << EL_DECOM.row(line_fault) << endl;
    //out_Elast_comp.close();

    // ========== save vx Darcy
    write_EVO("./output_data/EVO_vxD.txt", vxD_gko, line_fault, timesum, dt);
    //ofstream out_EVO_vxD("./output_data/EVO_vxD.txt", ios_base::app | ios_base::out);
    //out_EVO_vxD << timesum << "    " << dt << "    " << vxD.row(line_fault) << endl;
    //out_EVO_vxD.close();
    
    // ========== save time, dt, vmax
    ofstream out_data("./output_data/EVO_data.txt", ios_base::app | ios_base::out);
    out_data << setw(20) << timesum << setw(20) << dt << setw(20) << Vmax << setw(20) << ynlast + 1 << setw(20) << iterstep + 1 << endl;
    out_data.close();
    
    if (timestep % savestep == 0) {
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

        add_matrix(save_file_name, group_matrix, SIIB, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, OM0, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, OM, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, ARSF, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, BRSF, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, LRSF, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, RHO, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, ETA0, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, ETA1, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, ETA5, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, ETA00, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, IETAPLB, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, SXY0, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, YNY0, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, KKK, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, GGG, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, COHC, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, COHT, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, FRIC, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, FRIT, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, DILC, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, TTT, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, EIIB, matrix_names[j], dims1);
        j++;
        add_matrix(save_file_name, group_matrix, VSLIPB, matrix_names[j], dims1);



        const hsize_t dims2[2] = {Ny1, Nx1};
        const string matrix_names_plus[30] = {"pt", "vxs", "vys", "pf", "vxD", "vyD", "ETAB", "ETAB0", "ETAP", "ETAP0", "POR", "GGGP", "GGGB", "PTF0", "PT0", "PF0",
                                        "SXX0", "SYY0", "RHOX", "RHOFX", "ETADX", "PORX", "VX0", "VXF0", "RHOY", "RHOFY", "ETADY", "PORY", "VY0", "VYF0"}; // {"names"} has to be the same as in *matrix_plus
        j = 0;
        add_matrix(save_file_name, group_matrix,pt, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,vxs, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,vys, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,pf, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,vxD, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,vyD, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,ETAB, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,ETAB0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,ETAP, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,ETAP0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,POR, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,GGGP, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,GGGB, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,PTF0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,PT0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,PF0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,SXX0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,SYY0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,RHOX, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,RHOFX, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,ETADX, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,PORX, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,VX0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,VXF0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,RHOY, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,RHOFY, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,ETADY, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,PORY, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,VY0, matrix_names_plus[j], dims2);
        j++;
        add_matrix(save_file_name, group_matrix,VYF0, matrix_names_plus[j], dims2);

        
        const hsize_t dim1[1] = {num_timesteps};
        const string vector_names[6] = {"timesumcur", "dtcur", "maxvxsmod", "minvxsmod", "maxvysmod", "minvysmod"}; // {"names"} has to be the same as in *vec
        j = 0;
        add_vector(save_file_name, group_vector, timesumcur, vector_names[j], dim1);
        j++;
        add_vector(save_file_name, group_vector, dtcur, vector_names[j], dim1);
        j++;
        add_vector(save_file_name, group_vector, maxvxsmod, vector_names[j], dim1);
        j++;
        add_vector(save_file_name, group_vector, minvxsmod, vector_names[j], dim1);
        j++;
        add_vector(save_file_name, group_vector, maxvysmod, vector_names[j], dim1);
        j++;
        add_vector(save_file_name, group_vector, minvysmod, vector_names[j], dim1);


        const hsize_t dim2[1] = {marknum};
        const string vector_names_marker[5] = {"xm", "ym", "sxxm", "syym", "sxym"}; // {"names"} has to be the same as in *vec2
        j = 0;
        add_vector(save_file_name, group_vector, xm, vector_names_marker[j], dim2);
        j++;
        add_vector(save_file_name, group_vector, ym, vector_names_marker[j], dim2);
        j++;
        add_vector(save_file_name, group_vector, sxxm, vector_names_marker[j], dim2);
        j++;
        add_vector(save_file_name, group_vector, syym, vector_names_marker[j], dim2);
        j++;
        add_vector(save_file_name, group_vector, sxym, vector_names_marker[j], dim2);

        const hsize_t dim3[1] = {9};
        auto temp_gko = gko::matrix::Dense<double>::create(exec, gko::dim<2>(9, 1));
        double* temp = temp_gko->get_values();
        temp[0] = timesum;
        temp[1] = dt00;
        temp[2] = dtx;
        temp[3] = dty;
        temp[4] = dtlapusta;
        temp[5] = Vmax;
        temp[6] = maxvxy;
        temp[7] = dt;
        temp[8] = yndtdecrease;

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

        ofstream out_file;
        out_file.open("./input_data/StartingTimestep.txt");
        out_file << timestep << endl;
        out_file.close();
    }
}