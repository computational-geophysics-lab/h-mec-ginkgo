#include <string>
#include <eigen3/Eigen/Eigen>
#include <H5Cpp.h>
#include "constants.hpp"
#include <fstream>

using namespace std;


// function that checks if a integer value is bigger than 0 and lower than a set bound
int check_bounds(int k, const int bound) {
    if (k < 0) {
        k = 0;
    } else if (k > bound - 2) {
        k = bound - 2;
    }
    return k;
}

// function that sets the boundary values to values of one column/row further in
void copy_bounds(Eigen::MatrixXd& Temp) {
    Temp.col(0) = Temp.col(1);
    Temp.col(Nx) = Temp.col(Nx - 1);
    Temp.row(0) = Temp.row(1);
    Temp.row(Ny) = Temp.row(Ny - 1);
}

// function that removes already calculated results in some output files so that these are in order (of timesteps)
void remove_double_output(const string file_name, const int num) {
    ifstream file(file_name);
    string new_content = "";
    string line;
    for (int i = 0; i < num; i++) {
        getline(file, line);
        new_content += line + '\n';
    }
    file.close();
    ofstream file_out(file_name);
    file_out << new_content;
    file_out.close();
}