#include <string>
#include <ginkgo/ginkgo.hpp>
#include <H5Cpp.h>
#include <fstream>
#include "functions.hpp"

#include <iostream>
#include <vector>

using namespace std;
/* --------------Ginkgo-------------------Ginkgo----------------------
Changes:
 - Added function GinkgoLinSpaced to create a ginkgo vector with linear spacing, analog to Eigen::LinSpaced -> moved to constants.hpp, as I had troubles overloading the function because constants.hpp depended on functions.hpp, and functions.hpp depended on constants.hpp
 - the function copy_bounds takes Eigen::MatrixXd& Temp as argument -> change to ginkgo
Todo:
*/

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
void copy_bounds(std::unique_ptr<gko::matrix::Dense<double>>& matrix, int Nx, int Ny) {
    double* data = matrix->get_values(); // Access the underlying data
    auto leading_dim = matrix->get_stride(); // Stride (distance between rows)
    if(leading_dim !=Nx) {
        std::cout << "\nNy = " << Nx << " and the stride = " << leading_dim << "!\n";
        throw std::out_of_range("Ny of the matrix and its stride are not identical!");
    }
    // Set first and last columns based on neighboring columns
    for (int row = 0; row < Ny; ++row) {
        data[row * leading_dim + 0] = data[row * leading_dim + 1]; // First column
        data[row * leading_dim + (Nx - 1)] = data[row * leading_dim + (Nx - 2)]; // Last column
    }

    // Set first and last rows based on neighboring rows
    for (int col = 0; col < Nx; ++col) {
        data[0 * leading_dim + col] = data[1 * leading_dim + col]; // First row
        data[(Ny - 1) * leading_dim + col] = data[(Ny - 2) * leading_dim + col]; // Last row
    }
}
/*
// function that sets the boundary values to values of one column/row further in
void copy_bounds(Eigen::MatrixXd& Temp) {
    Temp.col(0) = Temp.col(1);
    Temp.col(Nx) = Temp.col(Nx - 1);
    Temp.row(0) = Temp.row(1);
    Temp.row(Ny) = Temp.row(Ny - 1);
}
*/

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


// Creates a random double between 0 and 1. If truly random numbers are needed, this needs to be improved
double randd() {
  return (double)rand() / (RAND_MAX + 1.0);
}

// Calculates a-b, then returns the difference, except if a and b are so close to each other that the difference is smaller by the reduction factor, then it returns 0 instead. This is useful if a=b=10^19 except for the last bit, which would return 1024 but is insignificant
double safe_difference(const double a, const double b, const double significant_digits) {
    const double diff = a - b;
    if (abs(diff) < significant_digits * abs(a)) {
        return 0.0;
    } else {
        return diff;
    }
}



void add_block(double* matrix, double block_matrix[4], double factor, int i, int j, int Nx, int Ny){
    matrix[i*Nx + j] += factor*block_matrix[0];
            matrix[i*Nx + j+1] += factor*block_matrix[1];
            matrix[(i+1)*Nx + j] += factor*block_matrix[2];
            matrix[(i+1)*Nx + j+1] += factor*block_matrix[3];
}

// Computes and returns the average of a (block_size_y x block_size_x) block starting at position i, j
double BlockSum(double* matrix, int i, int j, int block_size_y, int block_size_x, int Ny, int Nx){
    if (i+block_size_y > Ny) {
        std::cout << "The BlockSum is out of range: " << i+block_size_y << " is larger than " << Ny;
        throw std::out_of_range("Col index is out of bounds in BlockSum");
    }
    if (j+block_size_x > Nx) {
        std::cout << "The BlockSum is out of range: " << j+block_size_x << " is larger than " << Nx;
        throw std::out_of_range("Col index is out of bounds in BlockSum");
    }

    double sum = 0.0;
    for (int row=i; row < i+block_size_y;row++){
        for (int col=j; col < j + block_size_x; col++) {
            sum += matrix[row*Nx + col];
        }
    }
    return sum;
}

// Computes the minimal entry within a block
double BlockMin(double* matrix, int i, int j, int block_size_y, int block_size_x, int N_y, int N_x){
    double blockmin = matrix[i*N_x+j];

    for (int row=i; row<i+block_size_y; row++){
        for (int col=j; col<j+block_size_x; col++){
            if (matrix[row*N_x + col]<blockmin){
                blockmin = matrix[row*N_x + col];
            }
        }
    }
    return blockmin;
}

double findmax(double* matrix, int Nx, int Ny) {
    double max = 0.0;
    for (int l=0;l<Nx*Ny;l++) {
        if (matrix[l]>max) {
            max = matrix[l];
        }
    }
    return max;
}


double findmin(double* matrix, int Nx, int Ny) {
    double min = matrix[0];
    for (int l=0;l<Nx*Ny;l++) {
        if (matrix[l]<min) {
            min = matrix[l];
        }
    }
    return min;
}

void print_matrix(double* matrix, int Nx, int Ny) {
    std::cout << std::endl;
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            std::cout << matrix[i*Nx + j] << " ";
        }
        std::cout << "\n";
    }
}
void save_matrix(const std::string filename, double* matrix, int Nx, int Ny) {
    std::ofstream mat_output;
    mat_output.open(filename, std::ios_base::app |std::ios_base::out);

    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            mat_output << matrix[i*Nx + j] << " ";
        }
        mat_output << "\n";
    }
    mat_output.close();
}

// This is a improvised function that should not be part of the final version. It is used to set elements of the L_matrix
std::vector<int> read_ind_vector(const std::string &filepath) {
    // Open the file for reading
    ifstream fin(filepath);

    // Create an empty vector
    vector<char> v_c;

    // Read the contents of the file and store them in the
    // vector
    char c;
    while (fin >> c) {
        v_c.push_back(c);
    }

    // Close the file
    fin.close();

    if(v_c.size()!=269*6) {
        std::cerr << "The list" << filepath << " does not contain 269 elements!";
    }

    vector<int> indices;
    for (int i = 0; i < 269; i++) {
        int ind = v_c[i*6]*100000 + v_c[i*6+1]*10000 + v_c[i*6+2]*1000 + v_c[i*6+3]*100 + v_c[i*6+4]*10 + v_c[i*6+5] - '0' * 111111;
        indices.push_back(ind);
    }

    return indices;
}

std::vector<double> read_L_supposed(const std::string &filepath) {
        // make size at least as large as needed
        const int size = 130000;
        //int array[size];
        std::vector<double> array;

        ifstream file(filepath);

        int count = 0;
        double x;

        // check that array is not already full
        // and read integer from file,

        while (count < size && file >> x) {
            //array[count++] = x;
            array.push_back(x);
            count++;
        }


        // display the values stored in the array
        /*for (int i = 0; i < count; i++) {
            cout << array[i] <<' ';
        }*/
        std::cout << "\nThe count is = " << count;
        return array;
}