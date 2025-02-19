#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <H5Cpp.h>
#include <fstream>
#include <ginkgo/ginkgo.hpp>

// iostream and vector are only necessary for testing
#include <iostream>

#include "constants.hpp"
#include "ginkgo_exec.hpp"

using namespace std;



// lambda function that rounds towards 0
inline auto fix_towards_zero = [](const double temp) {
    return temp < 0. ? (int)ceil(temp) : (int)floor(temp);
};

// lambda function that checks if a value is between a min and max value
inline auto enforce_bounds = [](double val, const double min, const double max) {
    if (val < min) {
        val = min;
    } else if (val > max) {
        val = max;
    }
    return val;
};

// lambda function that squares each element of a 2x2 block. Needs the row size as input
inline auto square_block = [](double* matrix, int i, int j, int Nx) {
    return (pow(matrix[i * Nx + j], 2) + pow(matrix[i * Nx + j + 1], 2) + pow(matrix[(i + 1) * Nx + j], 2) + pow(matrix[(i + 1) * Nx + j + 1], 2));
};

// lambda function that computes: x / (x + y)
inline auto divplus = [](const double x, const double y) {
    return (x / (x + y));
};

double safe_difference(double a,  double b,  double significant_digits);
double randd();
double findmin(double* matrix, int Nx, int Ny);
double findmax(double* matrix, int Nx, int Ny);
int check_bounds(int k, int bound);
void copy_bounds(std::unique_ptr<gko::matrix::Dense<double>>& matrix, int Nx, int Ny);
void remove_double_output( string file_name, int num);

void add_block(double* matrix, double block_matrix[4], double factor, int i, int j, int Nx, int Ny);

// Computes and returns the average of a (block_size_y x block_size_x) block starting at position i, j
double BlockSum(double* matrix, int i, int j, int block_size_y, int block_size_x, int Ny, int Nx);


// Computes the minimal entry within a block
double BlockMin(double* matrix, int i, int j, int block_size_y, int block_size_x, int Ny, int Nx);

// Print a gko matrix
void print_matrix(double* matrix, int Nx, int Ny);
void save_matrix(const std::string filename, double* matrix, int Nx, int Ny);

// Only in testing, take out for final version:
std::vector<int> read_ind_vector(const std::string &filepath);
std::vector<double> read_L_supposed(const std::string &filepath);


#endif