#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <eigen3/Eigen/Eigen>
#include <H5Cpp.h>
#include "constants.hpp"
#include <fstream>

using namespace std;


// lambda function that rounds towards 0
auto fix_towards_zero = [](const double temp) {
    return temp < 0. ? (int)ceil(temp) : (int)floor(temp);
};

// lambda function that checks if a value is between a min and max value
auto enforce_bounds = [](double val, const double min, const double max) {
    if (val < min) {
        val = min;
    } else if (val > max) {
        val = max;
    }
    return val;
};

// lambda function that squares each element of a 2x2 block
auto square_block = [](const Eigen::Matrix2d mat) {
    return (pow(mat(0, 0), 2) + pow(mat(0, 1), 2) + pow(mat(1, 0), 2) + pow(mat(1, 1), 2));
};

// lambda function that computes: x / (x + y)
auto divplus = [](const double x, const double y) {
    return (x / (x + y));
};

int check_bounds(int k, const int bound);
void copy_bounds(Eigen::MatrixXd& Temp);
void remove_double_output(const string file_name, const int num);

#endif
