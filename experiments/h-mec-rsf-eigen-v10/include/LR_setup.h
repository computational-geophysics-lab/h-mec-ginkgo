#ifndef LR_SETUP_H
#define LR_SETUP_H

#include "constants.hpp"
#include "global_variables.hpp"
#include "functions.hpp"
#include "hdf5.hpp"

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/PardisoSupport>

using namespace std;

void LR_setup(Eigen::SparseMatrix<double> &L, double pfscale, double ptscale);

#endif