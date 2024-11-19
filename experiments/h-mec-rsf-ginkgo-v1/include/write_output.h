#ifndef WRITE_OUTPUT_H
#define WRITE_OUTPUT_H

#include "constants.hpp"
#include "global_variables.hpp"
#include "functions.hpp"
#include "hdf5.hpp"
#include "ginkgo_exec.hpp"

#include <iomanip>
#include <iostream>

using namespace std;

void write_output(int &timestep);

void write_row(std::string savefile_name, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& matrix, std::size_t row_index);
void write_vector(std::string savefile_name, const std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& vector);
void write_EVO(std::string filename, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& matrix, std::size_t row_index, double timesum, double dt);
void write_RSF_param(std::string filename, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& A_param, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& B_param, std::unique_ptr<gko::matrix::Dense<gko::default_precision>>& L_param, std::size_t row_index);

#endif
