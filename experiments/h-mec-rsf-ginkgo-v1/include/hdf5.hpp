#ifndef HDF5_HPP
#define HDF5_HPP

#include <string>
#include <H5Cpp.h>
#include <ginkgo/ginkgo.hpp>
#include "ginkgo_exec.hpp"

using namespace std;
using namespace H5;

void create_file(const string &filename);
void add_group(const string &filename, const string &groupname);
void add_matrix(const string &filename, const string &groupname, const double* data, const string &dataset_name, const hsize_t* dims);
void add_vector(const string &filename, const string &groupname, const double* vec_data, const string &dataset_name, const hsize_t* dims);
std::unique_ptr<gko::matrix::Dense<double>> read_matrix(const string &filename, const string &groupname, const string &dataset_name, const hsize_t* dims);
std::unique_ptr<gko::matrix::Dense<double>> read_vector(const string &filename, const string &groupname, const string &dataset_name, const hsize_t* dims);

#endif