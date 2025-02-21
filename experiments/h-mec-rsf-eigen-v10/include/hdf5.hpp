#ifndef HDF5_H
#define HDF5_H

#include <string>
#include <eigen3/Eigen/Eigen>
#include <H5Cpp.h>

using namespace std;
using namespace H5;

void create_file(const string &filename);
void add_group(const string &filename, const string &groupname);
void add_matrix(const string &filename, const string &groupname, const Eigen::MatrixXd& data, const string &dataset_name, const hsize_t* dims);
void add_vector(const string &filename, const string &groupname, const Eigen::VectorXd& data, const string &dataset_name, const hsize_t* dims);
Eigen::MatrixXd read_matrix(const string &filename, const string &groupname, const string &dataset_name, const hsize_t* dims);
Eigen::VectorXd read_vector(const string &filename, const string &groupname, const string &dataset_name, const hsize_t* dims);

#endif