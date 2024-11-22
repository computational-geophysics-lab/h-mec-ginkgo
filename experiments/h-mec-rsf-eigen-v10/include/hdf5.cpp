#include <string>
#include <eigen3/Eigen/Eigen>
#include <H5Cpp.h>

using namespace std;
using namespace H5;

/*
inputs:
    string filename
*/
void create_file(const string &filename) {
    H5File *file = new H5File(filename, H5F_ACC_TRUNC);

    delete file;
}

/*
inputs:
    string filename
    string groupname
*/
void add_group(const string &filename, const string &groupname) {
    H5File *file = new H5File(filename, H5F_ACC_RDWR);
    Group* group = new Group(file->createGroup("/" + groupname));

    delete group;
    delete file;
}

/*
inputs:
    string filename
    string groupname
    MatrixXd data -> pointer to matrix to store
    string dataset_name -> name the stored data
    hsize_t* dims -> 2 dimensional array with dimensions (x and y) of stored data
*/
void add_matrix(const string &filename, const string &groupname, const Eigen::MatrixXd& data, const string &dataset_name, const hsize_t* dims) {
    H5File *file = new H5File(filename, H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup(groupname));
    DataSpace *dataspace = new DataSpace(2, dims);
    DataSet *dataset = new DataSet(file->createDataSet("/" + groupname + "/" + dataset_name, PredType::NATIVE_DOUBLE, *dataspace));

    double out[dims[0]][dims[1]];
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            out[i][j] = data(i, j);
        }
    }

    dataset->write(out, PredType::NATIVE_DOUBLE);

    delete dataset;
    delete dataspace;
    delete group;
    delete file;
}

/*
inputs:
    string filename
    string groupname
    VectorXd data -> pointer to vector to store
    string dataset_name -> name the stored data
    hsize_t* dims -> 1 dimensional array with dimension (x) of stored data
*/
void add_vector(const string &filename, const string &groupname, const Eigen::VectorXd& data, const string &dataset_name, const hsize_t* dims) {
    H5File *file = new H5File(filename, H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup(groupname));
    DataSpace *dataspace = new DataSpace(1, dims);
    DataSet *dataset = new DataSet(file->createDataSet("/" + groupname + "/" + dataset_name, PredType::NATIVE_DOUBLE, *dataspace));

    double input[dims[0]];
    for (int i = 0; i < dims[0]; i++) {
        input[i] = data(i);
    }

    dataset->write(input, PredType::NATIVE_DOUBLE);

    delete dataset;
    delete dataspace;
    delete group;
    delete file;
}

/*
inputs:
    string filename
    string groupname
    string dataset_name -> name the stored data
    hsize_t* dims -> 2 dimensional array with dimensions (x and y) of stored data
output:
    MatrixXd data_out -> data saved in "dataset_name"
*/
Eigen::MatrixXd read_matrix(const string &filename, const string &groupname, const string &dataset_name, const hsize_t* dims) {
    Eigen::MatrixXd data_out(dims[0], dims[1]);
    double out[dims[0]][dims[1]];
    
    H5File *file = new H5File(filename, H5F_ACC_RDONLY);
    Group *group = new Group(file->openGroup(groupname));
    DataSet *dataset = new DataSet(group->openDataSet(dataset_name));
    DataSpace dataspace = dataset->getSpace();

    dataset->read(out, PredType::NATIVE_DOUBLE, dataspace);

    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            data_out(i, j) = out[i][j];
        }
    }

    delete dataset;
    delete group;
    delete file;

    return data_out;
}

/*
inputs:
    string filename
    string groupname
    string dataset_name -> name the stored data
    hsize_t* dims -> 1 dimensional array with dimension (x) of stored data
output:
    VectorXd data_out -> data saved in "dataset_name"
*/
Eigen::VectorXd read_vector(const string &filename, const string &groupname, const string &dataset_name, const hsize_t* dims) {
    Eigen::VectorXd data_out(dims[0]);
    double out[dims[0]];
    
    H5File *file = new H5File(filename, H5F_ACC_RDONLY);
    Group *group = new Group(file->openGroup(groupname));
    DataSet *dataset = new DataSet(group->openDataSet(dataset_name));
    DataSpace dataspace = dataset->getSpace();

    dataset->read(out, PredType::NATIVE_DOUBLE, dataspace);

    for (int i = 0; i < dims[0]; i++) {
        data_out(i) = out[i];
    }

    delete dataset;
    delete group;
    delete file;

    return data_out;
}