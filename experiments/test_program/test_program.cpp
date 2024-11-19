#include <iostream>
#include <ginkgo/ginkgo.hpp>

int main() {
        // Initialize the CUDA executor with the default device (device 0)
        auto omp_exec = gko::OmpExecutor::create();
        auto exec = gko::CudaExecutor::create(0, omp_exec);

        std::cout << "Entering main\n";
        int N = 100;

        // Create a Ginkgo dense matrix with dimensions N x 1 on the CUDA executor
        auto test_vec = gko::matrix::Dense<double>::create(exec, gko::dim<2>(N, 1));
        std::cout << "Created gko matrix\n";

        // Copy the matrix to the host for initialization
        auto test_vec_host = gko::matrix::Dense<double>::create(omp_exec, gko::dim<2>(N, 1));
        double* test_vec_entries_host = test_vec_host->get_values();

        // Initialize the matrix values
        for (int i = 0; i < N; ++i) {
            test_vec_entries_host[i] = i * 3.1;
        }

        // Copy the initialized values back to the device
        test_vec->copy_from(test_vec_host.get());

    return 0;
}
