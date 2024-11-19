#ifndef GINKGO_EXEC_H
#define GINKGO_EXEC_H

#include <ginkgo/ginkgo.hpp>

// ReferenceExecutor for debugging
inline const auto exec = gko::ReferenceExecutor::create();
inline const auto gpu_exec = exec;



// OmpExecutor for CPUs
//inline const auto exec = gko::OmpExecutor::create();

// CudaExecutor for Nvidia GPUs
//inline const auto gpu_exec = gko::CudaExecutor::create(0, gko::OmpExecutor::create());

// HipExecutor for AMD GPUs
// const auto exec = gko::HipExecutor::create(0, gko::OmpExecutor::create());

// DpcppExecutor for Intel GPUs
// const auto exec = gko::DpcppExecutor::create(0, gko::OmpExecutor::create());


#endif
