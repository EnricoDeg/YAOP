#include <iostream>
#include "cuda.h"
#include <cuda/std/mdspan>
#include "TKE_cuda.hpp"

namespace stdex = cuda::std;

TKE_cuda::TKE_cuda(int nproma, int nlevs, int nblocks)
    : TKE_backend(nproma, nlevs, nblocks) {
    cudaMalloc(&rho_up, nproma*nlevs*sizeof(double));
}

TKE_cuda::~TKE_cuda() {
    std::cout << "Finalizing TKE cuda... " << std::endl;
    cudaFree(rho_up);
}

void TKE_cuda::calc_impl() {

}
