#include <iostream>
#include "cuda.h"
#include <cuda/std/mdspan>
#include "TKE_cuda.hpp"

namespace stdex = cuda::std;

TKE_cuda::TKE_cuda(int nproma, int nlevs, int nblocks)
    : TKE_backend(nproma, nlevs, nblocks) {
    cudaMalloc(&rho_up, nproma*nlevs*sizeof(double));
}

void TKE_cuda::calc_impl() {

}
