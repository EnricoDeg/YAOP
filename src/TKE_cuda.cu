#include <iostream>

#include "cuda.h"
#include <cuda/std/mdspan>

#include "utils.hpp"
#include "TKE_cuda.hpp"
#include "cuda_check.hpp"

namespace cudastd = cuda::std;

constexpr auto dyn = cuda::std::dynamic_extent;
using ext2d_t = cuda::std::extents<size_t, dyn, dyn>;
using ext3d_t = cuda::std::extents<size_t, dyn, dyn, dyn>;
using mdspan_2d = cudastd::mdspan<double, ext2d_t>;
using mdspan_3d = cudastd::mdspan<double, ext3d_t>;

// this needs to be defined static or TKE_cuda.hpp can not 
// be included in C++ files (need to think about other solutions)
static mdspan_2d view_cuda_malloc(double *field, size_t dim1, size_t dim2);
static mdspan_3d view_cuda_malloc(double *field, size_t dim1, size_t dim2, size_t dim3);

// TKE internal memory views
static mdspan_2d rho_up_view;
static mdspan_3d tke_old_view;

// TKE interface memory views
static mdspan_3d tke_view;

// TKE CUDA kernels functions
__global__ void calc_impl_kernel();

TKE_cuda::TKE_cuda(int nproma, int nlevs, int nblocks, 
                   int block_size, int start_index, int end_index)
    : TKE_backend(nproma, nlevs, nblocks, block_size, start_index, end_index) {
    rho_up_view = view_cuda_malloc(rho_up, (size_t)nlevs, (size_t)nproma);
    tke_old_view = view_cuda_malloc(tke_old, (size_t)nproma, (size_t)(nlevs+1), (size_t)nblocks);
    is_view_init = false;
}

TKE_cuda::~TKE_cuda() {
    std::cout << "Finalizing TKE cuda... " << std::endl;
    check( cudaFree(rho_up) );
}

void TKE_cuda::calc_impl(int start_block, int end_block, double *tke) {
    if (!is_view_init) {
        mdspan_3d tke_view{ tke, ext3d_t{m_nblocks,m_nlevs,m_nproma} };
        is_view_init = true;
    }

    for (int jb=start_block; jb<=end_block; jb++) {
          int start_index, end_index;
        get_index_range(m_block_size, 0, m_nblocks-1, m_start_index, m_end_index,
                       jb, &start_index, &end_index);
    }

}

static mdspan_2d view_cuda_malloc(double *field, size_t dim1, size_t dim2) {

    check( cudaMalloc(&field, dim1*dim2*sizeof(double)) );
    mdspan_2d memview{ field, ext2d_t{dim1, dim2} };
    return memview;

}

static mdspan_3d view_cuda_malloc(double *field, size_t dim1, size_t dim2, size_t dim3) {

    check( cudaMalloc(&field, dim1*dim2*dim3*sizeof(double)) );
    mdspan_3d memview{ field, ext3d_t{dim1, dim2, dim3} };
    return memview;

}

__global__ void calc_impl_kernel() {

}

