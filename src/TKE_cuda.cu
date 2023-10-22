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
using mdspan_2d_double = cudastd::mdspan<double, ext2d_t>;
using mdspan_3d_double = cudastd::mdspan<double, ext3d_t>;
using mdspan_2d_int = cudastd::mdspan<int, ext2d_t>;

// this needs to be defined static or TKE_cuda.hpp can not 
// be included in C++ files (need to think about other solutions)
static mdspan_2d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2);
static mdspan_3d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2, size_t dim3);

// TKE internal memory views
static mdspan_2d_double rho_up_view;
static mdspan_3d_double tke_old_view;

// TKE interface memory views
static mdspan_3d_double tke_view;
static mdspan_2d_int dolic_c_view;

// TKE CUDA kernels functions
__global__ void calc_impl_kernel(int blockNo, int start_index, int end_index, 
                                 mdspan_2d_int dolic_c, mdspan_3d_double tke,
                                 mdspan_3d_double tke_old);

TKE_cuda::TKE_cuda(int nproma, int nlevs, int nblocks, 
                   int block_size, int start_index, int end_index)
    : TKE_backend(nproma, nlevs, nblocks, block_size, start_index, end_index) {
    // Initialize internal arrays
    rho_up_view = view_cuda_malloc(rho_up, (size_t)nlevs, (size_t)nproma);
    tke_old_view = view_cuda_malloc(tke_old, (size_t)nproma, (size_t)(nlevs+1), (size_t)nblocks);
    is_view_init = false;
}

TKE_cuda::~TKE_cuda() {
    std::cout << "Finalizing TKE cuda... " << std::endl;
    check( cudaFree(rho_up) );
    check( cudaFree(tke_old) );
}

void TKE_cuda::calc_impl(int start_block, int end_block, double *tke, int *dolic_c) {
    if (!is_view_init) {
        tke_view = mdspan_3d_double{ tke, ext3d_t{m_nblocks,m_nlevs,m_nproma} };
        dolic_c_view = mdspan_2d_int{ dolic_c, ext2d_t{m_nblocks,m_nproma} };
        is_view_init = true;
    }

    for (int jb=start_block; jb<=end_block; jb++) {
        int start_index, end_index;
        get_index_range(m_block_size, 0, m_nblocks-1, m_start_index, m_end_index,
                       jb, &start_index, &end_index);
        int threadsPerBlockI = 512;
        int blocksPerGridI = (end_index - start_index) / threadsPerBlockI + 1;
        dim3 blocksPerGrid(blocksPerGridI, 1, 1);
        dim3 threadsPerBlock(threadsPerBlockI, 1, 1);
        calc_impl_kernel<<<blocksPerGrid,threadsPerBlock>>>(jb, start_index, end_index, 
                                                            dolic_c_view, tke_view, 
                                                            tke_old_view);
    }

}

static mdspan_2d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2) {

    check( cudaMalloc(&field, dim1*dim2*sizeof(double)) );
    mdspan_2d_double memview{ field, ext2d_t{dim1, dim2} };
    return memview;

}

static mdspan_3d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2, size_t dim3) {

    check( cudaMalloc(&field, dim1*dim2*dim3*sizeof(double)) );
    mdspan_3d_double memview{ field, ext3d_t{dim1, dim2, dim3} };
    return memview;

}

__global__ void calc_impl_kernel(int blockNo, int start_index, int end_index, 
                                 mdspan_2d_int dolic_c, mdspan_3d_double tke, 
                                 mdspan_3d_double tke_old) {

    int jc = blockIdx.x * blockDim.x + threadIdx.x + start_index;
    if (jc <= end_index) {
        int levels = dolic_c(blockNo,jc);
        for (int level = 0; level < levels; level++) {
            tke_old(blockNo,level,jc) = tke(blockNo,level,jc);
            tke(blockNo,level,jc) = tke(blockNo,level,jc) + 1.0;
        }
    }

}

