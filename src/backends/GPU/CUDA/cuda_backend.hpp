/* Copyright (C) 2023  Enrico Degregori, Wilton Jaciel Loch
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef SRC_BACKENDS_GPU_CUDA_CUDA_BACKEND_HPP_
#define SRC_BACKENDS_GPU_CUDA_CUDA_BACKEND_HPP_

#include <cuda.h>
#include <cuda/std/mdspan>
#include "src/backends/GPU/CUDA/cuda_check.hpp"
#include "src/shared/interface/data_struct.hpp"

constexpr auto dyn = cuda::std::dynamic_extent;
using ext1d_d = cuda::std::extents<int, dyn>;
using ext2d_d = cuda::std::extents<int, dyn, dyn>;
using ext3d_d = cuda::std::extents<int, dyn, dyn, dyn>;
using ext1d_t = cuda::std::dextents<int, 1>;
using ext2d_t = cuda::std::dextents<int, 2>;
using ext3d_t = cuda::std::dextents<int, 3>;
using mdspan_1d_double = cuda::std::mdspan<double, ext1d_t>;
using mdspan_2d_double = cuda::std::mdspan<double, ext2d_t>;
using mdspan_3d_double = cuda::std::mdspan<double, ext3d_t>;
using mdspan_2d_int = cuda::std::mdspan<int, ext2d_t>;
using mdspan_3d_int = cuda::std::mdspan<int, ext3d_t>;

/*! \brief CUDA mdspan memory view policy.
 *
 *  It defines the policy to allocate/deallocate arrays and create CUDA mdspan objects.
 */
class cuda_mdspan_impl {
 public:
    /*! \brief Create a 1D mdspan object from double pointer.
     *
     */
    static mdspan_1d_double memview(double *data, int nlevs) {
        return mdspan_1d_double{ data, ext1d_d{nlevs} };
    }
    /*! \brief Create a 2D mdspan object from double pointer.
     *
     */
    static mdspan_2d_double memview(double *data, int nblocks, int nproma) {
        return mdspan_2d_double{ data, ext2d_d{nblocks, nproma} };
    }
    /*! \brief Create a 3D mdspan object from double pointer.
     *
     */
    static mdspan_3d_double memview(double *data, int nblocks, int nlevs, int nproma) {
        return mdspan_3d_double{ data, ext3d_d{nblocks, nlevs, nproma} };
    }
    /*! \brief Create a 2D mdspan object from int pointer.
     *
     */
    static mdspan_2d_int memview(int *data, int nblocks, int nproma) {
        return mdspan_2d_int{ data, ext2d_d{nblocks, nproma} };
    }
    /*! \brief Create a 3D mdspan object from int pointer.
     *
     */
    static mdspan_3d_int memview(int *data, int nblocks, int nlevs, int nproma) {
        return mdspan_3d_int{ data, ext3d_d{nblocks, nlevs, nproma} };
    }
    /*! \brief Allocate memory and create a 1D mdspan object from double pointer.
     *
     */
    static mdspan_1d_double memview_malloc(double *field, int dim1) {
        check( cudaMalloc(&field, dim1*sizeof(double)) );
        mdspan_1d_double memview{ field, ext1d_d{dim1} };
        return memview;
    }
    /*! \brief Allocate memory and create a 2D mdspan object from double pointer.
     *
     */
    static mdspan_2d_double memview_malloc(double *field, int dim1, int dim2) {
        check( cudaMalloc(&field, dim1*dim2*sizeof(double)) );
        mdspan_2d_double memview{ field, ext2d_d{dim1, dim2} };
        return memview;
    }
    /*! \brief Allocate memory and create a 3D mdspan object from double pointer.
     *
     */
    static mdspan_3d_double memview_malloc(double *field, int dim1, int dim2, int dim3) {
        check( cudaMalloc(&field, dim1*dim2*dim3*sizeof(double)) );
        mdspan_3d_double memview{ field, ext3d_d{dim1, dim2, dim3} };
        return memview;
    }
    /*! \brief Free memory from double pointer.
     *
     */
    static void memview_free(double *field) {
        check(cudaFree(field));
    }
};

/*! \brief CUDA kernel launch policy.
 *
 */
class cuda_launch_impl {
 public:
    /*! \brief It creates a CUDA kernel launch configuration object and launch the kernel with it.
     *
     */
    static void launch(int threadsPerBlock, int blocksPerGrid, void* func, void **args) {
        dim3 blocksPerGrid3(blocksPerGrid, 1, 1);
        dim3 threadsPerBlock3(threadsPerBlock, 1, 1);
        cudaLaunchConfig_t config = {0};
        config.gridDim = blocksPerGrid3;
        config.blockDim = threadsPerBlock3;
        check(cudaLaunchKernelExC(&config, func, args));
    }
};

namespace gpu_memview = cuda::std;
using gpu_memview_policy = cuda_mdspan_impl;
using gpu_launch_policy = cuda_launch_impl;

#endif  // SRC_BACKENDS_GPU_CUDA_CUDA_BACKEND_HPP_
