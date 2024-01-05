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

#ifndef SRC_BACKENDS_GPU_HIP_HIP_BACKEND_HPP_
#define SRC_BACKENDS_GPU_HIP_HIP_BACKEND_HPP_

#include <hip/hip_runtime.h>
#include <mdspan/mdspan.hpp>
#include "src/backends/GPU/HIP/hip_check.hpp"

constexpr auto dyn = Kokkos::dynamic_extent;
using ext1d_d = Kokkos::extents<int, dyn>;
using ext2d_d = Kokkos::extents<int, dyn, dyn>;
using ext3d_d = Kokkos::extents<int, dyn, dyn, dyn>;
using ext1d_t = Kokkos::dextents<int, 1>;
using ext2d_t = Kokkos::dextents<int, 2>;
using ext3d_t = Kokkos::dextents<int, 3>;

template<class T>
using mdspan_1d = Kokkos::mdspan<T, ext1d_t>;
template<class T>
using mdspan_2d = Kokkos::mdspan<T, ext2d_t>;
template<class T>
using mdspan_3d = Kokkos::mdspan<T, ext3d_t>;

/*! \brief HIP mdspan memory view policy.
 *
 *  It defines the policy to allocate/deallocate arrays and create Kokkos mdspan objects.
 */
template <class T>
class hip_mdspan_impl {
 public:
    /*! \brief Create a 1D mdspan object from given pointer.
     *
     */
    static mdspan_1d<T> memview(T *data, int nlevs) {
        return mdspan_1d<T>{ data, ext1d_d{nlevs} };
    }
    /*! \brief Create a 2D mdspan object from given pointer.
     *
     */
    static mdspan_2d<T> memview(T *data, int nblocks, int nproma) {
        return mdspan_2d<T>{ data, ext2d_d{nblocks, nproma} };
    }
    /*! \brief Create a 3D mdspan object from given pointer.
     *
     */
    static mdspan_3d<T> memview(T *data, int nblocks, int nlevs, int nproma) {
        return mdspan_3d<T>{ data, ext3d_d{nblocks, nlevs, nproma} };
    }
    /*! \brief Allocate memory and create a 1D mdspan object.
     *
     */
    static mdspan_1d<T> memview_malloc(T *field, int dim1) {
        check(hipMalloc(reinterpret_cast<void**>(field), dim1 * sizeof(T)));
        mdspan_1d<T> memview{ field, ext1d_d{dim1} };
        return memview;
    }
    /*! \brief Allocate memory and create a 2D mdspan object.
     *
     */
    static mdspan_2d<T> memview_malloc(T *field, int dim1, int dim2) {
        check(hipMalloc(reinterpret_cast<void**>(field), dim1 * dim2 * sizeof(T)));
        mdspan_2d<T> memview{ field, ext2d_d{dim1, dim2} };
        return memview;
    }
    /*! \brief Allocate memory and create a 3D mdspan object.
     *
     */
    static mdspan_3d<T> memview_malloc(T *field, int dim1, int dim2, int dim3) {
        check(hipMalloc(reinterpret_cast<void**>(field), dim1 * dim2 * dim3 * sizeof(T)));
        mdspan_3d<T> memview{ field, ext3d_d{dim1, dim2, dim3} };
        return memview;
    }
    /*! \brief Free memory from double pointer.
     *
     */
    static void memview_free(T *field) {
        check(hipFree(field));
    }
};

/*! \brief HIP kernel launch policy.
 *
 */
class hip_launch_impl {
 public:
    /*! \brief It calls the HIP function to launch the kernel function.
     *
     */
    static void launch(int threadsPerBlock, int blocksPerGrid, void* func, void **args) {
        dim3 blocksPerGrid3(blocksPerGrid, 1, 1);
        dim3 threadsPerBlock3(threadsPerBlock, 1, 1);
        check(hipLaunchKernelGGL(func, blocksPerGrid3, threadsPerBlock3, 0, 0, args));
    }
};

namespace memview_nms = Kokkos;
using gpu_launch_policy = hip_launch_impl;
template<class T>
using gpu_memview_policy = hip_mdspan_impl<T>;

#endif  // SRC_BACKENDS_GPU_HIP_HIP_BACKEND_HPP_
