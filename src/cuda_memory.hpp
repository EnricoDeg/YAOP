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

#ifndef SRC_CUDA_MEMORY_HPP_
#define SRC_CUDA_MEMORY_HPP_

#include <cuda.h>
#include <cuda/std/mdspan>
#include "src/cuda_check.hpp"
#include "src/data_struct.hpp"

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

class cuda_mdspan_impl {
 public:
    static mdspan_1d_double memview(double *data, int nlevs) {
        return mdspan_1d_double{ data, ext1d_d{nlevs} };
    }
    static mdspan_2d_double memview(double *data, int nblocks, int nproma) {
        return mdspan_2d_double{ data, ext2d_d{nblocks, nproma} };
    }
    static mdspan_3d_double memview(double *data, int nblocks, int nlevs, int nproma) {
        return mdspan_3d_double{ data, ext3d_d{nblocks, nlevs, nproma} };
    }
    static mdspan_2d_int memview(int *data, int nblocks, int nproma) {
        return mdspan_2d_int{ data, ext2d_d{nblocks, nproma} };
    }
    static mdspan_3d_int memview(int *data, int nblocks, int nlevs, int nproma) {
        return mdspan_3d_int{ data, ext3d_d{nblocks, nlevs, nproma} };
    }
    static mdspan_1d_double memview_malloc(double *field, int dim1) {
        check( cudaMalloc(&field, dim1*sizeof(double)) );
        mdspan_1d_double memview{ field, ext1d_d{dim1} };
        return memview;
    }
    static mdspan_2d_double memview_malloc(double *field, int dim1, int dim2) {
        check( cudaMalloc(&field, dim1*dim2*sizeof(double)) );
        mdspan_2d_double memview{ field, ext2d_d{dim1, dim2} };
        return memview;
    }
    static mdspan_3d_double memview_malloc(double *field, int dim1, int dim2, int dim3) {
        check( cudaMalloc(&field, dim1*dim2*dim3*sizeof(double)) );
        mdspan_3d_double memview{ field, ext3d_d{dim1, dim2, dim3} };
        return memview;
    }
    static void memview_free(double *field) {
        check(cudaFree(field));
    }
};

#endif  // SRC_CUDA_MEMORY_HPP_
