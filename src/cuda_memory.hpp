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
using ext1d_t = cuda::std::extents<size_t, dyn>;
using ext2d_t = cuda::std::extents<size_t, dyn, dyn>;
using ext3d_t = cuda::std::extents<size_t, dyn, dyn, dyn>;
using mdspan_1d_double = cuda::std::mdspan<double, ext1d_t>;
using mdspan_2d_double = cuda::std::mdspan<double, ext2d_t>;
using mdspan_3d_double = cuda::std::mdspan<double, ext3d_t>;
using mdspan_2d_int = cuda::std::mdspan<int, ext2d_t>;
using mdspan_3d_int = cuda::std::mdspan<int, ext3d_t>;

// TKE interface memory views
struct t_patch_view {
    mdspan_3d_double depth_CellInterface;
    mdspan_3d_double prism_center_dist_c;
    mdspan_3d_double inv_prism_center_dist_c;
    mdspan_3d_double prism_thick_c;
    mdspan_2d_int dolic_c;
    mdspan_2d_int dolic_e;
    mdspan_1d_double zlev_i;
    mdspan_3d_double wet_c;
    mdspan_3d_int edges_cell_idx;
    mdspan_3d_int edges_cell_blk;
};

struct t_ocean_state_view {
    mdspan_3d_double temp;
    mdspan_3d_double salt;
    mdspan_2d_double stretch_c;
    mdspan_2d_double eta_c;
    mdspan_3d_double p_vn_x1;
    mdspan_3d_double p_vn_x2;
    mdspan_3d_double p_vn_x3;
};

void fill_struct_view(struct t_patch_view *p_patch_view, struct t_patch *p_patch,
                             int nblocks, int nlevs, int nproma);


void fill_struct_view(struct t_ocean_state_view *ocean_state_view, struct t_ocean_state *ocean_state,
                             int nblocks, int nlevs, int nproma);

class cuda_mdspan_impl {
 public:
    static mdspan_2d_double memview(double *data, int nblocks, int nproma) {
        return mdspan_2d_double{ data, ext2d_t{nblocks, nproma} };
    }
    static mdspan_3d_double memview(double *data, int nblocks, int nlevs, int nproma) {
        return mdspan_3d_double{ data, ext3d_t{nblocks, nlevs, nproma} };
    }
    static mdspan_1d_double memview_malloc(double *field, int dim1) {
        check( cudaMalloc(&field, dim1*sizeof(double)) );
        mdspan_1d_double memview{ field, ext1d_t{static_cast<size_t>(dim1)} };
        return memview;
    }
    static mdspan_2d_double memview_malloc(double *field, int dim1, int dim2) {
        check( cudaMalloc(&field, dim1*dim2*sizeof(double)) );
        mdspan_2d_double memview{ field, ext2d_t{static_cast<size_t>(dim1),
                                                 static_cast<size_t>(dim2)} };
        return memview;
    }
    static mdspan_3d_double memview_malloc(double *field, size_t dim1, size_t dim2, size_t dim3) {
        check( cudaMalloc(&field, dim1*dim2*dim3*sizeof(double)) );
        mdspan_3d_double memview{ field, ext3d_t{static_cast<size_t>(dim1),
                                                 static_cast<size_t>(dim2),
                                                 static_cast<size_t>(dim3)} };
        return memview;
    }
    static void memview_free(double *field) {
        check(cudaFree(field));
    }
};

#endif  // SRC_CUDA_MEMORY_HPP_
