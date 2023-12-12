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

#ifndef SRC_BACKENDS_CPU_CPU_MEMORY_HPP_
#define SRC_BACKENDS_CPU_CPU_MEMORY_HPP_

#include <mdspan/mdspan.hpp>
#include "src/shared/interface/data_struct.hpp"

constexpr auto dyn = Kokkos::dynamic_extent;
using ext1d_d = Kokkos::extents<int, dyn>;
using ext2d_d = Kokkos::extents<int, dyn, dyn>;
using ext3d_d = Kokkos::extents<int, dyn, dyn, dyn>;
using ext1d_t = Kokkos::dextents<int, 1>;
using ext2d_t = Kokkos::dextents<int, 2>;
using ext3d_t = Kokkos::dextents<int, 3>;
using mdspan_1d_double = Kokkos::mdspan<double, ext1d_t>;
using mdspan_2d_double = Kokkos::mdspan<double, ext2d_t>;
using mdspan_3d_double = Kokkos::mdspan<double, ext3d_t>;
using mdspan_2d_int = Kokkos::mdspan<int, ext2d_t>;
using mdspan_3d_int = Kokkos::mdspan<int, ext3d_t>;

class cpu_mdspan_impl {
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
        field = reinterpret_cast<double *>(malloc(dim1 * sizeof(double)));
        mdspan_1d_double memview{ field, ext1d_d{dim1} };
        return memview;
    }
    static mdspan_2d_double memview_malloc(double *field, int dim1, int dim2) {
        field = reinterpret_cast<double *>(malloc(dim1 * dim2 * sizeof(double)));
        mdspan_2d_double memview{ field, ext2d_d{dim1, dim2} };
        return memview;
    }
    static mdspan_3d_double memview_malloc(double *field, int dim1, int dim2, int dim3) {
        field = reinterpret_cast<double *>(malloc(dim1 * dim2 * dim3 * sizeof(double)));
        mdspan_3d_double memview{ field, ext3d_d{dim1, dim2, dim3} };
        return memview;
    }
    static void memview_free(double *field) {
        free(field);
    }
};

#endif  // SRC_BACKENDS_CPU_CPU_MEMORY_HPP_
