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
#include "src/shared/assertion.hpp"

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

/*! \brief CPU mdspan memory view policy.
 *
 *  It defines the policy to allocate/deallocate arrays and create Kokkos mdspan objects.
 */
template <class T>
class cpu_mdspan_impl {
 public:
    /*! \brief Create a 1D mdspan object from given pointer.
     *
     *  If a NULL pointer is provided, an empty mdspan is created (size 0)
     */
    static mdspan_1d<T> memview(T *data, int nlevs) {
        YAOP_ASSERT(nlevs > 0);

        if (data == NULL)
            return mdspan_1d<T>{ data, ext1d_d{0} };
        else
            return mdspan_1d<T>{ data, ext1d_d{nlevs} };
    }
    /*! \brief Create a 2D mdspan object from pointer.
     *
     *  If a NULL pointer is provided, an empty mdspan is created (size 0)
     */
    static mdspan_2d<T> memview(T *data, int nblocks, int nproma) {
        YAOP_ASSERT(nblocks > 0);
        YAOP_ASSERT(nproma > 0);

        if (data == NULL)
            return mdspan_2d<T>{ data, ext2d_d{0, 0} };
        else
            return mdspan_2d<T>{ data, ext2d_d{nblocks, nproma} };
    }
    /*! \brief Create a 3D mdspan object from pointer.
     *
     *  If a NULL pointer is provided, an empty mdspan is created (size 0)
     */
    static mdspan_3d<T> memview(T *data, int nblocks, int nlevs, int nproma) {
        YAOP_ASSERT(nblocks > 0);
        YAOP_ASSERT(nlevs > 0);
        YAOP_ASSERT(nproma > 0);

        if (data == NULL)
            return mdspan_3d<T>{ data, ext3d_d{0, 0, 0} };
        else
            return mdspan_3d<T>{ data, ext3d_d{nblocks, nlevs, nproma} };
    }
    /*! \brief Allocate memory and create a 1D mdspan object.
     *
     */
    static mdspan_1d<T> memview_malloc(T *field, int dim1) {
        YAOP_ASSERT(dim1 >= 0);

        field = reinterpret_cast<T *>(malloc(dim1 * sizeof(T)));
        return mdspan_1d<T>{ field, ext1d_d{dim1} };
    }
    /*! \brief Allocate memory and create a 2D mdspan object from double pointer.
     *
     */
    static mdspan_2d<T> memview_malloc(T *field, int dim1, int dim2) {
        YAOP_ASSERT(dim1 >= 0);
        YAOP_ASSERT(dim2 >= 0);

        field = reinterpret_cast<T *>(malloc(dim1 * dim2 * sizeof(T)));
        return mdspan_2d<T>{ field, ext2d_d{dim1, dim2} };
    }
    /*! \brief Allocate memory and create a 3D mdspan object from double pointer.
     *
     */
    static mdspan_3d<T> memview_malloc(T *field, int dim1, int dim2, int dim3) {
        YAOP_ASSERT(dim1 >= 0);
        YAOP_ASSERT(dim2 >= 0);
        YAOP_ASSERT(dim3 >= 0);

        field = reinterpret_cast<T *>(malloc(dim1 * dim2 * dim3 * sizeof(T)));
        return mdspan_3d<T>{ field, ext3d_d{dim1, dim2, dim3} };
    }
    /*! \brief Free memory from double pointer.
     *
     */
    static void memview_free(T *field) {
        free(field);
    }
};

namespace cpu_memview = Kokkos;
template<class T>
using cpu_memview_policy = cpu_mdspan_impl<T>;

#endif  // SRC_BACKENDS_CPU_CPU_MEMORY_HPP_
