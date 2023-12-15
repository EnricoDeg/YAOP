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

#include <gtest/gtest.h>
#include "src/backends/CPU/cpu_kernels.hpp"

// Test single point with single vertical level
TEST(calculate_mxl_2, min_val_0D) {
    int nblocks = 1;
    int nproma = 1;
    int blockNo = 0;
    int start_index = 0;
    int end_index = 0;
    int max_levels = 1;
    double mxl_min = 0.01;

    int *dolic_c_ptr;
    double *Lmix_ptr;
    double *dzw_ptr;

    // Allocate memory and create memview objs
    mdspan_2d_int dolic_c = cpu_mdspan_impl::memview_malloc(dolic_c_ptr, nblocks, nproma);
    mdspan_3d_double tke_Lmix = cpu_mdspan_impl::memview_malloc(Lmix_ptr, nblocks, max_levels+1, nproma);
    mdspan_2d_double dzw_stretched = cpu_mdspan_impl::memview_malloc(dzw_ptr, max_levels, nproma);

    // Initialize arrays
    for (int jb = 0; jb < nblocks; jb++)
        for (int jc = 0; jc < nproma; jc++)
            dolic_c(jb, jc) = max_levels;

    for (int jb = 0; jb < nblocks; jb++)
        for (int level = 0; level < max_levels+1; level++)
            for (int jc = 0; jc < nproma; jc++)
                tke_Lmix(jb, level, jc) = 0.0;

    for (int level = 0; level < max_levels; level++)
        for (int jc = 0; jc < nproma; jc++)
            dzw_stretched(level, jc) = 1.0;

    // compute mixing length scale
    calculate_mxl_2(blockNo, start_index, end_index, max_levels, mxl_min,
                    dolic_c, tke_Lmix, dzw_stretched);

    // checks
    for (int jb = 0; jb < nblocks; jb++)
        for (int level = 0; level < max_levels+1; level++)
            for (int jc = 0; jc < nproma; jc++)
                ASSERT_EQ(tke_Lmix(jb, level, jc), mxl_min);
}
