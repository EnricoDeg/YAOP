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
TEST(calc_diffusivity, tke_Av_0D) {
    int nblocks = 1;
    int nproma = 1;
    int blockNo = 0;
    int start_index = 0;
    int end_index = 0;
    int max_levels = 0;
    int nlevs = 1;

    t_constant_tke p_constant_tke;
    p_constant_tke.KappaM_max = 1.0;
    p_constant_tke.c_k = 1.0;
    p_constant_tke.only_tke = true;
    p_constant_tke.use_Kappa_min = false;

    int *dolic_c_ptr = NULL;
    double *Lmix_ptr = NULL;
    double *sqrttke_ptr = NULL;
    double *Nsqr_ptr = NULL;
    double *Ssqr_ptr = NULL;
    double *tke_Av_ptr = NULL;
    double *tke_kv_ptr = NULL;
    double *tke_Pr_ptr = NULL;

    // Allocate memory and create memview objs
    mdspan_2d_int dolic_c = cpu_mdspan_impl::memview_malloc(dolic_c_ptr, nblocks, nproma);
    mdspan_3d_double tke_Lmix = cpu_mdspan_impl::memview_malloc(Lmix_ptr, nblocks, nlevs+1, nproma);
    mdspan_2d_double sqrttke = cpu_mdspan_impl::memview_malloc(sqrttke_ptr, nlevs+1, nproma);
    mdspan_2d_double Nsqr = cpu_mdspan_impl::memview_malloc(Nsqr_ptr, nlevs+1, nproma);
    mdspan_2d_double Ssqr = cpu_mdspan_impl::memview_malloc(Ssqr_ptr, nlevs+1, nproma);
    mdspan_3d_double tke_Av = cpu_mdspan_impl::memview_malloc(tke_Av_ptr, nblocks, nlevs+1, nproma);
    mdspan_2d_double tke_kv = cpu_mdspan_impl::memview_malloc(tke_kv_ptr, nlevs+1, nproma);
    mdspan_3d_double tke_Pr = cpu_mdspan_impl::memview_malloc(tke_Pr_ptr, nblocks, nlevs+1, nproma);

    // Initialize arrays
    for (int jb = 0; jb < nblocks; jb++)
        for (int jc = 0; jc < nproma; jc++)
            dolic_c(jb, jc) = max_levels;

    for (int jb = 0; jb < nblocks; jb++)
        for (int level = 0; level < nlevs+1; level++)
            for (int jc = 0; jc < nproma; jc++)
                tke_Lmix(jb, level, jc) = 0.0;

    for (int level = 0; level < nlevs+1; level++)
        for (int jc = 0; jc < nproma; jc++)
            sqrttke(level, jc) = 0.0;

    for (int jb = 0; jb < nblocks; jb++)
        for (int level = 0; level < nlevs+1; level++)
            for (int jc = 0; jc < nproma; jc++)
                tke_Av(jb, level, jc) = 1.0;

    // calculate diffusivities
    calc_diffusivity(blockNo, start_index, end_index, max_levels,
                     &p_constant_tke,
                     dolic_c, tke_Lmix, sqrttke, Nsqr, Ssqr,
                     tke_Av, tke_kv, tke_Pr);

    // checks
    ASSERT_EQ(tke_Av(0, 0, 0), 0.0);
    ASSERT_EQ(tke_Av(0, 1, 0), 1.0);
}
