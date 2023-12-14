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

#include <algorithm>
#include <iostream>
#include "src/backends/GPU/TKE_gpu.hpp"
#include "src/shared/utils.hpp"
#include "src/backends/GPU/gpu_kernels.hpp"

// Structures with memory views
static struct t_cvmix_view<gpu_memview::mdspan, gpu_memview::dextents> p_cvmix_view;
static struct t_patch_view<gpu_memview::mdspan, gpu_memview::dextents> p_patch_view;
static struct t_ocean_state_view<gpu_memview::mdspan, gpu_memview::dextents> ocean_state_view;
static struct t_atmo_fluxes_view<gpu_memview::mdspan, gpu_memview::dextents> atmos_fluxes_view;
static struct t_atmos_for_ocean_view<gpu_memview::mdspan, gpu_memview::dextents> p_as_view;
static struct t_sea_ice_view<gpu_memview::mdspan, gpu_memview::dextents> p_sea_ice_view;
static struct t_tke_internal_view<gpu_memview::mdspan, gpu_memview::dextents> p_internal_view;

TKE_gpu::TKE_gpu(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                   int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
                   int l_lc, double clc, double ReferencePressureIndbars, double pi)
    : TKE_backend(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
                  vert_cor_type, dtime, OceanReferenceDensity, grav,
                  l_lc, clc, ReferencePressureIndbars, pi) {
    // Allocate internal arrays memory and create memory views
    std::cout << "Initializing TKE (GPU)... " << std::endl;

    this->internal_fields_malloc<gpu_memview::mdspan, gpu_memview::dextents, gpu_memview_policy>
                                (&p_internal_view);
}

TKE_gpu::~TKE_gpu() {
    // Free internal arrays memory
    std::cout << "Finalizing TKE (GPU)... " << std::endl;

    this->internal_fields_free<gpu_memview_policy>();
}

void TKE_gpu::calc_impl(t_patch p_patch, t_cvmix p_cvmix,
                         t_ocean_state ocean_state, t_atmo_fluxes atmos_fluxes,
                         t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
                         int edges_block_size, int edges_start_block, int edges_end_block,
                         int edges_start_index, int edges_end_index, int cells_block_size,
                         int cells_start_block, int cells_end_block, int cells_start_index,
                         int cells_end_index) {
    // The pointer to the data should not change inside the time loop
    // structs view are filled only at the first time step
    if (!m_is_view_init) {
        this->fill_struct_memview<gpu_memview::mdspan, gpu_memview::dextents, gpu_memview_policy>
                                 (&p_cvmix_view, &p_cvmix, p_constant.nblocks, p_constant.nlevs, p_constant.nproma);
        this->fill_struct_memview<gpu_memview::mdspan, gpu_memview::dextents, gpu_memview_policy>
                                 (&p_patch_view, &p_patch, p_constant.nblocks, p_constant.nlevs, p_constant.nproma);
        this->fill_struct_memview<gpu_memview::mdspan, gpu_memview::dextents, gpu_memview_policy>
                                 (&ocean_state_view, &ocean_state, p_constant.nblocks, p_constant.nlevs,
                                  p_constant.nproma);
        this->fill_struct_memview<gpu_memview::mdspan, gpu_memview::dextents, gpu_memview_policy>
                                 (&atmos_fluxes_view, &atmos_fluxes, p_constant.nblocks, p_constant.nproma);
        this->fill_struct_memview<gpu_memview::mdspan, gpu_memview::dextents, gpu_memview_policy>
                                 (&p_as_view, &p_as, p_constant.nblocks, p_constant.nproma);
        this->fill_struct_memview<gpu_memview::mdspan, gpu_memview::dextents, gpu_memview_policy>
                                 (&p_sea_ice_view, &p_sea_ice, p_constant.nblocks, p_constant.nproma);
        m_is_view_init = true;
    }

    // over cells
    for (int jb = cells_start_block; jb <= cells_end_block; jb++) {
        int start_index, end_index;
        get_index_range(cells_block_size, cells_start_block, cells_end_block,
                        cells_start_index, cells_end_index, jb, &start_index, &end_index);
        int threadsPerBlockI = 512;  // too many registers used for 1024
        int blocksPerGridI = (end_index - start_index) / threadsPerBlockI + 1;
        void *args[] = {&jb, &start_index, &end_index,
                        &p_patch_view, &p_cvmix_view,
                        &ocean_state_view, &atmos_fluxes_view,
                        &p_as_view, &p_sea_ice_view,
                        &p_internal_view, &p_constant,
                        &p_constant_tke};
        this->launch_kernel<gpu_launch_policy>(threadsPerBlockI, blocksPerGridI,
                                              reinterpret_cast<void *>(calc_impl_cells), args);
    }

    // over edges
    for (int jb = edges_start_block; jb <= edges_end_block; jb++) {
        int start_index, end_index;
        get_index_range(edges_block_size, edges_start_block, edges_end_block,
                        edges_start_index, edges_end_index, jb, &start_index, &end_index);
        int threadsPerBlockI = 1024;
        int blocksPerGridI = (end_index - start_index) / threadsPerBlockI + 1;
        void *args[] = {&jb, &start_index, &end_index,
                        &p_patch_view, &p_cvmix_view,
                        &p_internal_view, &p_constant};
        this->launch_kernel<gpu_launch_policy>(threadsPerBlockI, blocksPerGridI,
                                              reinterpret_cast<void *>(calc_impl_edges), args);
    }
}
