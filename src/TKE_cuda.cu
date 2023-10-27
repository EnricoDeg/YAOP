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

#include <iostream>
#include "src/TKE_cuda.hpp"
#include "src/utils.hpp"
#include "src/cuda_memory.hpp"

// Structures with memory views
struct t_cvmix_view p_cvmix_view_l;
struct t_patch_view p_patch_view_l;
struct t_ocean_state_view ocean_state_view_l;
struct t_atmo_fluxes_view atmos_fluxes_view_l;
struct t_atmos_for_ocean_view p_as_view_l;
struct t_sea_ice_view p_sea_ice_view_l;
struct t_tke_internal_view p_internal_view_l;

// TKE CUDA kernels functions
__global__ void calc_impl_kernel(int blockNo, int start_index, int end_index,
                                 t_patch_view p_patch, t_cvmix_view p_cvmix,
                                 t_ocean_state_view ocean_state, t_atmo_fluxes_view atmos_fluxes,
                                 t_atmos_for_ocean_view p_as, t_sea_ice_view p_sea_ice,
                                 t_tke_internal_view p_internal);

TKE_cuda::TKE_cuda(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                   int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
                   int l_lc, double clc, double ReferencePressureIndbars, double pi)
    : TKE_backend(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
                  vert_cor_type, dtime, OceanReferenceDensity, grav,
                  l_lc, clc, ReferencePressureIndbars, pi) {
    // Allocate internal arrays memory and create memory views
    std::cout << "Initializing TKE cuda... " << std::endl;
    p_internal_view_l.tke_old = view_cuda_malloc(m_tke_old, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.rho_up = view_cuda_malloc(m_rho_up, static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    p_internal_view_l.rho_down = view_cuda_malloc(m_rho_down, static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    p_internal_view_l.forc_tke_surf_2D = view_cuda_malloc(m_forc_tke_surf_2D, static_cast<size_t>(nproma));
    p_internal_view_l.forc_rho_surf_2D = view_cuda_malloc(m_forc_rho_surf_2D, static_cast<size_t>(nproma));
    p_internal_view_l.bottom_fric_2D = view_cuda_malloc(m_bottom_fric_2D, static_cast<size_t>(nproma));
    p_internal_view_l.s_c = view_cuda_malloc(m_s_c, static_cast<size_t>(nproma));
    p_internal_view_l.dzw_stretched = view_cuda_malloc(m_dzw_stretched,
                                                       static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    p_internal_view_l.dzt_stretched = view_cuda_malloc(m_dzt_stretched,
                                                       static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.tke_Av = view_cuda_malloc(m_tke_Av, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.tke_kv = view_cuda_malloc(m_tke_kv, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.tke_iw_alpha_c = view_cuda_malloc(m_tke_iw_alpha_c,
                                                        static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.tke_iwe = view_cuda_malloc(m_tke_iwe, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.tke_iwe_forcing = view_cuda_malloc(m_tke_iwe_forcing, static_cast<size_t>(nlevs+1),
                                            static_cast<size_t>(nproma));
    p_internal_view_l.pressure = view_cuda_malloc(m_pressure, static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    p_internal_view_l.Nsqr = view_cuda_malloc(m_Nsqr, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.Ssqr = view_cuda_malloc(m_Ssqr, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    is_view_init = false;
}

TKE_cuda::~TKE_cuda() {
    // Free internal arrays memory
    std::cout << "Finalizing TKE cuda... " << std::endl;
    check(cudaFree(m_rho_up) );
    check(cudaFree(m_rho_down) );
    check(cudaFree(m_forc_tke_surf_2D) );
    check(cudaFree(m_forc_rho_surf_2D) );
    check(cudaFree(m_bottom_fric_2D) );
    check(cudaFree(m_s_c) );
    check(cudaFree(m_dzw_stretched) );
    check(cudaFree(m_dzt_stretched) );
    check(cudaFree(m_tke_old) );
    check(cudaFree(m_tke_Av) );
    check(cudaFree(m_tke_kv) );
    check(cudaFree(m_tke_iw_alpha_c) );
    check(cudaFree(m_tke_iwe) );
    check(cudaFree(m_tke_iwe_forcing) );
    check(cudaFree(m_pressure) );
    check(cudaFree(m_Nsqr) );
    check(cudaFree(m_Ssqr) );
}

void TKE_cuda::calc_impl(t_patch p_patch, t_cvmix p_cvmix,
                         t_ocean_state ocean_state, t_atmo_fluxes atmos_fluxes,
                         t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
                         int edges_block_size, int edges_start_block, int edges_end_block,
                         int edges_start_index, int edges_end_index, int cells_block_size,
                         int cells_start_block, int cells_end_block, int cells_start_index,
                         int cells_end_index) {
    if (!is_view_init) {
        fill_struct_view(&p_cvmix_view_l, &p_cvmix, m_nblocks, m_nlevs, m_nproma);
        fill_struct_view(&p_patch_view_l, &p_patch, m_nblocks, m_nlevs, m_nproma);
        fill_struct_view(&ocean_state_view_l, &ocean_state, m_nblocks, m_nlevs, m_nproma);
        fill_struct_view(&atmos_fluxes_view_l, &atmos_fluxes, m_nblocks, m_nlevs, m_nproma);
        fill_struct_view(&p_as_view_l, &p_as, m_nblocks, m_nlevs, m_nproma);
        fill_struct_view(&p_sea_ice_view_l, &p_sea_ice, m_nblocks, m_nlevs, m_nproma);
        is_view_init = true;
    }

    for (int jb = cells_start_block; jb <= cells_end_block; jb++) {
        int start_index, end_index;
        get_index_range(cells_block_size, cells_start_block, cells_end_block,
                        cells_start_index, cells_end_index, jb, &start_index, &end_index);
        int threadsPerBlockI = 512;
        int blocksPerGridI = (end_index - start_index) / threadsPerBlockI + 1;
        dim3 blocksPerGrid(blocksPerGridI, 1, 1);
        dim3 threadsPerBlock(threadsPerBlockI, 1, 1);
        calc_impl_kernel<<<blocksPerGrid, threadsPerBlock>>>(jb, start_index, end_index,
                                                            p_patch_view_l, p_cvmix_view_l,
                                                            ocean_state_view_l, atmos_fluxes_view_l,
                                                            p_as_view_l, p_sea_ice_view_l,
                                                            p_internal_view_l);
    }
}

__global__
void calc_impl_kernel(int blockNo, int start_index, int end_index, t_patch_view p_patch,
                      t_cvmix_view p_cvmix, t_ocean_state_view ocean_state,
                      t_atmo_fluxes_view atmos_fluxes,
                      t_atmos_for_ocean_view p_as, t_sea_ice_view p_sea_ice,
                      t_tke_internal_view p_internal) {
    int jc = blockIdx.x * blockDim.x + threadIdx.x + start_index;
    if (jc <= end_index) {
        int levels = p_patch.dolic_c(blockNo, jc);
        for (int level = 0; level < levels; level++) {
            p_internal.tke_old(level, jc) = p_cvmix.tke(blockNo, level, jc);
            p_cvmix.tke(blockNo, level, jc) = p_cvmix.tke(blockNo, level, jc) + 1.0;
        }
    }
}
