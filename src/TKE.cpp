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

#include "src/TKE.hpp"

#include <iostream>

#include "src/TKE_backend.hpp"
#ifdef CUDA
#include "src/TKE_gpu.hpp"
#else
#include "src/TKE_cpu.hpp"
#endif

struct TKE::Impl {
  TKE_backend::Ptr backend;
};

TKE::TKE(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
         int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
         int l_lc, double clc, double ReferencePressureIndbars, double pi)
    : m_impl(new Impl) {
    std::cout << "Initializing TKE... " << std::endl;
#ifdef CUDA
    m_impl->backend = TKE_backend::Ptr(new TKE_gpu(nproma, nlevs, nblocks,
                                       vert_mix_type, vmix_idemix_tke, vert_cor_type,
                                       dtime, OceanReferenceDensity, grav, l_lc, clc,
                                       ReferencePressureIndbars, pi));
#else
    m_impl->backend = TKE_backend::Ptr(new TKE_cpu(nproma, nlevs, nblocks,
                                       vert_mix_type, vmix_idemix_tke, vert_cor_type,
                                       dtime, OceanReferenceDensity, grav, l_lc, clc,
                                       ReferencePressureIndbars, pi));
#endif
    m_is_struct_init = false;
}

TKE::~TKE() {
    std::cout << "Finalizing TKE... " << std::endl;
    delete m_impl;
}

void TKE::calc(double *depth_CellInterface, double *prism_center_dist_c,
               double *inv_prism_center_dist_c, double *prism_thick_c,
               int *dolic_c, int *dolic_e, double *zlev_i, double *wet_c,
               int *edges_cell_idx, int *edges_cell_blk,
               double *temp, double *salt, double *stretch_c, double *eta_c,
               double *p_vn_x1, double *p_vn_x2, double *p_vn_x3,
               double *tke, double *tke_plc_in, double *hlc_in, double *wlc_in,
               double *u_stokes_in, double *a_veloc_v, double *a_temp_v, double *a_salt_v,
               double *iwe_Tdis, double *cvmix_dummy_1, double *cvmix_dummy_2,
               double *cvmix_dummy_3, double *tke_Tbpr, double *tke_Tspr,
               double *tke_Tdif, double *tke_Tdis, double *tke_Twin,
               double *tke_Tiwf, double *tke_Tbck, double *tke_Ttot,
               double *tke_Lmix, double *tke_Pr, double *stress_xw,
               double *stress_yw, double *fu10, double *concsum,
               int edges_block_size, int edges_start_block, int edges_end_block,
               int edges_start_index, int edges_end_index, int cells_block_size,
               int cells_start_block, int cells_end_block, int cells_start_index,
               int cells_end_index) {
    if (!m_is_struct_init) {
        fill_struct(&p_patch, depth_CellInterface, prism_center_dist_c,
                    inv_prism_center_dist_c, prism_thick_c, dolic_c, dolic_e,
                    zlev_i, wet_c, edges_cell_idx, edges_cell_blk);
        fill_struct(&p_cvmix, tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v,
                    a_temp_v, a_salt_v, iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2,
                    cvmix_dummy_3, tke_Tbpr, tke_Tspr, tke_Tdif, tke_Tdis, tke_Twin,
                    tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr);
        fill_struct(&ocean_state, temp, salt, stretch_c, eta_c, p_vn_x1, p_vn_x2, p_vn_x3);
        fill_struct(&atmos_fluxes, stress_xw, stress_yw);
        fill_struct(&p_as, fu10);
        fill_struct(&p_sea_ice, concsum);
        m_is_struct_init = true;
    }

    m_impl->backend->calc(p_patch, p_cvmix, ocean_state, atmos_fluxes, p_as, p_sea_ice,
                          edges_block_size, edges_start_block, edges_end_block,
                          edges_start_index, edges_end_index, cells_block_size,
                          cells_start_block, cells_end_block, cells_start_index,
                          cells_end_index);
}
