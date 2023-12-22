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

#include "src/YAOP.hpp"

#include <iostream>

#include "src/backends/TKE_backend.hpp"
#ifdef CUDA
#include "src/backends/GPU/TKE_gpu.hpp"
#else
#include "src/backends/CPU/TKE_cpu.hpp"
#endif

struct YAOP::Impl {
  TKE_backend<double>::Ptr backend_tke_dp;
  TKE_backend<float>::Ptr backend_tke_sp;
};

YAOP::YAOP(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
         int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
         int l_lc, double clc, double ReferencePressureIndbars, double pi)
    : m_impl(new Impl) {
    std::cout << "Initializing Ocean Physics Library ... " << std::endl;
#ifdef CUDA
    m_impl->backend_tke_dp = TKE_backend<double>::Ptr(new TKE_gpu(nproma, nlevs, nblocks,
                                       vert_mix_type, vmix_idemix_tke, vert_cor_type,
                                       dtime, OceanReferenceDensity, grav, l_lc, clc,
                                       ReferencePressureIndbars, pi));
#else
    m_impl->backend_tke_dp = TKE_backend<double>::Ptr(new TKE_cpu<double>(nproma, nlevs, nblocks,
                                       vert_mix_type, vmix_idemix_tke, vert_cor_type,
                                       dtime, OceanReferenceDensity, grav, l_lc, clc,
                                       ReferencePressureIndbars, pi));
#endif
    m_is_struct_init = false;
}

YAOP::YAOP(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
         int vert_cor_type, float dtime, float OceanReferenceDensity, float grav,
         int l_lc, float clc, float ReferencePressureIndbars, float pi)
    : m_impl(new Impl) {
    std::cout << "Initializing Ocean Physics Library ... " << std::endl;
#ifdef CUDA
    m_impl->backend_tke_sp = TKE_backend<float>::Ptr(new TKE_gpu(nproma, nlevs, nblocks,
                                       vert_mix_type, vmix_idemix_tke, vert_cor_type,
                                       dtime, OceanReferenceDensity, grav, l_lc, clc,
                                       ReferencePressureIndbars, pi));
#else
    m_impl->backend_tke_sp = TKE_backend<float>::Ptr(new TKE_cpu<float>(nproma, nlevs, nblocks,
                                       vert_mix_type, vmix_idemix_tke, vert_cor_type,
                                       dtime, OceanReferenceDensity, grav, l_lc, clc,
                                       ReferencePressureIndbars, pi));
#endif
    m_is_struct_init = false;
}

YAOP::~YAOP() {
    std::cout << "Finalizing Ocean Physics Library ... " << std::endl;
    delete m_impl;
}

void YAOP::calc_tke(double *depth_CellInterface, double *prism_center_dist_c,
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
    // This structures are filled during the first call to initialize the memory views.
    // Then they are left empty during the following calls, but passed anyway to the backend.
    // It allows to have the struct templated as well as the backend without having to template
    // the YAOP class.
    // There can be more elegant workarounds to achieve that but this is the simplest one.
    struct t_patch<double> p_patch;
    struct t_cvmix<double> p_cvmix;
    struct t_sea_ice<double> p_sea_ice;
    struct t_atmos_for_ocean<double> p_as;
    struct t_atmo_fluxes<double> atmos_fluxes;
    struct t_ocean_state<double> ocean_state;

    if (!m_is_struct_init) {
        fill_struct<double>(&p_patch, depth_CellInterface, prism_center_dist_c,
                    inv_prism_center_dist_c, prism_thick_c, dolic_c, dolic_e,
                    zlev_i, wet_c, edges_cell_idx, edges_cell_blk);
        fill_struct<double>(&p_cvmix, tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v,
                    a_temp_v, a_salt_v, iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2,
                    cvmix_dummy_3, tke_Tbpr, tke_Tspr, tke_Tdif, tke_Tdis, tke_Twin,
                    tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr);
        fill_struct<double>(&ocean_state, temp, salt, stretch_c, eta_c, p_vn_x1, p_vn_x2, p_vn_x3);
        fill_struct<double>(&atmos_fluxes, stress_xw, stress_yw);
        fill_struct<double>(&p_as, fu10);
        fill_struct<double>(&p_sea_ice, concsum);
        m_is_struct_init = true;
    }

    m_impl->backend_tke_dp->calc(p_patch, p_cvmix, ocean_state, atmos_fluxes, p_as, p_sea_ice,
                          edges_block_size, edges_start_block, edges_end_block,
                          edges_start_index, edges_end_index, cells_block_size,
                          cells_start_block, cells_end_block, cells_start_index,
                          cells_end_index);
}

void YAOP::calc_tke(float *depth_CellInterface, float *prism_center_dist_c,
                    float *inv_prism_center_dist_c, float *prism_thick_c,
                    int *dolic_c, int *dolic_e, float *zlev_i, float *wet_c,
                    int *edges_cell_idx, int *edges_cell_blk,
                    float *temp, float *salt, float *stretch_c, float *eta_c,
                    float *p_vn_x1, float *p_vn_x2, float *p_vn_x3,
                    float *tke, float *tke_plc_in, float *hlc_in, float *wlc_in,
                    float *u_stokes_in, float *a_veloc_v, float *a_temp_v, float *a_salt_v,
                    float *iwe_Tdis, float *cvmix_dummy_1, float *cvmix_dummy_2,
                    float *cvmix_dummy_3, float *tke_Tbpr, float *tke_Tspr,
                    float *tke_Tdif, float *tke_Tdis, float *tke_Twin,
                    float *tke_Tiwf, float *tke_Tbck, float *tke_Ttot,
                    float *tke_Lmix, float *tke_Pr, float *stress_xw,
                    float *stress_yw, float *fu10, float *concsum,
                    int edges_block_size, int edges_start_block, int edges_end_block,
                    int edges_start_index, int edges_end_index, int cells_block_size,
                    int cells_start_block, int cells_end_block, int cells_start_index,
                    int cells_end_index) {
    // This structures are filled during the first call to initialize the memory views.
    // Then they are left empty during the following calls, but passed anyway to the backend.
    // It allows to have the struct templated as well as the backend without having to template
    // the YAOP class.
    // There can be more elegant workarounds to achieve that but this is the simplest one.
    struct t_patch<float> p_patch;
    struct t_cvmix<float> p_cvmix;
    struct t_sea_ice<float> p_sea_ice;
    struct t_atmos_for_ocean<float> p_as;
    struct t_atmo_fluxes<float> atmos_fluxes;
    struct t_ocean_state<float> ocean_state;

    if (!m_is_struct_init) {
        fill_struct<float>(&p_patch, depth_CellInterface, prism_center_dist_c,
                    inv_prism_center_dist_c, prism_thick_c, dolic_c, dolic_e,
                    zlev_i, wet_c, edges_cell_idx, edges_cell_blk);
        fill_struct<float>(&p_cvmix, tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v,
                    a_temp_v, a_salt_v, iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2,
                    cvmix_dummy_3, tke_Tbpr, tke_Tspr, tke_Tdif, tke_Tdis, tke_Twin,
                    tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr);
        fill_struct<float>(&ocean_state, temp, salt, stretch_c, eta_c, p_vn_x1, p_vn_x2, p_vn_x3);
        fill_struct<float>(&atmos_fluxes, stress_xw, stress_yw);
        fill_struct<float>(&p_as, fu10);
        fill_struct<float>(&p_sea_ice, concsum);
        m_is_struct_init = true;
    }

    m_impl->backend_tke_sp->calc(p_patch, p_cvmix, ocean_state, atmos_fluxes, p_as, p_sea_ice,
                          edges_block_size, edges_start_block, edges_end_block,
                          edges_start_index, edges_end_index, cells_block_size,
                          cells_start_block, cells_end_block, cells_start_index,
                          cells_end_index);
}

void YAOP::calc_vertical_stability() {}

void YAOP::calc_pp() {}

void YAOP::calc_idemix() {}
