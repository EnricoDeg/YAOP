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
#include "src/TKE_cuda.hpp"
#include "src/utils.hpp"
#include "src/cuda_memory.hpp"
#include "src/constants_thermodyn.hpp"

struct t_constant {
    int vert_mix_type;
    int vmix_idemix_tke;
    int vert_cor_type;
    double dtime;
    double OceanReferenceDensity;
    double grav;
    int l_lc;
    double clc;
    double ReferencePressureIndbars;
    double pi;
    int nlevs;
};

struct t_constant_tke {
    double c_k;
    double c_eps;
    double cd;
    double alpha_tke;
    double clc;
    double mxl_min;
    double KappaM_min;
    double KappaH_min;
    double KappaM_max;
    double tke_surf_min;
    double tke_min;
    int tke_mxl_choice;
    int handle_old_vals;
    bool only_tke;
    bool use_Kappa_min;
    bool use_ubound_dirichlet;
    bool use_lbound_dirichlet;
};

// Structures with memory views
struct t_cvmix_view p_cvmix_view_l;
struct t_patch_view p_patch_view_l;
struct t_ocean_state_view ocean_state_view_l;
struct t_atmo_fluxes_view atmos_fluxes_view_l;
struct t_atmos_for_ocean_view p_as_view_l;
struct t_sea_ice_view p_sea_ice_view_l;
struct t_tke_internal_view p_internal_view_l;

// Structures with parameters
struct t_constant p_constant;
struct t_constant_tke p_constant_tke;

// TKE CUDA kernels functions
__global__ void calc_impl_kernel(int blockNo, int start_index, int end_index,
                                 t_patch_view p_patch, t_cvmix_view p_cvmix,
                                 t_ocean_state_view ocean_state, t_atmo_fluxes_view atmos_fluxes,
                                 t_atmos_for_ocean_view p_as, t_sea_ice_view p_sea_ice,
                                 t_tke_internal_view p_internal, t_constant p_constant,
                                 t_constant_tke p_constant_tke);

__device__
void integrate(int jc, int nlevels, int blockNo, t_patch_view p_patch, t_cvmix_view p_cvmix,
               t_tke_internal_view p_internal, t_constant p_constant,
               t_constant_tke p_constant_tke);

__device__
void solve_tridiag(int jc, int nlevels, int blockNo, mdspan_2d_double a,
                   mdspan_2d_double b, mdspan_2d_double c, mdspan_2d_double d,
                   mdspan_3d_double x, mdspan_2d_double cp, mdspan_2d_double dp);

__device__
double  calculate_density(double temp, double salt, double pressure);

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
    p_internal_view_l.a_dif = view_cuda_malloc(m_a_dif, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.b_dif = view_cuda_malloc(m_b_dif, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.c_dif = view_cuda_malloc(m_c_dif, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.a_tri = view_cuda_malloc(m_a_tri, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.b_tri = view_cuda_malloc(m_b_tri, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.c_tri = view_cuda_malloc(m_c_tri, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.d_tri = view_cuda_malloc(m_d_tri, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.mxl = view_cuda_malloc(m_mxl, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.sqrttke = view_cuda_malloc(m_sqrttke, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.KappaM_out = view_cuda_malloc(m_KappaM_out,
                                                    static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.Rinum = view_cuda_malloc(m_Rinum, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.prandtl = view_cuda_malloc(m_prandtl, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.KappaH_out = view_cuda_malloc(m_KappaH_out,
                                                    static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.forc = view_cuda_malloc(m_forc, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.K_diss_v = view_cuda_malloc(m_K_diss_v,
                                                  static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.P_diss_v = view_cuda_malloc(m_P_diss_v,
                                                  static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.ke = view_cuda_malloc(m_ke, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.cp = view_cuda_malloc(m_cp, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    p_internal_view_l.dp = view_cuda_malloc(m_dp, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    is_view_init = false;
}

TKE_cuda::~TKE_cuda() {
    // Free internal arrays memory
    std::cout << "Finalizing TKE cuda... " << std::endl;
    check(cudaFree(m_rho_up));
    check(cudaFree(m_rho_down));
    check(cudaFree(m_forc_tke_surf_2D));
    check(cudaFree(m_forc_rho_surf_2D));
    check(cudaFree(m_bottom_fric_2D));
    check(cudaFree(m_s_c));
    check(cudaFree(m_dzw_stretched));
    check(cudaFree(m_dzt_stretched));
    check(cudaFree(m_tke_old));
    check(cudaFree(m_tke_Av));
    check(cudaFree(m_tke_kv));
    check(cudaFree(m_tke_iw_alpha_c));
    check(cudaFree(m_tke_iwe));
    check(cudaFree(m_tke_iwe_forcing));
    check(cudaFree(m_pressure));
    check(cudaFree(m_Nsqr));
    check(cudaFree(m_Ssqr));
    check(cudaFree(m_a_dif));
    check(cudaFree(m_b_dif));
    check(cudaFree(m_c_dif));
    check(cudaFree(m_a_tri));
    check(cudaFree(m_b_tri));
    check(cudaFree(m_c_tri));
    check(cudaFree(m_d_tri));
    check(cudaFree(m_mxl));
    check(cudaFree(m_sqrttke));
    check(cudaFree(m_KappaM_out));
    check(cudaFree(m_Rinum));
    check(cudaFree(m_prandtl));
    check(cudaFree(m_KappaH_out));
    check(cudaFree(m_forc));
    check(cudaFree(m_ke));
    check(cudaFree(m_cp));
    check(cudaFree(m_dp));
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

        p_constant.vert_mix_type = m_vert_mix_type;
        p_constant.vmix_idemix_tke = m_vmix_idemix_tke;
        p_constant.vert_cor_type = m_vert_cor_type;
        p_constant.dtime = m_dtime;
        p_constant.OceanReferenceDensity = m_OceanReferenceDensity;
        p_constant.grav = m_grav;
        p_constant.l_lc = m_l_lc;
        p_constant.clc = m_clc;
        p_constant.ReferencePressureIndbars = m_ReferencePressureIndbars;
        p_constant.pi = m_pi;
        p_constant.nlevs = m_nlevs;

        p_constant_tke.c_k = 0.1;
        p_constant_tke.c_eps = 0.7;
        p_constant_tke.cd = 3.75;
        p_constant_tke.alpha_tke = 30.0;
        p_constant_tke.clc = 0.15;
        p_constant_tke.mxl_min = 1.0e-8;
        p_constant_tke.KappaM_min = 1.0e-4;
        p_constant_tke.KappaH_min = 1.0e-5;
        p_constant_tke.KappaM_max = 100.0;
        p_constant_tke.tke_surf_min = 1.0e-4;
        p_constant_tke.tke_min = 1.0e-6;
        p_constant_tke.tke_mxl_choice = 2;
        p_constant_tke.handle_old_vals = 1;
        p_constant_tke.only_tke = true;
        p_constant_tke.use_Kappa_min = false;
        p_constant_tke.use_ubound_dirichlet = false;
        p_constant_tke.use_lbound_dirichlet = false;

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
                                                            p_internal_view_l, p_constant,
                                                            p_constant_tke);
    }
    check(cudaDeviceSynchronize());
}

__global__
void calc_impl_kernel(int blockNo, int start_index, int end_index, t_patch_view p_patch,
                      t_cvmix_view p_cvmix, t_ocean_state_view ocean_state,
                      t_atmo_fluxes_view atmos_fluxes,
                      t_atmos_for_ocean_view p_as, t_sea_ice_view p_sea_ice,
                      t_tke_internal_view p_internal, t_constant p_constant,
                      t_constant_tke p_constant_tke) {
    int jc = blockIdx.x * blockDim.x + threadIdx.x + start_index;
    if (jc <= end_index) {
        int levels = p_constant.nlevs;
        // initialization
        for (int level = 0; level < levels+1; level++) {
            p_internal.tke_kv(level, jc) = 0.0;
            p_internal.tke_Av(level, jc) = 0.0;
            p_internal.tke_iw_alpha_c(level, jc) = 0.0;
            p_internal.tke_iwe(level, jc) = 0.0;
            if (p_constant.vert_mix_type == p_constant.vmix_idemix_tke) {
                p_internal.tke_iwe_forcing(level, jc) = -1.0 * p_cvmix.iwe_Tdis(blockNo, level, jc);
            } else {
                p_internal.tke_iwe_forcing(level, jc) = 0.0;
            }
        }
        p_internal.forc_rho_surf_2D(jc) = 0.0;
        p_internal.bottom_fric_2D(jc) = 0.0;
        if (p_constant.vert_cor_type == 1)
            p_internal.s_c(jc) = ocean_state.stretch_c(blockNo, jc);
        else
            p_internal.s_c(jc) = 1.0;
        p_internal.pressure(0, jc) = 1.0;
        for (int level = 1; level < levels; level++)
            p_internal.pressure(level, jc) = p_patch.zlev_i(level) * p_constant.ReferencePressureIndbars;
        for (int level = 0; level < levels; level++)
            p_internal.dzw_stretched(level, jc) = p_patch.prism_thick_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);

        for (int level = 0; level < levels + 1; level++)
            p_internal.dzt_stretched(level, jc) = p_patch.prism_center_dist_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);

        // pre-integration
        levels = p_patch.dolic_c(blockNo, jc);
        if (levels > 0) {
            for (int level = 0; level < levels+1; level++) {
                p_internal.tke_old(level, jc) = p_cvmix.tke(blockNo, level, jc);
                p_internal.Nsqr(level, jc) = 0.0;
                p_internal.Ssqr(level, jc) = 0.0;
            }
            for (int level = 0; level < levels; level++) {
                p_internal.rho_up(level, jc) = 0.0;
                p_internal.rho_down(level, jc) = 0.0;
            }
            double tau_abs = (1.0 - p_sea_ice.concsum(blockNo, jc))
                             * sqrt(pow(atmos_fluxes.stress_xw(blockNo, jc), 2.0)
                                  + pow(atmos_fluxes.stress_yw(blockNo, jc), 2.0));
            p_internal.forc_tke_surf_2D(jc) = tau_abs / p_constant.OceanReferenceDensity;

            for (int level = 0; level < levels-1; level++)
                p_internal.rho_up(level, jc) = calculate_density(ocean_state.temp(blockNo, level, jc),
                                                                 ocean_state.salt(blockNo, level, jc),
                                                                 p_internal.pressure(level+1, jc));
            for (int level = 1; level < levels; level++)
                p_internal.rho_down(level, jc) = calculate_density(ocean_state.temp(blockNo, level, jc),
                                                                   ocean_state.salt(blockNo, level, jc),
                                                                   p_internal.pressure(level, jc));
            for (int level = 1; level < levels; level++)
                p_internal.Nsqr(level, jc) = p_constant.grav / p_constant.OceanReferenceDensity *
                                             (p_internal.rho_down(level, jc) - p_internal.rho_up(level-1, jc));


            for (int level = 1; level < levels; level++)
                p_internal.Ssqr(level, jc) = pow((ocean_state.p_vn_x1(blockNo, level-1, jc) -
                                             ocean_state.p_vn_x1(blockNo, level, jc) ) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc) , 2.0) +
                                             pow((ocean_state.p_vn_x2(blockNo, level-1, jc) -
                                             ocean_state.p_vn_x2(blockNo, level, jc) ) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc) , 2.0) +
                                             pow((ocean_state.p_vn_x3(blockNo, level-1, jc) -
                                             ocean_state.p_vn_x3(blockNo, level, jc) ) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc) , 2.0);

            // integration
            integrate(jc, p_constant.nlevs, blockNo, p_patch, p_cvmix, p_internal, p_constant, p_constant_tke);
        }

        // post-integration
    }
}

__device__
void integrate(int jc, int nlevels, int blockNo, t_patch_view p_patch, t_cvmix_view p_cvmix,
               t_tke_internal_view p_internal, t_constant p_constant,
               t_constant_tke p_constant_tke) {
    // Initialize diagnostics
    for (int level = 0; level < nlevels+1; level++) {
        p_cvmix.tke_Tbpr(blockNo, level, jc) = 0.0;
        p_cvmix.tke_Tspr(blockNo, level, jc) = 0.0;
        p_cvmix.tke_Tdif(blockNo, level, jc) = 0.0;
        p_cvmix.tke_Tdis(blockNo, level, jc) = 0.0;
        p_cvmix.tke_Twin(blockNo, level, jc) = 0.0;
        p_cvmix.tke_Tiwf(blockNo, level, jc) = 0.0;
        p_cvmix.tke_Tbck(blockNo, level, jc) = 0.0;
        p_cvmix.tke_Ttot(blockNo, level, jc) = 0.0;
        p_internal.a_dif(level, jc) = 0.0;
        p_internal.b_dif(level, jc) = 0.0;
        p_internal.c_dif(level, jc) = 0.0;
        p_internal.a_tri(level, jc) = 0.0;
        p_internal.b_tri(level, jc) = 0.0;
        p_internal.c_tri(level, jc) = 0.0;
    }

    // calculate mixing length scale
    for (int level = 0; level < nlevels+1; level++) {
        p_internal.sqrttke(level, jc) = sqrt(max(0.0, p_internal.tke_old(level, jc)));
        p_internal.mxl(level, jc) = sqrt(2.0) * p_internal.sqrttke(level, jc) /
                                    sqrt(max(1.0e-12, p_internal.Nsqr(level, jc)));
    }

    if (p_constant_tke.tke_mxl_choice == 2) {
        p_internal.mxl(0, jc) = 0.0;
        p_internal.mxl(nlevels, jc) = 0.0;
        for (int level = 1; level < nlevels; level++)
            p_internal.mxl(level, jc) = min(p_internal.mxl(level, jc),
                                        p_internal.mxl(level-1, jc) + p_internal.dzw_stretched(level-1, jc));
        p_internal.mxl(nlevels-1, jc) = min(p_internal.mxl(nlevels-1, jc),
                                        p_constant_tke.mxl_min +  p_internal.dzw_stretched(nlevels-1, jc));
        for (int level = nlevels-2; level > 0; level--)
            p_internal.mxl(level, jc) = min(p_internal.mxl(level, jc),
                                            p_internal.mxl(level+1, jc) +  p_internal.dzw_stretched(level, jc));
        for (int level = 0; level < nlevels+1; level++)
            p_internal.mxl(level, jc) = max(p_internal.mxl(level, jc), p_constant_tke.mxl_min);
    } else if (p_constant_tke.tke_mxl_choice == 3) {
        // TODO(EnricoDeg): not default
    } else {
        // Error
    }

    // calculate diffusivities
    for (int level = 0; level < nlevels+1; level++) {
        p_internal.KappaM_out(level, jc) = min(p_constant_tke.KappaM_max,
                                               p_constant_tke.c_k * p_internal.mxl(level, jc) *
                                               p_internal.sqrttke(level, jc));
        if (!p_constant_tke.only_tke)
            p_internal.Rinum(level, jc) = min(p_internal.Rinum(level, jc),
                                              p_internal.KappaM_out(level, jc) * p_internal.Nsqr(level, jc) /
                                              max(1.0e-12, p_internal.tke_iw_alpha_c(level, jc) *
                                              pow(p_internal.tke_iwe(level, jc), 2.0)));
        p_internal.prandtl(level, jc) = max(1.0, min(10.0, 6.6 * p_internal.Rinum(level, jc)));
        p_internal.KappaH_out(level, jc) = p_internal.KappaM_out(level, jc) / p_internal.prandtl(level, jc);
        if (p_constant_tke.use_Kappa_min) {
            p_internal.KappaM_out(level, jc) = max(p_constant_tke.KappaM_min, p_internal.KappaM_out(level, jc));
            p_internal.KappaH_out(level, jc) = max(p_constant_tke.KappaH_min, p_internal.KappaH_out(level, jc));
        }
    }

    // tke forcing
    // forcing by shear and buoycancy production
    p_internal.K_diss_v(0, jc) = p_internal.Ssqr(0, jc) * p_internal.KappaM_out(0, jc);
    p_internal.P_diss_v(0, jc) = p_internal.forc_rho_surf_2D(jc) * p_constant.grav * p_constant.OceanReferenceDensity;
    for (int level = 1; level < nlevels+1; level++) {
        p_internal.K_diss_v(level, jc) = p_internal.Ssqr(level, jc) * p_internal.KappaM_out(level, jc);
        p_internal.P_diss_v(level, jc) = p_internal.Nsqr(level, jc) * p_internal.KappaH_out(level, jc);
    }

    for (int level = 0; level < nlevels+1; level++) {
        p_internal.forc(level, jc) = p_internal.K_diss_v(level, jc) - p_internal.P_diss_v(level, jc);
        // additional langmuir turbulence term
        if (p_constant.l_lc)
            p_internal.forc(level, jc) += p_cvmix.tke_plc(blockNo, level, jc);
        // forcing by internal wave dissipation
        if (!p_constant_tke.only_tke)
            p_internal.forc(level, jc) += p_internal.tke_iwe_forcing(level, jc);
    }

    // vertical diffusion and dissipation is solved implicitely




    // solve the tri-diag matrix
//    solve_tridiag(jc, nlevels, blockNo, p_internal.a_tri, p_internal.b_tri, p_internal.c_tri,
//                  p_internal.d_tri, p_cvmix.tke, p_internal.cp, p_internal.dp);
}

__device__
void solve_tridiag(int jc, int nlevels, int blockNo, mdspan_2d_double a,
                   mdspan_2d_double b, mdspan_2d_double c, mdspan_2d_double d,
                   mdspan_3d_double x, mdspan_2d_double cp, mdspan_2d_double dp) {
    // initialize c-prime and d-prime
    cp(0, jc) = c(0, jc) / b(0, jc);
    dp(0, jc) = d(0, jc) / b(0, jc);
    // solve for vectors c-prime and d-prime
    for (int level = 1; level < nlevels+1; level++) {
        double m = b(level, jc) - cp(level-1, jc) * a(level, jc);
        double fxa = 1.0 / m;
        cp(level, jc) = c(level, jc) * fxa;
        dp(level, jc) = (d(level, jc) - dp(level-1, jc) * a(level, jc))*fxa;
    }
    // initialize solution x
    x(blockNo, nlevels, jc) = dp(nlevels, jc);
    // solve for x from the vectors c-prime and d-prime
    for (int level = nlevels-1; level >=0; level--)
        x(blockNo, level, jc) = dp(level, jc) - cp(level, jc) * x(blockNo, level+1, jc);
}

__device__
double  calculate_density(double temp, double salt, double pressure) {
    double dvs, fne, fst, qn3;
    double qnq, qvs, s, s3h;
    double t, denom, s__2;
    double rho;

    // This is the adisit part, that transforms potential in in-situ temperature
    qnq = -pressure * (-a_a3 + pressure * a_c3);
    qn3 = -pressure * a_a4;
    qvs = (pressure * (a_b1 - a_d * pressure)) *
          (salt - z_sref) + pressure * (a_a1 + pressure * (a_c1 - a_e1 * pressure));

    dvs = (a_b2 * pressure) * (salt - z_sref) +
           1.0 + pressure * (-a_a2 + pressure * (a_c2 - a_e2 * pressure));

    t   = (temp + qvs) / dvs;
    fne = - qvs + t * (dvs + t * (qnq + t * qn3)) - temp;

    fst = dvs + t * (2.0 * qnq + 3.0 * qn3 * t);

    t    = t - fne / fst;
    s    = max(salt, 0.0);
    s__2 = pow(s, 2.0);
    s3h  = s * sqrt(s);

    rho = r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))
        + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))
        + r_d0 * s__2 + s3h * (r_c0 + t * (r_c1 + r_c2 * t));

    denom = 1.0 - pressure / (pressure * (r_h0 + t *
            (r_h1 + t * (r_h2 + t * r_h3))
            + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))
            + r_aj0 * s3h + (r_ak0 + t * (r_ak1 + t * r_ak2)
            + s * (r_am0 + t * (r_am1 + t * r_am2))) * pressure)
            + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))
            + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))
            + s3h * (r_g0 + t * (r_g1 + r_g2 * t)));

    return rho/denom;
}
