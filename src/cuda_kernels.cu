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
#include "src/cuda_kernels.hpp"
#include "src/constants_thermodyn.hpp"

__global__
void calc_impl_cells(int blockNo, int start_index, int end_index, t_patch_view p_patch,
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
            p_internal.tke_Av(blockNo, level, jc) = 0.0;
            if (p_constant.vert_mix_type == p_constant.vmix_idemix_tke) {
                p_cvmix.tke_Tiwf(blockNo, level, jc) = -1.0 * p_cvmix.iwe_Tdis(blockNo, level, jc);
            } else {
                p_cvmix.tke_Tiwf(blockNo, level, jc) = 0.0;
            }
        }
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
            double tau_abs = (1.0 - p_sea_ice.concsum(blockNo, jc))
                             * sqrt(pow(atmos_fluxes.stress_xw(blockNo, jc), 2.0)
                                  + pow(atmos_fluxes.stress_yw(blockNo, jc), 2.0));
            p_internal.forc_tke_surf_2D(jc) = tau_abs / p_constant.OceanReferenceDensity;

            for (int level = 1; level < levels; level++) {
                double rho_down = calculate_density(ocean_state.temp(blockNo, level, jc),
                                                    ocean_state.salt(blockNo, level, jc),
                                                    p_internal.pressure(level, jc));
                double rho_up = calculate_density(ocean_state.temp(blockNo, level-1, jc),
                                                  ocean_state.salt(blockNo, level-1, jc),
                                                  p_internal.pressure(level, jc));
                p_internal.Nsqr(level, jc) = p_constant.grav / p_constant.OceanReferenceDensity *
                                             (rho_down - rho_up);
            }


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

            //  write tke vert. diffusivity to vert tracer diffusivities
            for (int level = 0; level < p_constant.nlevs+1; level++) {
                p_cvmix.a_temp_v(blockNo, level, jc) = p_internal.tke_kv(level, jc);
                p_cvmix.a_salt_v(blockNo, level, jc) = p_internal.tke_kv(level, jc);
            }
        }
    }
}

__device__
void integrate(int jc, int nlevels, int blockNo, t_patch_view p_patch, t_cvmix_view p_cvmix,
               t_tke_internal_view p_internal, t_constant p_constant,
               t_constant_tke p_constant_tke) {
    double tke_surf, diff_surf_forc, tke_bott, diff_bott_forc;

    // Initialize diagnostics
    for (int level = 0; level < nlevels+1; level++) {
        p_cvmix.tke_Twin(blockNo, level, jc) = 0.0;
    }

    // calculate mixing length scale
    for (int level = 0; level < nlevels+1; level++) {
        p_internal.sqrttke(level, jc) = sqrt(max(0.0, p_internal.tke_old(level, jc)));
        p_cvmix.tke_Lmix(blockNo, level, jc) = sqrt(2.0) * p_internal.sqrttke(level, jc) /
                                    sqrt(max(1.0e-12, p_internal.Nsqr(level, jc)));
    }

    if (p_constant_tke.tke_mxl_choice == 2) {
        p_cvmix.tke_Lmix(blockNo, 0, jc) = 0.0;
        p_cvmix.tke_Lmix(blockNo, nlevels, jc) = 0.0;
        for (int level = 1; level < nlevels; level++)
            p_cvmix.tke_Lmix(blockNo, level, jc) = min(p_cvmix.tke_Lmix(blockNo, level, jc),
                                        p_cvmix.tke_Lmix(blockNo, level-1, jc) + p_internal.dzw_stretched(level-1, jc));
        p_cvmix.tke_Lmix(blockNo, nlevels-1, jc) = min(p_cvmix.tke_Lmix(blockNo, nlevels-1, jc),
                                        p_constant_tke.mxl_min +  p_internal.dzw_stretched(nlevels-1, jc));
        for (int level = nlevels-2; level > 0; level--)
            p_cvmix.tke_Lmix(blockNo, level, jc) = min(p_cvmix.tke_Lmix(blockNo, level, jc),
                                        p_cvmix.tke_Lmix(blockNo, level+1, jc) +  p_internal.dzw_stretched(level, jc));
        for (int level = 0; level < nlevels+1; level++)
            p_cvmix.tke_Lmix(blockNo, level, jc) = max(p_cvmix.tke_Lmix(blockNo, level, jc), p_constant_tke.mxl_min);
    } else if (p_constant_tke.tke_mxl_choice == 3) {
        // TODO(EnricoDeg): not default
    } else {
        // Error
    }

    // calculate diffusivities
    for (int level = 0; level < nlevels+1; level++) {
        p_internal.tke_Av(blockNo, level, jc) = min(p_constant_tke.KappaM_max,
                                               p_constant_tke.c_k * p_cvmix.tke_Lmix(blockNo, level, jc) *
                                               p_internal.sqrttke(level, jc));
        if (!p_constant_tke.only_tke)
            p_internal.Rinum(level, jc) = min(p_internal.Rinum(level, jc),
                                              p_internal.tke_Av(blockNo, level, jc) * p_internal.Nsqr(level, jc) /
                                              1.0e-12);
        p_cvmix.tke_Pr(blockNo, level, jc) = max(1.0, min(10.0, 6.6 * p_internal.Rinum(level, jc)));
        p_internal.tke_kv(level, jc) = p_internal.tke_Av(blockNo, level, jc) / p_cvmix.tke_Pr(blockNo, level, jc);
        if (p_constant_tke.use_Kappa_min) {
            p_internal.tke_Av(blockNo, level, jc) = max(p_constant_tke.KappaM_min,
                                                        p_internal.tke_Av(blockNo, level, jc));
            p_internal.tke_kv(level, jc) = max(p_constant_tke.KappaH_min, p_internal.tke_kv(level, jc));
        }
    }

    // tke forcing
    // forcing by shear and buoycancy production
    p_cvmix.tke_Tspr(blockNo, 0, jc) = p_internal.Ssqr(0, jc) * p_internal.tke_Av(blockNo, 0, jc);
    p_cvmix.tke_Tbpr(blockNo, 0, jc) = 0.0;
    for (int level = 1; level < nlevels+1; level++) {
        p_cvmix.tke_Tspr(blockNo, level, jc) = p_internal.Ssqr(level, jc) * p_internal.tke_Av(blockNo, level, jc);
        p_cvmix.tke_Tbpr(blockNo, level, jc) = p_internal.Nsqr(level, jc) * p_internal.tke_kv(level, jc);
    }

    for (int level = 0; level < nlevels+1; level++) {
        p_internal.forc(level, jc) = p_cvmix.tke_Tspr(blockNo, level, jc) - p_cvmix.tke_Tbpr(blockNo, level, jc);
        // additional langmuir turbulence term
        if (p_constant.l_lc)
            p_internal.forc(level, jc) += p_cvmix.tke_plc(blockNo, level, jc);
        // forcing by internal wave dissipation
        if (!p_constant_tke.only_tke)
            p_internal.forc(level, jc) += p_cvmix.tke_Tiwf(blockNo, level, jc);
    }

    // vertical diffusion and dissipation is solved implicitely
    for (int level = 0; level < nlevels; level++) {
        int kp1 = min(level+1, nlevels-1);
        int kk = max(level, 1);
        p_internal.ke(level, jc) = 0.5 * p_constant_tke.alpha_tke *
                                   (p_internal.tke_Av(blockNo, kp1, jc) + p_internal.tke_Av(blockNo, kk, jc));
    }
    p_internal.ke(nlevels, jc) = 0.0;

    // a is upper diagonal of matrix
    // b is main diagonal of matrix
    // c is lower diagonal of matrix
    p_internal.a_dif(0, jc) = 0.0;
    p_internal.c_dif(0, jc) = p_internal.ke(0, jc) /
                              (p_internal.dzt_stretched(0, jc) * p_internal.dzw_stretched(0, jc));

    for (int level = 1; level < nlevels; level++) {
        p_internal.a_dif(level, jc) = p_internal.ke(level-1, jc) /
                                      (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level-1, jc));
        p_internal.b_dif(level, jc) = p_internal.ke(level-1, jc) /
                                      (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level-1, jc)) +
                                      p_internal.ke(level, jc) /
                                      (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level, jc));
        p_internal.c_dif(level, jc) = p_internal.ke(level, jc) /
                              (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level, jc));
    }

    p_internal.a_dif(nlevels, jc) = p_internal.ke(nlevels-1, jc) /
                                    (p_internal.dzt_stretched(nlevels, jc) * p_internal.dzw_stretched(nlevels-1, jc));
    p_internal.c_dif(nlevels, jc) = 0.0;

    // copy tke_old
    for (int level = 0; level < nlevels+1; level++)
        p_internal.tke_upd(level, jc) = p_internal.tke_old(level, jc);

    // upper boundary condition
    if (p_constant_tke.use_ubound_dirichlet) {
        p_internal.sqrttke(0, jc) = 0.0;
        p_internal.forc(0, jc) = 0.0;
        tke_surf = max(p_constant_tke.tke_surf_min,
                       p_constant_tke.cd * p_internal.forc_tke_surf_2D(jc));
        p_internal.tke_upd(0, jc) = tke_surf;
        diff_surf_forc = p_internal.a_dif(1, jc) * tke_surf;
        p_internal.forc(1, jc) += diff_surf_forc;
        p_internal.a_dif(1, jc) = 0.0;
        p_internal.b_dif(0, jc) = 0.0;
        p_internal.c_dif(0, jc) = 0.0;
    } else {
        p_internal.forc(0, jc) += (p_constant_tke.cd * pow(p_internal.forc_tke_surf_2D(jc), 1.5)) /
                                  p_internal.dzt_stretched(0, jc);
        p_internal.b_dif(0, jc) = p_internal.ke(0, jc) /
                                  (p_internal.dzt_stretched(0, jc) * p_internal.dzw_stretched(0, jc));
        diff_surf_forc = 0.0;
    }

    // lower boundary condition
    if (p_constant_tke.use_lbound_dirichlet) {
        p_internal.sqrttke(nlevels, jc) = 0.0;
        p_internal.forc(nlevels, jc) = 0.0;
        tke_bott = p_constant_tke.tke_min;
        p_internal.tke_upd(nlevels, jc) = tke_bott;
        diff_bott_forc = p_internal.c_dif(nlevels-1, jc) * tke_bott;
        p_internal.forc(nlevels-1, jc) += diff_bott_forc;
        p_internal.c_dif(nlevels-1, jc) = 0.0;
        p_internal.b_dif(nlevels, jc) = 0.0;
        p_internal.a_dif(nlevels, jc) = 0.0;
    } else {
        p_internal.b_dif(nlevels, jc) = p_internal.ke(nlevels-1, jc) /
                                        (p_internal.dzt_stretched(nlevels, jc) *
                                        p_internal.dzw_stretched(nlevels-1, jc));
        diff_bott_forc = 0.0;
    }

    // construct tridiagonal matrix to solve diffusion and dissipation implicitely
    for (int level = 0; level < nlevels+1; level++) {
        p_internal.a_tri(level, jc) = - p_constant.dtime * p_internal.a_dif(level, jc);
        p_internal.b_tri(level, jc) = 1.0 + p_constant.dtime * p_internal.b_dif(level, jc);
        p_internal.c_tri(level, jc) = - p_constant.dtime * p_internal.c_dif(level, jc);
    }

    for (int level = 1; level < nlevels; level++)
        p_internal.b_tri(level, jc) += p_constant.dtime * p_constant_tke.c_eps *
                                       p_internal.sqrttke(level, jc) / p_cvmix.tke_Lmix(blockNo, level, jc);

    // d is r.h.s. of implicite equation (d: new tke with only explicite tendencies included)
    for (int level = 0; level < nlevels+1; level++)
        p_internal.d_tri(level, jc) = p_internal.tke_upd(level, jc) +
                                      p_constant.dtime * p_internal.forc(level, jc);

    // solve the tri-diag matrix
    solve_tridiag(jc, nlevels, blockNo, p_internal.a_tri, p_internal.b_tri, p_internal.c_tri,
                  p_internal.d_tri, p_cvmix.tke, p_internal.cp, p_internal.dp);

    // diagnose implicit tendencies (only for diagnostics)
    // vertical diffusion of TKE
    for (int level = 1; level < nlevels; level++)
        p_cvmix.tke_Tdif(blockNo, level, jc) = p_internal.a_dif(level, jc) * p_cvmix.tke(blockNo, level-1, jc) -
                                               p_internal.b_dif(level, jc) * p_cvmix.tke(blockNo, level, jc) +
                                               p_internal.c_dif(level, jc) * p_cvmix.tke(blockNo, level+1, jc);

    p_cvmix.tke_Tdif(blockNo, 0, jc) = - p_internal.b_dif(0, jc) * p_cvmix.tke(blockNo, 0, jc) +
                                         p_internal.c_dif(0, jc) * p_cvmix.tke(blockNo, 1, jc);
    p_cvmix.tke_Tdif(blockNo, nlevels, jc) = p_internal.a_dif(nlevels, jc) * p_cvmix.tke(blockNo, nlevels-1, jc) -
                                             p_internal.b_dif(nlevels, jc) * p_cvmix.tke(blockNo, nlevels, jc);
    p_cvmix.tke_Tdif(blockNo, 1, jc) += diff_surf_forc;
    p_cvmix.tke_Tdif(blockNo, nlevels-1, jc) += diff_bott_forc;

    // flux out of first box due to diffusion with Dirichlet boundary value of TKE
    // (tke_surf=tke_upd(0)) and TKE of box below (tke_new(1))
    if (p_constant_tke.use_ubound_dirichlet)
        p_cvmix.tke_Tdif(blockNo, 0, jc) = - p_internal.ke(0, jc) / p_internal.dzw_stretched(0, jc) /
                                           p_internal.dzt_stretched(0, jc) *
                                           (tke_surf - p_cvmix.tke(blockNo, 1, jc));

    if (p_constant_tke.use_lbound_dirichlet)
        p_cvmix.tke_Tdif(blockNo, nlevels, jc) = p_internal.ke(nlevels-1, jc) /
                                                 p_internal.dzw_stretched(nlevels-1, jc) /
                                                 p_internal.dzt_stretched(nlevels, jc) *
                                                 (p_cvmix.tke(blockNo, nlevels-1, jc) - tke_bott);

    // dissipation of TKE
    p_cvmix.tke_Tdis(blockNo, 0, jc) = 0.0;
    p_cvmix.tke_Tdis(blockNo, nlevels, jc) = 0.0;
    for (int level = 1; level < nlevels; level++)
        p_cvmix.tke_Tdis(blockNo, level, jc) = - p_constant_tke.c_eps / p_cvmix.tke_Lmix(blockNo, level, jc) *
                                                 p_internal.sqrttke(level, jc) * p_cvmix.tke(blockNo, level, jc);

    // Part 5: reset tke to bounding values
    // copy of unrestored tke to diagnose energy input by restoring
    for (int level = 0; level < nlevels+1; level++)
        p_internal.tke_unrest(level, jc) = p_cvmix.tke(blockNo, level, jc);

    // restrict values of TKE to tke_min, if IDEMIX is not used
    if (p_constant_tke.only_tke) {
        for (int level = 0; level < nlevels+1; level++) {
            p_cvmix.tke(blockNo, level, jc) = max(p_cvmix.tke(blockNo, level, jc), p_constant_tke.tke_min);
        }
    }

    // Part 6: Assign diagnostic variables
    for (int level = 0; level < nlevels+1; level++) {
        p_cvmix.tke_Tbpr(blockNo, level, jc) *= -1.0;
        p_cvmix.tke_Tbck(blockNo, level, jc) = (p_cvmix.tke(blockNo, level, jc) -
                                                p_internal.tke_unrest(level, jc)) /
                                                p_constant.dtime;
    }

    if (p_constant_tke.use_ubound_dirichlet) {
        p_cvmix.tke_Twin(blockNo, 0, jc) = (p_cvmix.tke(blockNo, 0, jc) - p_internal.tke_old(0, jc)) /
                                           p_constant.dtime - p_cvmix.tke_Tdif(blockNo, 0, jc);
        p_cvmix.tke_Tbck(blockNo, 0, jc) = 0.0;
    } else {
        p_cvmix.tke_Twin(blockNo, 0, jc) = (p_constant_tke.cd * pow(p_internal.forc_tke_surf_2D(jc), 1.5)) /
                                           p_internal.dzt_stretched(0, jc);
    }

    if (p_constant_tke.use_lbound_dirichlet) {
        p_cvmix.tke_Twin(blockNo, nlevels, jc) = (p_cvmix.tke(blockNo, nlevels, jc) -
                                                  p_internal.tke_old(nlevels, jc)) /
                                                  p_constant.dtime -
                                                  p_cvmix.tke_Tdif(blockNo, nlevels, jc);
        p_cvmix.tke_Tbck(blockNo, nlevels, jc) = 0.0;
    } else {
        p_cvmix.tke_Twin(blockNo, nlevels, jc) = 0.0;
    }

    for (int level = 0; level < nlevels+1; level++) {
        p_cvmix.tke_Ttot(blockNo, level, jc) = (p_cvmix.tke(blockNo, level, jc) -
                                                p_internal.tke_old(level, jc)) / p_constant.dtime;
    }

    // the rest is for debugging
    for (int level = 0; level < nlevels+1; level++) {
        p_cvmix.cvmix_dummy_1(blockNo, level, jc) = p_internal.tke_kv(level, jc);
        p_cvmix.cvmix_dummy_2(blockNo, level, jc) = p_internal.tke_Av(blockNo, level, jc);
        p_cvmix.cvmix_dummy_3(blockNo, level, jc) = p_internal.Nsqr(level, jc);
    }
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

__global__
void calc_impl_edges(int blockNo, int start_index, int end_index, t_patch_view p_patch,
                      t_cvmix_view p_cvmix, t_tke_internal_view p_internal, t_constant p_constant) {
    int je = blockIdx.x * blockDim.x + threadIdx.x + start_index;
    if (je <= end_index) {
        int levels = p_constant.nlevs;
        for (int level = 0; level < levels+1; level++)
            p_cvmix.a_veloc_v(blockNo, level, je) = 0.0;

        levels = p_patch.dolic_e(blockNo, je);
        int cell_1_idx = p_patch.edges_cell_idx(0, blockNo, je);
        int cell_1_block = p_patch.edges_cell_blk(0, blockNo, je);
        int cell_2_idx = p_patch.edges_cell_idx(1, blockNo, je);
        int cell_2_block = p_patch.edges_cell_blk(1, blockNo, je);
        for (int level = 1; level < levels; level++)
            p_cvmix.a_veloc_v(blockNo, level, je) = 0.5 * (p_internal.tke_Av(cell_1_block, level, cell_1_idx) +
                                                  p_internal.tke_Av(cell_2_block, level, cell_2_idx));
    }
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
