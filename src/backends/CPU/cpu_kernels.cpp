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
#include <cmath>
#include "src/backends/CPU/cpu_kernels.hpp"
#include "src/backends/kernels.hpp"
#include "src/shared/constants/constants_thermodyn.hpp"

using std::max;
using std::min;

void calc_impl_cells(int blockNo, int start_index, int end_index,
                     t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
                     t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
                     t_ocean_state_view<cpu_memview::mdspan, cpu_memview::dextents> ocean_state,
                     t_atmo_fluxes_view<cpu_memview::mdspan, cpu_memview::dextents> atmos_fluxes,
                     t_atmos_for_ocean_view<cpu_memview::mdspan, cpu_memview::dextents> p_as,
                     t_sea_ice_view<cpu_memview::mdspan, cpu_memview::dextents> p_sea_ice,
                     t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
                     t_constant p_constant,
                     t_constant_tke p_constant_tke) {
    // initialization
    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_internal.tke_kv(level, jc) = 0.0;
            p_internal.tke_Av(blockNo, level, jc) = 0.0;
            if (p_constant.vert_mix_type == p_constant.vmix_idemix_tke) {
                p_cvmix.tke_Tiwf(blockNo, level, jc) = -1.0 * p_cvmix.iwe_Tdis(blockNo, level, jc);
            } else {
                p_cvmix.tke_Tiwf(blockNo, level, jc) = 0.0;
            }
            p_internal.dzt_stretched(level, jc) = p_patch.prism_center_dist_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);
        }
    }

    for (int level = 0; level < p_constant.nlevs; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.dzw_stretched(level, jc) = p_patch.prism_thick_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);

    // pre-integration
    for (int level = 0; level < p_constant.nlevs+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.tke_old(level, jc) = p_cvmix.tke(blockNo, level, jc);

    for (int jc = start_index; jc <= end_index; jc++) {
        double tau_abs = (1.0 - p_sea_ice.concsum(blockNo, jc))
                          * sqrt((atmos_fluxes.stress_xw(blockNo, jc) *
                                  atmos_fluxes.stress_xw(blockNo, jc))
                               + (atmos_fluxes.stress_yw(blockNo, jc) *
                                  atmos_fluxes.stress_yw(blockNo, jc)));
        p_internal.forc_tke_surf_2D(jc) = tau_abs / p_constant.OceanReferenceDensity;
    }

    // compute max level on block (maxval fortran function)
    int max_levels = 0;
    for (int jc = start_index; jc <= end_index; jc++)
        if (p_patch.dolic_c(blockNo, jc) > max_levels)
            max_levels = p_patch.dolic_c(blockNo, jc);

    // Loop over internal interfaces, surface (jk=1) and bottom (jk=kbot+1) excluded
    for (int level = 1; level < max_levels; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < p_patch.dolic_c(blockNo, jc)) {
                double rho_up = calculate_density(ocean_state.temp(blockNo, level-1, jc),
                                                  ocean_state.salt(blockNo, level-1, jc),
                                                  p_patch.zlev_i(level) *
                                                  p_constant.ReferencePressureIndbars);
                double rho_down = calculate_density(ocean_state.temp(blockNo, level, jc),
                                                    ocean_state.salt(blockNo, level, jc),
                                                    p_patch.zlev_i(level) *
                                                    p_constant.ReferencePressureIndbars);
                p_internal.Nsqr(level, jc) = p_constant.grav / p_constant.OceanReferenceDensity *
                                             (rho_down - rho_up) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc);
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
            }
        }
    }

    for (int level = 0; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.dzw_stretched(level, jc) = p_patch.prism_thick_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);
    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.dzt_stretched(level, jc) = p_patch.prism_center_dist_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);
}

void integrate(int blockNo, int start_index, int end_index,
               t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
               t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
               t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
               t_constant p_constant,
               t_constant_tke p_constant_tke) {
    double tke_surf, diff_surf_forc, tke_bott, diff_bott_forc;

    // compute max level on block (maxval fortran function)
    int max_levels = 0;
    for (int jc = start_index; jc <= end_index; jc++)
        if (p_patch.dolic_c(blockNo, jc) > max_levels)
            max_levels = p_patch.dolic_c(blockNo, jc);

    // Initialize diagnostics and calculate mixing length scale
    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_cvmix.tke_Twin(blockNo, level, jc) = 0.0;
            p_internal.sqrttke(level, jc) = sqrt(max(0.0, p_internal.tke_old(level, jc)));
            p_cvmix.tke_Lmix(blockNo, level, jc) = sqrt(2.0) * p_internal.sqrttke(level, jc) /
                                    sqrt(max(1.0e-12, p_internal.Nsqr(level, jc)));
        }
    }

    if (p_constant_tke.tke_mxl_choice == 2) {
        calculate_mxl_2(blockNo, start_index, end_index, max_levels, p_constant_tke.mxl_min,
                        p_patch.dolic_c, p_cvmix.tke_Lmix, p_internal.dzw_stretched);
    } else if (p_constant_tke.tke_mxl_choice == 3) {
        // TODO(EnricoDeg): not default
    } else {
        // Error
    }

    // calculate diffusivities
    calculate_diffusivity(blockNo, start_index, end_index, max_levels, &p_constant_tke,
                          p_patch.dolic_c, p_cvmix.tke_Lmix, p_internal.sqrttke,
                          p_internal.Nsqr, p_internal.Ssqr,
                          p_internal.tke_Av, p_internal.tke_kv, p_cvmix.tke_Pr);

    // tke forcing
    forcing(blockNo, start_index, end_index, max_levels, p_constant.l_lc, p_constant_tke.only_tke,
            p_patch.dolic_c, p_internal.Ssqr, p_internal.Nsqr, p_internal.tke_Av,
            p_internal.tke_kv, p_cvmix.tke_Tspr, p_cvmix.tke_Tbpr,
            p_cvmix.tke_plc, p_cvmix.tke_Tiwf, p_internal.forc);

    // vertical dissipation and diffusion solved implicitly
    // c is lower diagonal of matrix
    for (int level = 0; level < max_levels; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < p_patch.dolic_c(blockNo, jc)) {
                int dolic = p_patch.dolic_c(blockNo, jc);
                int kp1 = min(level+1, dolic-1);
                int kk = max(level, 1);
                p_internal.ke(level, jc) = 0.5 * p_constant_tke.alpha_tke *
                                       (p_internal.tke_Av(blockNo, kp1, jc) + p_internal.tke_Av(blockNo, kk, jc));
                p_internal.c_dif(level, jc) = p_internal.ke(level, jc) /
                              (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level, jc));
            }
        }
    }

    // not part of the diffusion matrix, thus value is arbitrary (set to zero)
    for (int jc = start_index; jc <= end_index; jc++)
        if (p_patch.dolic_c(blockNo, jc) >= 0)
            p_internal.c_dif(p_patch.dolic_c(blockNo, jc), jc) = 0.0;

    // b is main diagonal of matrix
    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < p_patch.dolic_c(blockNo, jc))
                p_internal.b_dif(level, jc) = p_internal.ke(level-1, jc) /
                                      (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level-1, jc)) +
                                      p_internal.ke(level, jc) /
                                      (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level, jc));
    // a is upper diagonal of matrix
    for (int level = 1; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < p_patch.dolic_c(blockNo, jc) + 1)
                p_internal.a_dif(level, jc) = p_internal.ke(level-1, jc) /
                                      (p_internal.dzt_stretched(level, jc) * p_internal.dzw_stretched(level-1, jc));

    // not part of the diffusion matrix, thus value is arbitrary (set to zero)
    for (int jc = start_index; jc <= end_index; jc++)
        if (p_patch.dolic_c(blockNo, jc) >= 0)
            p_internal.a_dif(0, jc) = 0.0;

    // copy tke_old
    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.tke_upd(level, jc) = p_internal.tke_old(level, jc);

    // upper boundary condition
    if (p_constant_tke.use_ubound_dirichlet) {
        for (int jc = start_index; jc <= end_index; jc++) {
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
        }
    } else {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_internal.forc(0, jc) += (p_constant_tke.cd * pow(p_internal.forc_tke_surf_2D(jc), 1.5)) /
                                      p_internal.dzt_stretched(0, jc);
            p_internal.b_dif(0, jc) = p_internal.ke(0, jc) /
                                      (p_internal.dzt_stretched(0, jc) * p_internal.dzw_stretched(0, jc));
            diff_surf_forc = 0.0;
        }
    }

    // lower boundary condition
    if (p_constant_tke.use_lbound_dirichlet) {
        for (int jc = start_index; jc <= end_index; jc++) {
            int dolic = p_patch.dolic_c(blockNo, jc);
            p_internal.sqrttke(dolic, jc) = 0.0;
            p_internal.forc(dolic, jc) = 0.0;
            tke_bott = p_constant_tke.tke_min;
            p_internal.tke_upd(dolic, jc) = tke_bott;
            diff_bott_forc = p_internal.c_dif(dolic-1, jc) * tke_bott;
            p_internal.forc(dolic-1, jc) += diff_bott_forc;
            p_internal.c_dif(dolic-1, jc) = 0.0;
            p_internal.b_dif(dolic, jc) = 0.0;
            p_internal.a_dif(dolic, jc) = 0.0;
        }
    } else {
        for (int jc = start_index; jc <= end_index; jc++) {
            int dolic = p_patch.dolic_c(blockNo, jc);
            p_internal.b_dif(dolic, jc) = p_internal.ke(dolic-1, jc) /
                                            (p_internal.dzt_stretched(dolic, jc) *
                                            p_internal.dzw_stretched(dolic-1, jc));
            diff_bott_forc = 0.0;
        }
    }

    // construct tridiagonal matrix to solve diffusion and dissipation implicitely
    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_internal.a_tri(level, jc) = - p_constant.dtime * p_internal.a_dif(level, jc);
            p_internal.b_tri(level, jc) = 1.0 + p_constant.dtime * p_internal.b_dif(level, jc);
            p_internal.c_tri(level, jc) = - p_constant.dtime * p_internal.c_dif(level, jc);
        }
    }

    for (int level = 0; level < max_levels+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < p_patch.dolic_c(blockNo, jc) + 1) {
                p_internal.d_tri(level, jc) = p_internal.tke_upd(level, jc) +
                                      p_constant.dtime * p_internal.forc(level, jc);
            }
        }
    }

    // solve the tri-diag matrix
    solve_tridiag(blockNo, start_index, end_index, max_levels, p_patch.dolic_c,
                  p_internal.a_tri, p_internal.b_tri, p_internal.c_tri,
                  p_internal.d_tri, p_cvmix.tke, p_internal.cp, p_internal.dp);
}

inline void calculate_mxl_2(int blockNo, int start_index, int end_index, int max_levels, double mxl_min,
                            mdspan_2d_int dolic_c, mdspan_3d_double tke_Lmix, mdspan_2d_double dzw_stretched) {
    for (int jc = start_index; jc <= end_index; jc++) {
        if (dolic_c(blockNo, jc) > 0) {
            tke_Lmix(blockNo, 0, jc) = 0.0;
            tke_Lmix(blockNo, dolic_c(blockNo, jc), jc) = 0.0;
        }
    }

    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                tke_Lmix(blockNo, level, jc) = min(tke_Lmix(blockNo, level, jc),
                         tke_Lmix(blockNo, level-1, jc) + dzw_stretched(level-1, jc));

    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc) > 0) {
            int dolic = dolic_c(blockNo, jc);
            tke_Lmix(blockNo, dolic-1, jc) = min(tke_Lmix(blockNo, dolic-1, jc),
                                             mxl_min + dzw_stretched(dolic-1, jc));
        }

    for (int level = max_levels-2; level > 0; level--)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) - 1)
                tke_Lmix(blockNo, level, jc) = min(tke_Lmix(blockNo, level, jc),
                         tke_Lmix(blockNo, level+1, jc) +  dzw_stretched(level, jc));

    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) + 1)
                tke_Lmix(blockNo, level, jc) = max(tke_Lmix(blockNo, level, jc), mxl_min);
}

inline void calculate_diffusivity(int blockNo, int start_index, int end_index, int max_levels,
                                  t_constant_tke *p_constant_tke,
                                  mdspan_2d_int dolic_c, mdspan_3d_double tke_Lmix, mdspan_2d_double sqrttke,
                                  mdspan_2d_double Nsqr, mdspan_2d_double Ssqr,
                                  mdspan_3d_double tke_Av, mdspan_2d_double tke_kv, mdspan_3d_double tke_Pr) {
    for (int level = 0; level < max_levels+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < dolic_c(blockNo, jc)) {
                tke_Av(blockNo, level, jc) = min(p_constant_tke->KappaM_max,
                                                 p_constant_tke->c_k * tke_Lmix(blockNo, level, jc) *
                                                 sqrttke(level, jc));
                tke_Pr(blockNo, level, jc) = Nsqr(level, jc) / max(Ssqr(level, jc), 1.0e-12);
                if (!p_constant_tke->only_tke)
                    tke_Pr(blockNo, level, jc) = min(tke_Pr(blockNo, level, jc),
                                                     tke_Av(blockNo, level, jc) * Nsqr(level, jc) / 1.0e-12);
                tke_Pr(blockNo, level, jc) = max(1.0, min(10.0, 6.6 * tke_Pr(blockNo, level, jc)));
                tke_kv(level, jc) = tke_Av(blockNo, level, jc) / tke_Pr(blockNo, level, jc);
                if (p_constant_tke->use_Kappa_min) {
                    tke_Av(blockNo, level, jc) = max(p_constant_tke->KappaM_min, tke_Av(blockNo, level, jc));
                    tke_kv(level, jc) = max(p_constant_tke->KappaH_min, tke_kv(level, jc));
                }
            }
        }
    }
}

inline void forcing(int blockNo, int start_index, int end_index, int max_levels, bool l_lc, bool only_tke,
                    mdspan_2d_int dolic_c, mdspan_2d_double Ssqr, mdspan_2d_double Nsqr, mdspan_3d_double tke_Av,
                    mdspan_2d_double tke_kv, mdspan_3d_double tke_Tspr, mdspan_3d_double tke_Tbpr,
                    mdspan_3d_double tke_plc, mdspan_3d_double tke_Tiwf, mdspan_2d_double forc) {
    for (int level = 0; level < max_levels+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < dolic_c(blockNo, jc) + 1) {
                // forcing by shear and buoycancy production
                tke_Tspr(blockNo, level, jc) = Ssqr(level, jc) * tke_Av(blockNo, level, jc);
                tke_Tbpr(blockNo, level, jc) = Nsqr(level, jc) * tke_kv(level, jc);
                if (level == 0) tke_Tbpr(blockNo, 0, jc) = 0.0;

                forc(level, jc) = tke_Tspr(blockNo, level, jc) - tke_Tbpr(blockNo, level, jc);
                // additional langmuir turbulence term
                if (l_lc)
                    forc(level, jc) += tke_plc(blockNo, level, jc);
                // forcing by internal wave dissipation
                if (!only_tke)
                    forc(level, jc) += tke_Tiwf(blockNo, level, jc);
            }
        }
    }
}

inline void solve_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d_int dolic_c,
                          mdspan_2d_double a, mdspan_2d_double b, mdspan_2d_double c, mdspan_2d_double d,
                          mdspan_3d_double x, mdspan_2d_double cp, mdspan_2d_double dp) {
    // initialize a-prime (cp) and d-prime
    for (int jc = start_index; jc <= end_index; jc++) {
        if (dolic_c(blockNo, jc)+1 > 0) {
            int dolic = dolic_c(blockNo, jc);
            cp(dolic, jc) = c(dolic, jc) / b(dolic, jc);
            dp(dolic, jc) = d(dolic, jc) / b(dolic, jc);
        }
    }

    // solve for vectors a-prime and d-prime
    for (int level = max_levels-1; level >= 0; level--) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level <= dolic_c(blockNo, jc)) {
                double fxa = 1.0 / (b(level, jc) - cp(level+1, jc) * c(level, jc));
            }
        }
    }

    // initialize x
    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc)+1 > 0)
            x(blockNo, 0, jc) = dp(0, jc);

    // solve for x from the vectors a-prime and d-prime
    for (int level = 1; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level <= dolic_c(blockNo, jc))
                x(blockNo, level, jc) = dp(level, jc) - cp(level, jc) * x(blockNo, level-1, jc);
}

void calc_impl_edges(int blockNo, int start_index, int end_index,
                     t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
                     t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
                     t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
                     t_constant p_constant) {
}
