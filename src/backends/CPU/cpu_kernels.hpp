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

#ifndef SRC_BACKENDS_CPU_CPU_KERNELS_HPP_
#define SRC_BACKENDS_CPU_CPU_KERNELS_HPP_

#include <algorithm>
#include <cmath>
#include "src/backends/CPU/cpu_memory.hpp"
#include "src/shared/interface/memview_struct.hpp"
#include "src/shared/interface/data_struct.hpp"

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
                     t_constant_tke p_constant_tke);

void calc_impl_edges(int blockNo, int start_index, int end_index,
                     t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
                     t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
                     t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
                     t_constant p_constant);

void integrate(int blockNo, int start_index, int end_index,
               t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
               t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
               t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
               t_constant p_constant,
               t_constant_tke p_constant_tke);

template <class T>
inline
void calc_mxl_2(int blockNo, int start_index, int end_index, int max_levels, T mxl_min,
                mdspan_2d<int> dolic_c, mdspan_3d<T> tke_Lmix, mdspan_2d<T> dzw_stretched) {
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

template <class T>
inline
void calc_diffusivity(int blockNo, int start_index, int end_index, int max_levels,
                      t_constant_tke *p_constant_tke,
                      mdspan_2d<int> dolic_c, mdspan_3d<T> tke_Lmix, mdspan_2d<T> sqrttke,
                      mdspan_2d<T> Nsqr, mdspan_2d<T> Ssqr,
                      mdspan_3d<T> tke_Av, mdspan_2d<T> tke_kv, mdspan_3d<T> tke_Pr) {
    for (int level = 0; level < max_levels+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < dolic_c(blockNo, jc) + 1) {
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

template <class T>
inline
void calc_forcing(int blockNo, int start_index, int end_index, int max_levels, bool l_lc, bool only_tke,
                  mdspan_2d<int> dolic_c, mdspan_2d<T> Ssqr, mdspan_2d<T> Nsqr, mdspan_3d<T> tke_Av,
                  mdspan_2d<T> tke_kv, mdspan_3d<T> tke_Tspr, mdspan_3d<T> tke_Tbpr,
                  mdspan_3d<T> tke_plc, mdspan_3d<T> tke_Tiwf, mdspan_2d<T> forc) {
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

template <class T>
inline
void build_diffusion_dissipation_tridiag(int blockNo, int start_index, int end_index, int max_levels,
                                         mdspan_2d<int> dolic_c, T alpha_tke,
                                         mdspan_3d<T> tke_Av, mdspan_2d<T> dzt_stretched,
                                         mdspan_2d<T> dzw_stretched,
                                         mdspan_2d<T> ke, mdspan_2d<T> a_dif, mdspan_2d<T> b_dif,
                                         mdspan_2d<T> c_dif) {
    // c is lower diagonal of matrix
    for (int level = 0; level < max_levels; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < dolic_c(blockNo, jc)) {
                int kp1 = min(level+1, dolic_c(blockNo, jc)-1);
                int kk = max(level, 1);
                ke(level, jc) = 0.5 * alpha_tke * (tke_Av(blockNo, kp1, jc) + tke_Av(blockNo, kk, jc));
                c_dif(level, jc) = ke(level, jc) / (dzt_stretched(level, jc) * dzw_stretched(level, jc));
            }
        }
    }

    // not part of the diffusion matrix, thus value is arbitrary (set to zero)
    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc) >= 0)
            c_dif(dolic_c(blockNo, jc), jc) = 0.0;

    // b is main diagonal of matrix
    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                b_dif(level, jc) = ke(level-1, jc) / (dzt_stretched(level, jc) * dzw_stretched(level-1, jc)) +
                                   ke(level, jc) / (dzt_stretched(level, jc) * dzw_stretched(level, jc));

    // a is upper diagonal of matrix
    for (int level = 1; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) + 1)
                a_dif(level, jc) = ke(level-1, jc) / (dzt_stretched(level, jc) * dzw_stretched(level-1, jc));

    // not part of the diffusion matrix, thus value is arbitrary (set to zero)
    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc) >= 0)
            a_dif(0, jc) = 0.0;
}

template <class T>
inline
void build_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                   T dtime, T c_eps, int nlevs,
                   mdspan_2d<T> a_dif, mdspan_2d<T> b_dif, mdspan_2d<T> c_dif,
                   mdspan_2d<T> sqrttke, mdspan_3d<T> tke_Lmix, mdspan_2d<T> tke_upd,
                   mdspan_2d<T> forc,
                   mdspan_2d<T> a_tri, mdspan_2d<T> b_tri, mdspan_2d<T> c_tri,
                   mdspan_2d<T> d_tri) {
    for (int level = 0; level < nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            a_tri(level, jc) = - dtime * a_dif(level, jc);
            b_tri(level, jc) = 1.0 + dtime * b_dif(level, jc);
            c_tri(level, jc) = - dtime * c_dif(level, jc);
        }
    }

    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                b_tri(level, jc) = b_tri(level, jc) + dtime * c_eps * sqrttke(level, jc) /
                                   tke_Lmix(blockNo, level, jc);

    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) + 1)
                d_tri(level, jc) = tke_upd(level, jc) + dtime * forc(level, jc);
}

template <class T>
inline
void solve_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                   mdspan_2d<T> a, mdspan_2d<T> b, mdspan_2d<T> c, mdspan_2d<T> d,
                   mdspan_3d<T> x, mdspan_2d<T> cp, mdspan_2d<T> dp) {
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
                T fxa = 1.0 / (b(level, jc) - cp(level+1, jc) * c(level, jc));
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

template <class T>
inline
void tke_vertical_diffusion(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                            T diff_surf_forc, T diff_bott_forc,
                            mdspan_2d<T> a_dif, mdspan_2d<T> b_dif, mdspan_2d<T> c_dif,
                            mdspan_3d<T> tke, mdspan_3d<T> tke_Tdif) {
    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                 tke_Tdif(blockNo, level, jc) = a_dif(level, jc) * tke(blockNo, level-1, jc) -
                                                b_dif(level, jc) * tke(blockNo, level, jc) +
                                                c_dif(level, jc) * tke(blockNo, level+1, jc);

    for (int jc = start_index; jc <= end_index; jc++) {
        if (dolic_c(blockNo, jc) > 0) {
            int dolic = dolic_c(blockNo, jc);
            tke_Tdif(blockNo, 0, jc) = - b_dif(0, jc) * tke(blockNo, 0, jc) +
                                         c_dif(0, jc) * tke(blockNo, 1, jc);
            tke_Tdif(blockNo, 1, jc) += diff_surf_forc;
            tke_Tdif(blockNo, dolic-1, jc) += diff_bott_forc;
            tke_Tdif(blockNo, dolic, jc) = a_dif(dolic, jc) * tke(blockNo, dolic-1, jc) -
                                           b_dif(dolic, jc) * tke(blockNo, dolic, jc);
        }
    }
}

inline
void tke_vertical_diffusion_ub_dirichlet(int blockNo, int start_index, int end_index, mdspan_2d<int> dolic_c,
                                         double tke_surf, mdspan_2d<double> ke, mdspan_2d<double> dzw_stretched,
                                         mdspan_2d<double> dzt_stretched, mdspan_3d<double> tke,
                                         mdspan_3d<double> tke_Tdif);

inline
void tke_vertical_diffusion_lb_dirichlet(int blockNo, int start_index, int end_index, mdspan_2d<int> dolic_c,
                                         double tke_bott, mdspan_2d<double> ke, mdspan_2d<double> dzw_stretched,
                                         mdspan_2d<double> dzt_stretched, mdspan_3d<double> tke,
                                         mdspan_3d<double> tke_Tdif);

inline
void tke_vertical_dissipation(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                              int nlevs, double c_eps, mdspan_3d<double> tke_Lmix, mdspan_2d<double> sqrttke,
                              mdspan_3d<double> tke, mdspan_3d<double> tke_Tdis);

#endif  // SRC_BACKENDS_CPU_CPU_KERNELS_HPP_
