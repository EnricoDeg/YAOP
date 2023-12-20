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

inline
void build_diffusion_dissipation_tridiag(int blockNo, int start_index, int end_index, int max_levels,
                                         mdspan_2d<int> dolic_c, double alpha_tke,
                                         mdspan_3d<double> tke_Av, mdspan_2d<double> dzt_stretched,
                                         mdspan_2d<double> dzw_stretched,
                                         mdspan_2d<double> ke, mdspan_2d<double> a_dif, mdspan_2d<double> b_dif,
                                         mdspan_2d<double> c_dif);

inline
void build_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                   double dtime, double c_eps, int nlevs,
                   mdspan_2d<double> a_dif, mdspan_2d<double> b_dif, mdspan_2d<double> c_dif,
                   mdspan_2d<double> sqrttke, mdspan_3d<double> tke_Lmix, mdspan_2d<double> tke_upd,
                   mdspan_2d<double> forc,
                   mdspan_2d<double> a_tri, mdspan_2d<double> b_tri, mdspan_2d<double> c_tri,
                   mdspan_2d<double> d_tri);

inline
void solve_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                   mdspan_2d<double> a, mdspan_2d<double> b, mdspan_2d<double> c, mdspan_2d<double> d,
                   mdspan_3d<double> x, mdspan_2d<double> cp, mdspan_2d<double> dp);

inline
void tke_vertical_diffusion(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                            double diff_surf_forc, double diff_bott_forc,
                            mdspan_2d<double> a_dif, mdspan_2d<double> b_dif, mdspan_2d<double> c_dif,
                            mdspan_3d<double> tke, mdspan_3d<double> tke_Tdif);

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
