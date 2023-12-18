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

#include "src/backends/CPU/cpu_memory.hpp"
#include "src/shared/interface/memview_struct.hpp"

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

inline void calculate_mxl_2(int blockNo, int start_index, int end_index, int max_levels, double mxl_min,
                            mdspan_2d_int dolic_c, mdspan_3d_double tke_Lmix, mdspan_2d_double dzw_stretched);

inline void calculate_diffusivity(int blockNo, int start_index, int end_index, int max_levels,
                                  t_constant_tke *p_constant_tke,
                                  mdspan_2d_int dolic_c, mdspan_3d_double tke_Lmix, mdspan_2d_double sqrttke,
                                  mdspan_2d_double Nsqr, mdspan_2d_double Ssqr,
                                  mdspan_3d_double tke_Av, mdspan_2d_double tke_kv, mdspan_3d_double tke_Pr);

inline void forcing(int blockNo, int start_index, int end_index, int max_levels, bool l_lc, bool only_tke,
                    mdspan_2d_int dolic_c, mdspan_2d_double Ssqr, mdspan_2d_double Nsqr, mdspan_3d_double tke_Av,
                    mdspan_2d_double tke_kv, mdspan_3d_double tke_Tspr, mdspan_3d_double tke_Tbpr,
                    mdspan_3d_double tke_plc, mdspan_3d_double tke_Tiwf, mdspan_2d_double forc);

inline void solve_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d_int dolic_c,
                          mdspan_2d_double a, mdspan_2d_double b, mdspan_2d_double c, mdspan_2d_double d,
                          mdspan_3d_double x, mdspan_2d_double cp, mdspan_2d_double dp);

inline void vertical_diffusion(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d_int dolic_c,
                               double diff_surf_forc, double diff_bott_forc,
                               mdspan_2d_double a_dif, mdspan_2d_double b_dif, mdspan_2d_double c_dif,
                               mdspan_3d_double tke, mdspan_3d_double tke_Tdif);

inline void vertical_diffusion_ub_dirichlet(int blockNo, int start_index, int end_index, mdspan_2d_int dolic_c,
                                            double tke_surf, mdspan_2d_double ke, mdspan_2d_double dzw_stretched,
                                            mdspan_2d_double dzt_stretched, mdspan_3d_double tke,
                                            mdspan_3d_double tke_Tdif);

inline void vertical_diffusion_lb_dirichlet(int blockNo, int start_index, int end_index, mdspan_2d_int dolic_c,
                                            double tke_bott, mdspan_2d_double ke, mdspan_2d_double dzw_stretched,
                                            mdspan_2d_double dzt_stretched, mdspan_3d_double tke,
                                            mdspan_3d_double tke_Tdif);

#endif  // SRC_BACKENDS_CPU_CPU_KERNELS_HPP_
