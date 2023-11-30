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

#ifndef SRC_CUDA_KERNELS_HPP_
#define SRC_CUDA_KERNELS_HPP_

#include "src/cuda_memory.hpp"
#include "src/memview_struct.hpp"

__global__ void calc_impl_cells(int blockNo, int start_index, int end_index,
                                 t_patch_view p_patch,
                                 t_cvmix_view<mdspan_2d_double, mdspan_3d_double> p_cvmix,
                                 t_ocean_state_view ocean_state,
                                 t_atmo_fluxes_view<mdspan_2d_double> atmos_fluxes,
                                 t_atmos_for_ocean_view<mdspan_2d_double> p_as,
                                 t_sea_ice_view<mdspan_2d_double> p_sea_ice,
                                 t_tke_internal_view<mdspan_1d_double, mdspan_2d_double, mdspan_3d_double> p_internal,
                                 t_constant p_constant,
                                 t_constant_tke p_constant_tke);

__global__
void calc_impl_edges(int blockNo, int start_index, int end_index, t_patch_view p_patch,
                     t_cvmix_view<mdspan_2d_double, mdspan_3d_double> p_cvmix,
                     t_tke_internal_view<mdspan_1d_double, mdspan_2d_double, mdspan_3d_double> p_internal,
                     t_constant p_constant);

__device__
void integrate(int jc, int blockNo, t_patch_view p_patch,
               t_cvmix_view<mdspan_2d_double, mdspan_3d_double> p_cvmix,
               t_tke_internal_view<mdspan_1d_double, mdspan_2d_double, mdspan_3d_double> p_internal,
               t_constant p_constant,
               t_constant_tke p_constant_tke);

__device__
void solve_tridiag(int jc, int nlevels, int blockNo, mdspan_2d_double a,
                   mdspan_2d_double b, mdspan_2d_double c, mdspan_2d_double d,
                   mdspan_3d_double x, mdspan_2d_double cp, mdspan_2d_double dp);

__device__
double  calculate_density(double temp, double salt, double pressure);

#endif  // SRC_CUDA_KERNELS_HPP_
