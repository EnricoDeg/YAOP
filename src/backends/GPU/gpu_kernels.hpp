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

#ifndef SRC_BACKENDS_GPU_GPU_KERNELS_HPP_
#define SRC_BACKENDS_GPU_GPU_KERNELS_HPP_

#include "src/backends/GPU/CUDA/cuda_backend.hpp"
#include "src/shared/interface/memview_struct.hpp"


/*! \brief TKE computation on the cells.
*
*/
__global__ void calc_impl_cells(int blockNo, int start_index, int end_index,
                                 t_patch_view<gpu_memview::mdspan, gpu_memview::dextents> p_patch,
                                 t_cvmix_view<gpu_memview::mdspan, gpu_memview::dextents> p_cvmix,
                                 t_ocean_state_view<gpu_memview::mdspan, gpu_memview::dextents> ocean_state,
                                 t_atmo_fluxes_view<gpu_memview::mdspan, gpu_memview::dextents> atmos_fluxes,
                                 t_atmos_for_ocean_view<gpu_memview::mdspan, gpu_memview::dextents> p_as,
                                 t_sea_ice_view<gpu_memview::mdspan, gpu_memview::dextents> p_sea_ice,
                                 t_tke_internal_view<gpu_memview::mdspan, gpu_memview::dextents> p_internal,
                                 t_constant p_constant,
                                 t_constant_tke p_constant_tke);

/*! \brief TKE computation on the edges.
*
*/
__global__
void calc_impl_edges(int blockNo, int start_index, int end_index,
                     t_patch_view<gpu_memview::mdspan, gpu_memview::dextents> p_patch,
                     t_cvmix_view<gpu_memview::mdspan, gpu_memview::dextents> p_cvmix,
                     t_tke_internal_view<gpu_memview::mdspan, gpu_memview::dextents> p_internal,
                     t_constant p_constant);

/*! \brief TKE integration for each vertical level.
*
*/
__device__
void integrate(int jc, int blockNo,
               t_patch_view<gpu_memview::mdspan, gpu_memview::dextents> p_patch,
               t_cvmix_view<gpu_memview::mdspan, gpu_memview::dextents> p_cvmix,
               t_tke_internal_view<gpu_memview::mdspan, gpu_memview::dextents> p_internal,
               t_constant p_constant,
               t_constant_tke p_constant_tke);

/*! \brief Solve tridiagonal system.
*
*/
__device__
void solve_tridiag(int jc, int nlevels, int blockNo, mdspan_2d_double a,
                   mdspan_2d_double b, mdspan_2d_double c, mdspan_2d_double d,
                   mdspan_3d_double x, mdspan_2d_double cp, mdspan_2d_double dp);

/*! \brief Compute pointwise density from temperature, salinity and pressure.
*
*/
__device__
double  calculate_density(double temp, double salt, double pressure);

#endif  // SRC_BACKENDS_GPU_GPU_KERNELS_HPP_
