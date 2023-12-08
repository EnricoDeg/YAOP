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

#include <iostream>
#include "src/TKE_cpu.hpp"
#include "src/cpu_memory.hpp"
#include "src/utils.hpp"

// Structures with memory views
static struct t_cvmix_view<Kokkos::mdspan, Kokkos::dextents> p_cvmix_view;
static struct t_patch_view<Kokkos::mdspan, Kokkos::dextents> p_patch_view;
static struct t_ocean_state_view<Kokkos::mdspan, Kokkos::dextents> ocean_state_view;
static struct t_atmo_fluxes_view<Kokkos::mdspan, Kokkos::dextents> atmos_fluxes_view;
static struct t_atmos_for_ocean_view<Kokkos::mdspan, Kokkos::dextents> p_as_view;
static struct t_sea_ice_view<Kokkos::mdspan, Kokkos::dextents> p_sea_ice_view;
static struct t_tke_internal_view<Kokkos::mdspan, Kokkos::dextents> p_internal_view;

TKE_cpu::TKE_cpu(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                   int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
                   int l_lc, double clc, double ReferencePressureIndbars, double pi)
    : TKE_backend(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
                  vert_cor_type, dtime, OceanReferenceDensity, grav,
                  l_lc, clc, ReferencePressureIndbars, pi) {
    // Allocate internal arrays memory and create memory views
    std::cout << "Initializing TKE cpu... " << std::endl;

    this->internal_fields_malloc<Kokkos::mdspan, Kokkos::dextents, cpu_mdspan_impl>
                                (&p_internal_view);
    is_view_init = false;
}

TKE_cpu::~TKE_cpu() {
    // Free internal arrays memory
    std::cout << "Finalizing TKE cpu... " << std::endl;

    this->internal_fields_free<cpu_mdspan_impl>();
}

void TKE_cpu::calc_impl(t_patch p_patch, t_cvmix p_cvmix,
                         t_ocean_state ocean_state, t_atmo_fluxes atmos_fluxes,
                         t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
                         int edges_block_size, int edges_start_block, int edges_end_block,
                         int edges_start_index, int edges_end_index, int cells_block_size,
                         int cells_start_block, int cells_end_block, int cells_start_index,
                         int cells_end_index) {
}
