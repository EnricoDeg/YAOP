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

#include "src/backends/TKE_backend.hpp"
#include <iostream>

TKE_backend::TKE_backend(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                         int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
                         int l_lc, double clc, double ReferencePressureIndbars, double pi) {
    // Fill structures with parameters
    p_constant.nproma = nproma;
    p_constant.nblocks = nblocks;
    p_constant.vert_mix_type = vert_mix_type;
    p_constant.vmix_idemix_tke = vmix_idemix_tke;
    p_constant.vert_cor_type = vert_cor_type;
    p_constant.dtime = dtime;
    p_constant.OceanReferenceDensity = OceanReferenceDensity;
    p_constant.grav = grav;
    p_constant.l_lc = l_lc;
    p_constant.clc = clc;
    p_constant.ReferencePressureIndbars = ReferencePressureIndbars;
    p_constant.pi = pi;
    p_constant.nlevs = nlevs;

    // Internal parameters are set for now to default values
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

    m_is_view_init = false;
}

void TKE_backend::calc(t_patch<double> p_patch, t_cvmix<double> p_cvmix,
                       t_ocean_state<double> ocean_state, t_atmo_fluxes<double> atmos_fluxes,
                       t_atmos_for_ocean<double> p_as, t_sea_ice<double> p_sea_ice,
                       int edges_block_size, int edges_start_block, int edges_end_block,
                       int edges_start_index, int edges_end_index, int cells_block_size,
                       int cells_start_block, int cells_end_block, int cells_start_index,
                       int cells_end_index) {
    this->calc_impl(p_patch, p_cvmix, ocean_state, atmos_fluxes, p_as, p_sea_ice,
                    edges_block_size, edges_start_block, edges_end_block,
                    edges_start_index, edges_end_index, cells_block_size,
                    cells_start_block, cells_end_block, cells_start_index,
                    cells_end_index);
}
