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

#include "src/TKE.hpp"

#include <iostream>

#include "src/TKE_backend.hpp"
#include "src/TKE_cuda.hpp"
#include "src/data_struct.hpp"

struct TKE::Impl {
  TKE_backend::Ptr backend;
};

TKE::TKE(int nproma, int nlevs, int nblocks,
         int block_size, int start_index, int end_index)
    : m_impl(new Impl) {
    std::cout << "Initializing TKE... " << std::endl;
    m_impl->backend = TKE_backend::Ptr(new TKE_cuda(nproma, nlevs, nblocks,
                                                    block_size, start_index, end_index));
}

TKE::~TKE() {
    std::cout << "Finalizing TKE... " << std::endl;
    delete m_impl;
}

void TKE::calc(int start_block, int end_block,
               double *depth_CellInterface, double *prism_center_dist_c,
               double *inv_prism_center_dist_c, double *prism_thick_c,
               int *dolic_c, int *dolic_e, double *zlev_i, double *wet_c,
               int *edges_cell_idx, int *edges_cell_blk,
               double *temp, double *salt, double *stretch_c, double *eta_c,
               double *tke, double *tke_plc_in, double *hlc_in, double *wlc_in,
               double *u_stokes_in, double *a_veloc_v, double *a_temp_v, double *a_salt_v,
               double *iwe_Tdis, double *cvmix_dummy_1, double *cvmix_dummy_2,
               double *cvmix_dummy_3, double *tke_Tbpr, double *tke_Tspr,
               double *tke_Tdif, double *tke_Tdis, double *tke_Twin,
               double *tke_Tiwf, double *tke_Tbck, double *tke_Ttot,
               double *tke_Lmix, double *tke_Pr, double *stress_xw,
               double *stress_yw, double *fu10, double *concsum) {
    struct t_patch p_patch;
    p_patch.dolic_c = dolic_c;
    struct t_cvmix p_cvmix;
    p_cvmix.tke = tke;
    struct t_sea_ice p_sea_ice;
    struct t_atmos_for_ocean p_as;
    struct t_atmo_fluxes atmos_fluxes;
    struct t_ocean_state ocean_state;

    m_impl->backend->calc(start_block, end_block, p_patch, p_cvmix, ocean_state,
                          atmos_fluxes, p_as, p_sea_ice);
}
