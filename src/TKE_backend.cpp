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

#include "src/TKE_backend.hpp"
#include <iostream>

TKE_backend::TKE_backend(int nproma, int nlevs, int nblocks)
    : m_nproma(nproma), m_nlevs(nlevs), m_nblocks(nblocks) {
}

void TKE_backend::calc(struct t_patch p_patch, struct t_cvmix p_cvmix,
                       struct t_ocean_state ocean_state, struct t_atmo_fluxes atmos_fluxes,
                       struct t_atmos_for_ocean p_as, struct t_sea_ice p_sea_ice,
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
