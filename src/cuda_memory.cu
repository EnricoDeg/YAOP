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

#include "src/cuda_memory.hpp"

void fill_struct_view(struct t_patch_view *p_patch_view, struct t_patch *p_patch,
                             int nblocks, int nlevs, int nproma) {
    p_patch_view->depth_CellInterface = mdspan_3d_double{ p_patch->depth_CellInterface,
                                                          ext3d_t{nblocks, nlevs+1, nproma} };
    p_patch_view->prism_center_dist_c = mdspan_3d_double{ p_patch->prism_center_dist_c,
                                                          ext3d_t{nblocks, nlevs+1, nproma} };
    p_patch_view->inv_prism_center_dist_c = mdspan_3d_double{ p_patch->inv_prism_center_dist_c,
                                                              ext3d_t{nblocks, nlevs+1, nproma} };
    p_patch_view->prism_thick_c = mdspan_3d_double{ p_patch->prism_thick_c, ext3d_t{nblocks, nlevs, nproma} };
    p_patch_view->dolic_c = mdspan_2d_int{ p_patch->dolic_c, ext2d_t{nblocks, nproma} };
    p_patch_view->dolic_e = mdspan_2d_int{ p_patch->dolic_e, ext2d_t{nblocks, nproma} };
    p_patch_view->zlev_i = mdspan_1d_double{ p_patch->zlev_i, ext1d_t{nlevs} };
    p_patch_view->wet_c = mdspan_3d_double{ p_patch->wet_c, ext3d_t{nblocks, nlevs, nproma} };
    p_patch_view->edges_cell_idx = mdspan_3d_int{ p_patch->edges_cell_idx, ext3d_t{2, nlevs, nproma} };
    p_patch_view->edges_cell_blk = mdspan_3d_int{ p_patch->edges_cell_blk, ext3d_t{2, nlevs, nproma} };
}


void fill_struct_view(struct t_ocean_state_view *ocean_state_view, struct t_ocean_state *ocean_state,
                             int nblocks, int nlevs, int nproma) {
    ocean_state_view->temp = mdspan_3d_double{ ocean_state->temp, ext3d_t{nblocks, nlevs, nproma} };
    ocean_state_view->salt = mdspan_3d_double{ ocean_state->salt, ext3d_t{nblocks, nlevs, nproma} };
    ocean_state_view->stretch_c = mdspan_2d_double{ ocean_state->stretch_c, ext2d_t{nblocks, nproma} };
    ocean_state_view->eta_c = mdspan_2d_double{ ocean_state->eta_c, ext2d_t{nblocks, nproma} };
    ocean_state_view->p_vn_x1 = mdspan_3d_double{ ocean_state->p_vn_x1, ext3d_t{nblocks, nlevs, nproma} };
    ocean_state_view->p_vn_x2 = mdspan_3d_double{ ocean_state->p_vn_x2, ext3d_t{nblocks, nlevs, nproma} };
    ocean_state_view->p_vn_x3 = mdspan_3d_double{ ocean_state->p_vn_x3, ext3d_t{nblocks, nlevs, nproma} };
}
