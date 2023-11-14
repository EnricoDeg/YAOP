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


void fill_struct_view(struct t_cvmix_view *p_cvmix_view, struct t_cvmix *p_cvmix,
                             int nblocks, int nlevs, int nproma) {
    p_cvmix_view->tke = mdspan_3d_double{ p_cvmix->tke, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_plc = mdspan_3d_double{ p_cvmix->tke_plc, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->hlc = mdspan_2d_double{ p_cvmix->hlc, ext2d_t{nblocks, nproma} };
    p_cvmix_view->wlc = mdspan_3d_double{ p_cvmix->wlc, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->u_stokes = mdspan_2d_double{ p_cvmix->u_stokes, ext2d_t{nblocks, nproma} };
    p_cvmix_view->a_veloc_v = mdspan_3d_double{ p_cvmix->a_veloc_v, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->a_temp_v = mdspan_3d_double{ p_cvmix->a_temp_v, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->a_salt_v = mdspan_3d_double{ p_cvmix->a_salt_v, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->iwe_Tdis = mdspan_3d_double{ p_cvmix->iwe_Tdis, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->cvmix_dummy_1 = mdspan_3d_double{ p_cvmix->cvmix_dummy_1, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->cvmix_dummy_2 = mdspan_3d_double{ p_cvmix->cvmix_dummy_2, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->cvmix_dummy_3 = mdspan_3d_double{ p_cvmix->cvmix_dummy_3, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Tbpr = mdspan_3d_double{ p_cvmix->tke_Tbpr, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Tspr = mdspan_3d_double{ p_cvmix->tke_Tspr, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Tdif = mdspan_3d_double{ p_cvmix->tke_Tdif, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Tdis = mdspan_3d_double{ p_cvmix->tke_Tdis, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Twin = mdspan_3d_double{ p_cvmix->tke_Twin, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Tiwf = mdspan_3d_double{ p_cvmix->tke_Tiwf, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Tbck = mdspan_3d_double{ p_cvmix->tke_Tbck, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Ttot = mdspan_3d_double{ p_cvmix->tke_Ttot, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Lmix = mdspan_3d_double{ p_cvmix->tke_Lmix, ext3d_t{nblocks, nlevs+1, nproma} };
    p_cvmix_view->tke_Pr = mdspan_3d_double{ p_cvmix->tke_Pr, ext3d_t{nblocks, nlevs+1, nproma} };
}

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

void fill_struct_view(struct t_atmo_fluxes_view *atmos_fluxes_view, struct t_atmo_fluxes *atmos_fluxes,
                             int nblocks, int nlevs, int nproma) {
    atmos_fluxes_view->stress_xw = mdspan_2d_double{ atmos_fluxes->stress_xw, ext2d_t{nblocks, nproma} };
    atmos_fluxes_view->stress_yw = mdspan_2d_double{ atmos_fluxes->stress_yw, ext2d_t{nblocks, nproma} };
}

void fill_struct_view(struct t_atmos_for_ocean_view *p_as_view, struct t_atmos_for_ocean *p_as,
                             int nblocks, int nlevs, int nproma) {
    p_as_view->fu10 = mdspan_2d_double{ p_as->fu10, ext2d_t{nblocks, nproma} };
}

void fill_struct_view(struct t_sea_ice_view *p_sea_ice_view, struct t_sea_ice *p_sea_ice,
                             int nblocks, int nlevs, int nproma) {
    p_sea_ice_view->concsum = mdspan_2d_double{ p_sea_ice->concsum, ext2d_t{nblocks, nproma} };
}

mdspan_1d_double view_cuda_malloc(double *field, size_t dim1) {
    check( cudaMalloc(&field, dim1*sizeof(double)) );
    mdspan_1d_double memview{ field, ext1d_t{dim1} };
    return memview;
}

mdspan_2d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2) {
    check( cudaMalloc(&field, dim1*dim2*sizeof(double)) );
    mdspan_2d_double memview{ field, ext2d_t{dim1, dim2} };
    return memview;
}

mdspan_3d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2, size_t dim3) {
    check( cudaMalloc(&field, dim1*dim2*dim3*sizeof(double)) );
    mdspan_3d_double memview{ field, ext3d_t{dim1, dim2, dim3} };
    return memview;
}
