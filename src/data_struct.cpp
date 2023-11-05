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

#include "src/data_struct.hpp"

void fill_struct(struct t_patch *p_patch, double *depth_CellInterface, double *prism_center_dist_c,
                 double *inv_prism_center_dist_c, double *prism_thick_c, int *dolic_c, int *dolic_e,
                 double *zlev_i, double *wet_c, int *edges_cell_idx, int *edges_cell_blk) {
    p_patch->depth_CellInterface = depth_CellInterface;
    p_patch->prism_center_dist_c = prism_center_dist_c;
    p_patch->inv_prism_center_dist_c = inv_prism_center_dist_c;
    p_patch->prism_thick_c = prism_thick_c;
    p_patch->dolic_c = dolic_c;
    p_patch->dolic_e = dolic_e;
    p_patch->zlev_i = zlev_i;
    p_patch->wet_c = wet_c;
    p_patch->edges_cell_idx = edges_cell_idx;
    p_patch->edges_cell_blk = edges_cell_blk;
}

void fill_struct(struct t_cvmix * p_cvmix, double *tke, double *tke_plc, double *hlc, double *wlc,
                 double *u_stokes, double *a_veloc_v, double *a_temp_v, double *a_salt_v, double *iwe_Tdis,
                 double *cvmix_dummy_1, double *cvmix_dummy_2, double *cvmix_dummy_3, double *tke_Tbpr,
                 double *tke_Tspr, double *tke_Tdif, double *tke_Tdis, double *tke_Twin, double *tke_Tiwf,
                 double *tke_Tbck, double *tke_Ttot, double *tke_Lmix, double *tke_Pr) {
    p_cvmix->tke = tke;
    p_cvmix->tke_plc = tke_plc;
    p_cvmix->hlc = hlc;
    p_cvmix->wlc = wlc;
    p_cvmix->u_stokes = u_stokes;
    p_cvmix->a_veloc_v = a_veloc_v;
    p_cvmix->a_temp_v = a_temp_v;
    p_cvmix->a_salt_v = a_salt_v;
    p_cvmix->iwe_Tdis = iwe_Tdis;
    p_cvmix->cvmix_dummy_1 = cvmix_dummy_1;
    p_cvmix->cvmix_dummy_2 = cvmix_dummy_2;
    p_cvmix->cvmix_dummy_3 = cvmix_dummy_3;
    p_cvmix->tke_Tbpr = tke_Tbpr;
    p_cvmix->tke_Tspr = tke_Tspr;
    p_cvmix->tke_Tdif = tke_Tdif;
    p_cvmix->tke_Tdis = tke_Tdis;
    p_cvmix->tke_Twin = tke_Twin;
    p_cvmix->tke_Tiwf = tke_Tiwf;
    p_cvmix->tke_Tbck = tke_Tbck;
    p_cvmix->tke_Ttot = tke_Ttot;
    p_cvmix->tke_Lmix = tke_Lmix;
    p_cvmix->tke_Pr = tke_Pr;
}

void fill_struct(struct t_ocean_state * ocean_state, double *temp, double *salt, double *stretch_c,
                 double *eta_c, double *p_vn_x1, double *p_vn_x2, double *p_vn_x3) {
    ocean_state->temp = temp;
    ocean_state->salt = salt;
    ocean_state->stretch_c = stretch_c;
    ocean_state->eta_c = eta_c;
    ocean_state->p_vn_x1 = p_vn_x1;
    ocean_state->p_vn_x2 = p_vn_x2;
    ocean_state->p_vn_x3 = p_vn_x3;
}

void fill_struct(struct t_atmo_fluxes *atmo_fluxes, double *stress_xw, double *stress_yw) {
    atmo_fluxes->stress_xw = stress_xw;
    atmo_fluxes->stress_yw = stress_yw;
}

void fill_struct(struct t_atmos_for_ocean *p_as, double *fu10) {
    p_as->fu10 = fu10;
}

void fill_struct(struct t_sea_ice *p_sea_ice, double *concsum) {
    p_sea_ice->concsum = concsum;
}
