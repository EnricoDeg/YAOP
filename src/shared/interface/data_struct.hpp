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

#ifndef SRC_SHARED_INTERFACE_DATA_STRUCT_HPP_
#define SRC_SHARED_INTERFACE_DATA_STRUCT_HPP_

// TKE constants
template <class T>
struct t_constant {
    int nproma;
    int nblocks;
    int vert_mix_type;
    int vmix_idemix_tke;
    int vert_cor_type;
    T dtime;
    T OceanReferenceDensity;
    T grav;
    int l_lc;
    T clc;
    T ReferencePressureIndbars;
    T pi;
    int nlevs;
};

template <class T>
struct t_constant_tke {
    T c_k;
    T c_eps;
    T cd;
    T alpha_tke;
    T clc;
    T mxl_min;
    T KappaM_min;
    T KappaH_min;
    T KappaM_max;
    T tke_surf_min;
    T tke_min;
    int tke_mxl_choice;
    int handle_old_vals;
    bool only_tke;
    bool use_Kappa_min;
    bool use_ubound_dirichlet;
    bool use_lbound_dirichlet;
};

template <class T>
struct t_patch {
    T *depth_CellInterface;
    T *prism_center_dist_c;
    T *inv_prism_center_dist_c;
    T *prism_thick_c;
    int *dolic_c;
    int *dolic_e;
    T *zlev_i;
    T *wet_c;
    int *edges_cell_idx;
    int *edges_cell_blk;
};

template <class T>
struct t_cvmix {
    T *tke;
    T *tke_plc;
    T *hlc;
    T *wlc;
    T *u_stokes;
    T *a_veloc_v;
    T *a_temp_v;
    T *a_salt_v;
    T *iwe_Tdis;
    T *cvmix_dummy_1;
    T *cvmix_dummy_2;
    T *cvmix_dummy_3;
    T *tke_Tbpr;
    T *tke_Tspr;
    T *tke_Tdif;
    T *tke_Tdis;
    T *tke_Twin;
    T *tke_Tiwf;
    T *tke_Tbck;
    T *tke_Ttot;
    T *tke_Lmix;
    T *tke_Pr;
};

template <class T>
struct t_ocean_state {
    T *temp;
    T *salt;
    T *stretch_c;
    T *eta_c;
    T *p_vn_x1;
    T *p_vn_x2;
    T *p_vn_x3;
};

template <class T>
struct t_atmo_fluxes {
    T *stress_xw;
    T *stress_yw;
};

struct t_atmos_for_ocean {
    double *fu10;
};

struct t_sea_ice {
    double *concsum;
};

/*! \brief Fill grid info data struct from array pointers.
*
*/
template <class T>
void fill_struct(t_patch<T> *p_patch, T *depth_CellInterface, T *prism_center_dist_c,
                 T *inv_prism_center_dist_c, T *prism_thick_c, int *dolic_c, int *dolic_e,
                 T *zlev_i, T *wet_c, int *edges_cell_idx, int *edges_cell_blk) {
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

/*! \brief Fill cvmix data struct from array pointers.
*
*/
template <class T>
void fill_struct(t_cvmix<T> * p_cvmix, T *tke, T *tke_plc, T *hlc, T *wlc,
                 T *u_stokes, T *a_veloc_v, T *a_temp_v, T *a_salt_v, T *iwe_Tdis,
                 T *cvmix_dummy_1, T *cvmix_dummy_2, T *cvmix_dummy_3, T *tke_Tbpr,
                 T *tke_Tspr, T *tke_Tdif, T *tke_Tdis, T *tke_Twin, T *tke_Tiwf,
                 T *tke_Tbck, T *tke_Ttot, T *tke_Lmix, T *tke_Pr) {
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

/*! \brief Fill ocean state data struct from array pointers.
*
*/
template <class T>
void fill_struct(t_ocean_state<T> * ocean_state, T *temp, T *salt, T *stretch_c,
                 T *eta_c, T *p_vn_x1, T *p_vn_x2, T *p_vn_x3) {
    ocean_state->temp = temp;
    ocean_state->salt = salt;
    ocean_state->stretch_c = stretch_c;
    ocean_state->eta_c = eta_c;
    ocean_state->p_vn_x1 = p_vn_x1;
    ocean_state->p_vn_x2 = p_vn_x2;
    ocean_state->p_vn_x3 = p_vn_x3;
}

/*! \brief Fill atmosphere fluxes data struct from array pointers.
*
*/
template <class T>
void fill_struct(t_atmo_fluxes<T> *atmo_fluxes, T *stress_xw, T *stress_yw) {
    atmo_fluxes->stress_xw = stress_xw;
    atmo_fluxes->stress_yw = stress_yw;
}

/*! \brief Fill atmosphere data struct from array pointers.
*
*/
void fill_struct(struct t_atmos_for_ocean *p_as, double *fu10);

/*! \brief Fill sea ice data struct from array pointers.
*
*/
void fill_struct(struct t_sea_ice *p_sea_ice, double *concsum);

#endif  // SRC_SHARED_INTERFACE_DATA_STRUCT_HPP_
