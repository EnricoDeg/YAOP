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

struct t_patch {
    double *depth_CellInterface;
    double *prism_center_dist_c;
    double *inv_prism_center_dist_c;
    double *prism_thick_c;
    int *dolic_c;
    int *dolic_e;
    double *zlev_i;
    double *wet_c;
    int *edges_cell_idx;
    int *edges_cell_blk;
};

struct t_cvmix {
    double *tke;
    double *tke_plc;
    double *hlc;
    double *wlc;
    double *u_stokes;
    double *a_veloc_v;
    double *a_temp_v;
    double *a_salt_v;
    double *iwe_Tdis;
    double *cvmix_dummy_1;
    double *cvmix_dummy_2;
    double *cvmix_dummy_3;
    double *tke_Tbpr;
    double *tke_Tspr;
    double *tke_Tdif;
    double *tke_Tdis;
    double *tke_Twin;
    double *tke_Tiwf;
    double *tke_Tbck;
    double *tke_Ttot;
    double *tke_Lmix;
    double *tke_Pr;
};

struct t_ocean_state {
    double *temp;
    double *salt;
    double *stretch_c;
    double *eta_c;
    double *p_vn_x1;
    double *p_vn_x2;
    double *p_vn_x3;
};

struct t_atmo_fluxes {
    double *stress_xw;
    double *stress_yw;
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
void fill_struct(struct t_patch *p_patch, double *depth_CellInterface, double *prism_center_dist_c,
                 double *inv_prism_center_dist_c, double *prism_thick_c, int *dolic_c, int *dolic_e,
                 double *zlev_i, double *wet_c, int *edges_cell_idx, int *edges_cell_blk);

/*! \brief Fill cvmix data struct from array pointers.
*
*/
void fill_struct(struct t_cvmix * p_cvmix, double *tke, double *tke_plc, double *hlc, double *wlc,
                 double *u_stokes, double *a_veloc_v, double *a_temp_v, double *a_salt_v, double *iwe_Tdis,
                 double *cvmix_dummy_1, double *cvmix_dummy_2, double *cvmix_dummy_3, double *tke_Tbpr,
                 double *tke_Tspr, double *tke_Tdif, double *tke_Tdis, double *tke_Twin, double *tke_Tiwf,
                 double *tke_Tbck, double *tke_Ttot, double *tke_Lmix, double *tke_Pr);

/*! \brief Fill ocean state data struct from array pointers.
*
*/
void fill_struct(struct t_ocean_state * ocean_state, double *temp, double *salt, double *stretch_c,
                 double *eta_c, double *p_vn_x1, double *p_vn_x2, double *p_vn_x3);

/*! \brief Fill atmosphere fluxes data struct from array pointers.
*
*/
void fill_struct(struct t_atmo_fluxes *atmo_fluxes, double *stress_xw, double *stress_yw);

/*! \brief Fill atmosphere data struct from array pointers.
*
*/
void fill_struct(struct t_atmos_for_ocean *p_as, double *fu10);

/*! \brief Fill sea ice data struct from array pointers.
*
*/
void fill_struct(struct t_sea_ice *p_sea_ice, double *concsum);

#endif  // SRC_SHARED_INTERFACE_DATA_STRUCT_HPP_
