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

class t_patch_base {
 public:
    t_patch_base() {}
    ~t_patch_base() {}
};

template <class T>
class t_patch : public t_patch_base {
 public:
    void fill(T *depth_CellInterface_i, T *prism_center_dist_c_i, T *inv_prism_center_dist_c_i,
              T *prism_thick_c_i, int *dolic_c_i, int *dolic_e_i, T *zlev_i_i, T *wet_c_i,
              int *edges_cell_idx_i, int *edges_cell_blk_i) {
        this->depth_CellInterface = depth_CellInterface_i;
        this->prism_center_dist_c = prism_center_dist_c_i;
        this->inv_prism_center_dist_c = inv_prism_center_dist_c_i;
        this->prism_thick_c = prism_thick_c_i;
        this->dolic_c = dolic_c_i;
        this->dolic_e = dolic_e_i;
        this->zlev_i = zlev_i_i;
        this->wet_c = wet_c_i;
        this->edges_cell_idx = edges_cell_idx_i;
        this->edges_cell_blk = edges_cell_blk_i;
    }
    T * get_depth_CellInterface() { return this->depth_CellInterface; }
    T * get_prism_center_dist_c() { return this->prism_center_dist_c; }
    T * get_inv_prism_center_dist_c() { return this->inv_prism_center_dist_c; }
    T * get_prism_thick_c() { return this->prism_thick_c; }
    int * get_dolic_c() { return this->dolic_c; }
    int * get_dolic_e() { return this->dolic_e; }
    T * get_zlev_i() { return this->zlev_i; }
    T * get_wet_c() { return this->wet_c; }
    int * get_edges_cell_idx() { return this->edges_cell_idx; }
    int * get_edges_cell_blk() { return this->edges_cell_blk; }

 protected:
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

class t_cvmix_base {
 public:
    t_cvmix_base() {}
    ~t_cvmix_base() {}
};

template <class T>
class t_cvmix : public t_cvmix_base {
 public:
    void fill(T *tke_i, T *tke_plc_i, T *hlc_i, T *wlc_i, T *u_stokes_i, T *a_veloc_v_i,
              T *a_temp_v_i, T *a_salt_v_i, T *iwe_Tdis_i, T *cvmix_dummy_1_i,
              T *cvmix_dummy_2_i, T *cvmix_dummy_3_i, T *tke_Tbpr_i, T *tke_Tspr_i,
              T *tke_Tdif_i, T *tke_Tdis_i, T *tke_Twin_i, T *tke_Tiwf_i,
              T *tke_Tbck_i, T *tke_Ttot_i, T *tke_Lmix_i, T *tke_Pr_i) {
        this->tke = tke_i;
        this->tke_plc = tke_plc_i;
        this->hlc = hlc_i;
        this->wlc = wlc_i;
        this->u_stokes = u_stokes_i;
        this->a_veloc_v = a_veloc_v_i;
        this->a_temp_v = a_temp_v_i;
        this->a_salt_v = a_salt_v_i;
        this->iwe_Tdis = iwe_Tdis_i;
        this->cvmix_dummy_1 = cvmix_dummy_1_i;
        this->cvmix_dummy_2 = cvmix_dummy_2_i;
        this->cvmix_dummy_3 = cvmix_dummy_3_i;
        this->tke_Tbpr = tke_Tbpr_i;
        this->tke_Tspr = tke_Tspr_i;
        this->tke_Tdif = tke_Tdif_i;
        this->tke_Tdis = tke_Tdis_i;
        this->tke_Twin = tke_Twin_i;
        this->tke_Tiwf = tke_Tiwf_i;
        this->tke_Tbck = tke_Tbck_i;
        this->tke_Ttot = tke_Ttot_i;
        this->tke_Lmix = tke_Lmix_i;
        this->tke_Pr = tke_Pr_i;
    }
    T * get_tke() { return this->tke; }
    T * get_tke_plc() { return this->tke_plc; }
    T * get_hlc() { return this->hlc; }
    T * get_wlc() { return this->wlc; }
    T * get_u_stokes() { return this->u_stokes; }
    T * get_a_veloc_v() { return this->a_veloc_v; }
    T * get_a_temp_v() { return this->a_temp_v; }
    T * get_a_salt_v() { return this->a_salt_v; }
    T * get_iwe_Tdis() { return this->iwe_Tdis; }
    T * get_cvmix_dummy_1() { return this->cvmix_dummy_1; }
    T * get_cvmix_dummy_2() { return this->cvmix_dummy_2; }
    T * get_cvmix_dummy_3() { return this->cvmix_dummy_3; }
    T * get_tke_Tbpr() { return this->tke_Tbpr; }
    T * get_tke_Tspr() { return this->tke_Tspr; }
    T * get_tke_Tdif() { return this->tke_Tdif; }
    T * get_tke_Tdis() { return this->tke_Tdis; }
    T * get_tke_Twin() { return this->tke_Twin; }
    T * get_tke_Tiwf() { return this->tke_Tiwf; }
    T * get_tke_Tbck() { return this->tke_Tbck; }
    T * get_tke_Ttot() { return this->tke_Ttot; }
    T * get_tke_Lmix() { return this->tke_Lmix; }
    T * get_tke_Pr() { return this->tke_Pr; }

 protected:
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

class t_ocean_state_base {
 public:
    t_ocean_state_base() {}
    ~t_ocean_state_base() {}
};

template <class T>
class t_ocean_state : public t_ocean_state_base {
 public:
    void fill(T *temp_i, T *salt_i, T *stretch_c_i, T *eta_c_i,
              T *p_vn_x1_i, T *p_vn_x2_i, T *p_vn_x3_i) {
        this->temp = temp_i;
        this->salt = salt_i;
        this->stretch_c = stretch_c_i;
        this->eta_c = eta_c_i;
        this->p_vn_x1 = p_vn_x1_i;
        this->p_vn_x2 = p_vn_x2_i;
        this->p_vn_x3 = p_vn_x3_i;
    }
    T * get_temp() { return this->temp; }
    T * get_salt() { return this->salt; }
    T * get_stretch_c() { return this->stretch_c; }
    T * get_eta_c() { return this->eta_c; }
    T * get_p_vn_x1() { return this->p_vn_x1; }
    T * get_p_vn_x2() { return this->p_vn_x2; }
    T * get_p_vn_x3() { return this->p_vn_x3; }

 protected:
    T *temp;
    T *salt;
    T *stretch_c;
    T *eta_c;
    T *p_vn_x1;
    T *p_vn_x2;
    T *p_vn_x3;
};

class t_atmo_fluxes_base {
 public:
    t_atmo_fluxes_base() {}
    ~t_atmo_fluxes_base() {}
};

template <class T>
class t_atmo_fluxes : public t_atmo_fluxes_base {
 public:
    void fill(T *stress_xw_i, T* stress_yw_i) {
        this->stress_xw = stress_xw_i;
        this->stress_yw = stress_yw_i;
    }
    T * get_stress_xw() { return this->stress_xw; }
    T * get_stress_yw() { return this->stress_yw; }

 protected:
    T *stress_xw;
    T *stress_yw;
};

class t_atmos_for_ocean_base {
 public:
    t_atmos_for_ocean_base() {}
    ~t_atmos_for_ocean_base() {}
};

template <class T>
class t_atmos_for_ocean : public t_atmos_for_ocean_base {
 public:
    void fill(T *fu10_i) { this->fu10 = fu10_i; }
    T * get_fu10() { return this->fu10; }

 protected:
    T *fu10;
};

class t_sea_ice_base {
 public:
    t_sea_ice_base() {}
    ~t_sea_ice_base() {}
};

template <class T>
class t_sea_ice : public t_sea_ice_base {
 public:
    void fill(T *concsum_i) { this->concsum = concsum_i; }
    T * get_concsum() { return this->concsum; }

 protected:
    T *concsum;
};

#endif  // SRC_SHARED_INTERFACE_DATA_STRUCT_HPP_
