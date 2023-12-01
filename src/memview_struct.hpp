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

#ifndef SRC_MEMVIEW_STRUCT_HPP_
#define SRC_MEMVIEW_STRUCT_HPP_

template <typename T2d>
struct t_sea_ice_view {
    T2d concsum;
};

template <typename T2d>
struct t_atmos_for_ocean_view {
    T2d fu10;
};

template <typename T2d, typename T3d>
struct t_cvmix_view {
    T3d tke;
    T3d tke_plc;
    T2d hlc;
    T3d wlc;
    T2d u_stokes;
    T3d a_veloc_v;
    T3d a_temp_v;
    T3d a_salt_v;
    T3d iwe_Tdis;
    T3d cvmix_dummy_1;
    T3d cvmix_dummy_2;
    T3d cvmix_dummy_3;
    T3d tke_Tbpr;
    T3d tke_Tspr;
    T3d tke_Tdif;
    T3d tke_Tdis;
    T3d tke_Twin;
    T3d tke_Tiwf;
    T3d tke_Tbck;
    T3d tke_Ttot;
    T3d tke_Lmix;
    T3d tke_Pr;
};

template <typename T2d>
struct t_atmo_fluxes_view {
    T2d stress_xw;
    T2d stress_yw;
};

template <typename T2d, typename T3d>
struct t_ocean_state_view {
    T3d temp;
    T3d salt;
    T2d stretch_c;
    T2d eta_c;
    T3d p_vn_x1;
    T3d p_vn_x2;
    T3d p_vn_x3;
};

template <typename T1d, typename T2d, typename T3d>
struct t_tke_internal_view {
    T1d forc_tke_surf_2D;
    T2d dzw_stretched;
    T2d dzt_stretched;
    T2d tke_old;
    T3d tke_Av;
    T2d tke_kv;
    T2d Nsqr;
    T2d Ssqr;
    T2d a_dif;
    T2d b_dif;
    T2d c_dif;
    T2d a_tri;
    T2d b_tri;
    T2d c_tri;
    T2d d_tri;
    T2d sqrttke;
    T2d forc;
    T2d ke;
    T2d cp;
    T2d dp;
    T2d tke_upd;
    T2d tke_unrest;
};


#endif  // SRC_MEMVIEW_STRUCT_HPP_
