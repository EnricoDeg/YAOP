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

#ifndef SRC_SHARED_INTERFACE_MEMVIEW_STRUCT_HPP_
#define SRC_SHARED_INTERFACE_MEMVIEW_STRUCT_HPP_

template <class T,
          template <class ...> class memview,
          template <class, size_t> class dext>
struct t_patch_view {
    memview<T, dext<int, 3>> depth_CellInterface;
    memview<T, dext<int, 3>> prism_center_dist_c;
    memview<T, dext<int, 3>> inv_prism_center_dist_c;
    memview<T, dext<int, 3>> prism_thick_c;
    memview<int, dext<int, 2>> dolic_c;
    memview<int, dext<int, 2>> dolic_e;
    memview<T, dext<int, 1>> zlev_i;
    memview<T, dext<int, 3>> wet_c;
    memview<int, dext<int, 3>> edges_cell_idx;
    memview<int, dext<int, 3>> edges_cell_blk;
};

template <class T,
          template <class ...> class memview,
          template <class, size_t> class dext>
struct t_sea_ice_view {
    memview<T, dext<int, 2>> concsum;
};

template <class T,
          template <class ...> class memview,
          template <class, size_t> class dext>
struct t_atmos_for_ocean_view {
    memview<T, dext<int, 2>> fu10;
};

template <class T,
          template <class ...> class memview,
          template <class, size_t> class dext>
struct t_cvmix_view {
    memview<T, dext<int, 3>> tke;
    memview<T, dext<int, 3>> tke_plc;
    memview<T, dext<int, 2>> hlc;
    memview<T, dext<int, 3>> wlc;
    memview<T, dext<int, 2>> u_stokes;
    memview<T, dext<int, 3>> a_veloc_v;
    memview<T, dext<int, 3>> a_temp_v;
    memview<T, dext<int, 3>> a_salt_v;
    memview<T, dext<int, 3>> iwe_Tdis;
    memview<T, dext<int, 3>> cvmix_dummy_1;
    memview<T, dext<int, 3>> cvmix_dummy_2;
    memview<T, dext<int, 3>> cvmix_dummy_3;
    memview<T, dext<int, 3>> tke_Tbpr;
    memview<T, dext<int, 3>> tke_Tspr;
    memview<T, dext<int, 3>> tke_Tdif;
    memview<T, dext<int, 3>> tke_Tdis;
    memview<T, dext<int, 3>> tke_Twin;
    memview<T, dext<int, 3>> tke_Tiwf;
    memview<T, dext<int, 3>> tke_Tbck;
    memview<T, dext<int, 3>> tke_Ttot;
    memview<T, dext<int, 3>> tke_Lmix;
    memview<T, dext<int, 3>> tke_Pr;
};

template <class T,
          template <class ...> class memview,
          template <class, size_t> class dext>
struct t_atmo_fluxes_view {
    memview<T, dext<int, 2>> stress_xw;
    memview<T, dext<int, 2>> stress_yw;
};

template <class T,
          template <class ...> class memview,
          template <class, size_t> class dext>
struct t_ocean_state_view {
    memview<T, dext<int, 3>> temp;
    memview<T, dext<int, 3>> salt;
    memview<T, dext<int, 2>> stretch_c;
    memview<T, dext<int, 2>> eta_c;
    memview<T, dext<int, 3>> p_vn_x1;
    memview<T, dext<int, 3>> p_vn_x2;
    memview<T, dext<int, 3>> p_vn_x3;
};

template <class T,
          template <class ...> class memview,
          template <class, size_t> class dext>
struct t_tke_internal_view {
    memview<T, dext<int, 1>> forc_tke_surf_2D;
    memview<T, dext<int, 2>> dzw_stretched;
    memview<T, dext<int, 2>> dzt_stretched;
    memview<T, dext<int, 2>> tke_old;
    memview<T, dext<int, 3>> tke_Av;
    memview<T, dext<int, 2>> tke_kv;
    memview<T, dext<int, 2>> Nsqr;
    memview<T, dext<int, 2>> Ssqr;
    memview<T, dext<int, 2>> a_dif;
    memview<T, dext<int, 2>> b_dif;
    memview<T, dext<int, 2>> c_dif;
    memview<T, dext<int, 2>> a_tri;
    memview<T, dext<int, 2>> b_tri;
    memview<T, dext<int, 2>> c_tri;
    memview<T, dext<int, 2>> d_tri;
    memview<T, dext<int, 2>> sqrttke;
    memview<T, dext<int, 2>> forc;
    memview<T, dext<int, 2>> ke;
    memview<T, dext<int, 2>> cp;
    memview<T, dext<int, 2>> dp;
    memview<T, dext<int, 2>> tke_upd;
    memview<T, dext<int, 2>> tke_unrest;
};

#endif  // SRC_SHARED_INTERFACE_MEMVIEW_STRUCT_HPP_
