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

template <template <class ...> class memview,
          template <class, size_t> class dext>
struct t_patch_view {
    memview<double, dext<int, 3>> depth_CellInterface;
    memview<double, dext<int, 3>> prism_center_dist_c;
    memview<double, dext<int, 3>> inv_prism_center_dist_c;
    memview<double, dext<int, 3>> prism_thick_c;
    memview<int, dext<int, 2>> dolic_c;
    memview<int, dext<int, 2>> dolic_e;
    memview<double, dext<int, 1>> zlev_i;
    memview<double, dext<int, 3>> wet_c;
    memview<int, dext<int, 3>> edges_cell_idx;
    memview<int, dext<int, 3>> edges_cell_blk;
};

template <template <class ...> class memview,
          template <class, size_t> class dext>
struct t_sea_ice_view {
    memview<double, dext<int, 2>> concsum;
};

template <template <class ...> class memview,
          template <class, size_t> class dext>
struct t_atmos_for_ocean_view {
    memview<double, dext<int, 2>> fu10;
};

template <template <class ...> class memview,
          template <class, size_t> class dext>
struct t_cvmix_view {
    memview<double, dext<int, 3>> tke;
    memview<double, dext<int, 3>> tke_plc;
    memview<double, dext<int, 2>> hlc;
    memview<double, dext<int, 3>> wlc;
    memview<double, dext<int, 2>> u_stokes;
    memview<double, dext<int, 3>> a_veloc_v;
    memview<double, dext<int, 3>> a_temp_v;
    memview<double, dext<int, 3>> a_salt_v;
    memview<double, dext<int, 3>> iwe_Tdis;
    memview<double, dext<int, 3>> cvmix_dummy_1;
    memview<double, dext<int, 3>> cvmix_dummy_2;
    memview<double, dext<int, 3>> cvmix_dummy_3;
    memview<double, dext<int, 3>> tke_Tbpr;
    memview<double, dext<int, 3>> tke_Tspr;
    memview<double, dext<int, 3>> tke_Tdif;
    memview<double, dext<int, 3>> tke_Tdis;
    memview<double, dext<int, 3>> tke_Twin;
    memview<double, dext<int, 3>> tke_Tiwf;
    memview<double, dext<int, 3>> tke_Tbck;
    memview<double, dext<int, 3>> tke_Ttot;
    memview<double, dext<int, 3>> tke_Lmix;
    memview<double, dext<int, 3>> tke_Pr;
};

template <template <class ...> class memview,
          template <class, size_t> class dext>
struct t_atmo_fluxes_view {
    memview<double, dext<int, 2>> stress_xw;
    memview<double, dext<int, 2>> stress_yw;
};

template <template <class ...> class memview,
          template <class, size_t> class dext>
struct t_ocean_state_view {
    memview<double, dext<int, 3>> temp;
    memview<double, dext<int, 3>> salt;
    memview<double, dext<int, 2>> stretch_c;
    memview<double, dext<int, 2>> eta_c;
    memview<double, dext<int, 3>> p_vn_x1;
    memview<double, dext<int, 3>> p_vn_x2;
    memview<double, dext<int, 3>> p_vn_x3;
};

template <template <class ...> class memview,
          template <class, size_t> class dext>
struct t_tke_internal_view {
    memview<double, dext<int, 1>> forc_tke_surf_2D;
    memview<double, dext<int, 2>> dzw_stretched;
    memview<double, dext<int, 2>> dzt_stretched;
    memview<double, dext<int, 2>> tke_old;
    memview<double, dext<int, 3>> tke_Av;
    memview<double, dext<int, 2>> tke_kv;
    memview<double, dext<int, 2>> Nsqr;
    memview<double, dext<int, 2>> Ssqr;
    memview<double, dext<int, 2>> a_dif;
    memview<double, dext<int, 2>> b_dif;
    memview<double, dext<int, 2>> c_dif;
    memview<double, dext<int, 2>> a_tri;
    memview<double, dext<int, 2>> b_tri;
    memview<double, dext<int, 2>> c_tri;
    memview<double, dext<int, 2>> d_tri;
    memview<double, dext<int, 2>> sqrttke;
    memview<double, dext<int, 2>> forc;
    memview<double, dext<int, 2>> ke;
    memview<double, dext<int, 2>> cp;
    memview<double, dext<int, 2>> dp;
    memview<double, dext<int, 2>> tke_upd;
    memview<double, dext<int, 2>> tke_unrest;
};

#endif  // SRC_MEMVIEW_STRUCT_HPP_
