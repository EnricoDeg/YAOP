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
