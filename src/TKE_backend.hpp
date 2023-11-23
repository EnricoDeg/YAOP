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

#ifndef SRC_TKE_BACKEND_HPP_
#define SRC_TKE_BACKEND_HPP_

#include <memory>
#include "src/data_struct.hpp"
#include "src/memview_struct.hpp"

class TKE_backend {
 public:
    typedef std::shared_ptr<TKE_backend> Ptr;
    TKE_backend(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
                int l_lc, double clc, double ReferencePressureIndbars, double pi);
    virtual ~TKE_backend() = default;
    void calc(t_patch p_patch, t_cvmix p_cvmix,
              t_ocean_state ocean_state, t_atmo_fluxes atmos_fluxes,
              t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
              int edges_block_size, int edges_start_block, int edges_end_block,
              int edges_start_index, int edges_end_index, int cells_block_size,
              int cells_start_block, int cells_end_block, int cells_start_index,
              int cells_end_index);

 protected:
    virtual void calc_impl(t_patch p_patch, t_cvmix p_cvmix,
                           t_ocean_state ocean_state, t_atmo_fluxes atmos_fluxes,
                           t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
                           int edges_block_size, int edges_start_block, int edges_end_block,
                           int edges_start_index, int edges_end_index, int cells_block_size,
                           int cells_start_block, int cells_end_block, int cells_start_index,
                           int cells_end_index) = 0;

    template <typename memview, typename memview_policy>
    void fill_struct_memview(t_sea_ice_view<memview> *p_sea_ice_view, struct t_sea_ice *p_sea_ice,
                             int nblocks, int nproma) {
        p_sea_ice_view->concsum = memview_policy::memview_2d_impl(p_sea_ice->concsum, nblocks, nproma);
    }

 protected:
    // Structures with parameters
    struct t_constant p_constant;
    struct t_constant_tke p_constant_tke;

    double *m_tke_old;
    double *m_tke_Av;
    double *m_tke_kv;
    double *m_forc_tke_surf_2D;
    double *m_dzw_stretched;
    double *m_dzt_stretched;
    double *m_Nsqr;
    double *m_Ssqr;
    double *m_a_dif;
    double *m_b_dif;
    double *m_c_dif;
    double *m_a_tri;
    double *m_b_tri;
    double *m_c_tri;
    double *m_d_tri;
    double *m_sqrttke;
    double *m_forc;
    double *m_ke;
    double *m_cp;
    double *m_dp;
    double *m_tke_upd;
    double *m_tke_unrest;
};

#endif  // SRC_TKE_BACKEND_HPP_
