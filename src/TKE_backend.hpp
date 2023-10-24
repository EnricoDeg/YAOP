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

class TKE_backend {
 public:
    typedef std::shared_ptr<TKE_backend> Ptr;
    TKE_backend(int nproma, int nlevs, int nblocks,
                int block_size, int start_index, int end_index);
    virtual ~TKE_backend() = default;
    void calc(int start_block, int end_block, struct t_patch p_patch, struct t_cvmix p_cvmix,
              struct t_ocean_state ocean_state, struct t_atmo_fluxes atmos_fluxes,
              struct t_atmos_for_ocean p_as, struct t_sea_ice p_sea_ice);

 protected:
    virtual void calc_impl(int start_block, int end_block, struct t_patch p_patch, struct t_cvmix p_cvmix,
                           struct t_ocean_state ocean_state, struct t_atmo_fluxes atmos_fluxes,
                           struct t_atmos_for_ocean p_as, struct t_sea_ice p_sea_ice) = 0;

 protected:
    int m_nproma;
    int m_nlevs;
    int m_nblocks;
    int m_block_size;
    int m_start_index;
    int m_end_index;

    double *m_rho_up;
    double *m_rho_down;
    double *m_tke_old;
    double *m_tke_Av;
    double *m_tke_kv;
    double *m_tke_iw_alpha_c;
    double *m_tke_iwe;
    double *m_tke_iwe_forcing;
    double *m_forc_tke_surf_2D;
    double *m_forc_rho_surf_2D;
    double *m_bottom_fric_2D;
    double *m_s_c;
    double *m_dzw_stretched;
    double *m_dzt_stretched;
    double *m_pressure;
    double *m_Nsqr;
    double *m_Ssqr;
};

#endif  // SRC_TKE_BACKEND_HPP_
