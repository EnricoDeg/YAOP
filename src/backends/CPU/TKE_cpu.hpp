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

#ifndef SRC_BACKENDS_CPU_TKE_CPU_HPP_
#define SRC_BACKENDS_CPU_TKE_CPU_HPP_

#include "src/backends/TKE_backend.hpp"

/*! \brief TKE CPU backend class, derived from TKE_backend.
 *
 */
class TKE_cpu : public TKE_backend {
 public:
    /*! \brief TKE_cpu class constructor.
    *
    *   It calls the TKE_backend constructor and the internal_fields_malloc method with
    *   specific memory view and policy.
    */
    TKE_cpu(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
             int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
             int l_lc, double clc, double ReferencePressureIndbars, double pi);

    /*! \brief TKE_cpu class destructor.
    *
    *   It calls internal_fields_free with specific memory view policy.
    */
    ~TKE_cpu();

 protected:
    /*! \brief CPU implementation of TKE.
    *
    *   It fills the memory view structures during the first call and then compute the
    *   turbulent kinetic energy vertical scheme.
    */
    void calc_impl(struct t_patch p_patch, struct t_cvmix p_cvmix,
                   struct t_ocean_state ocean_state, struct t_atmo_fluxes atmos_fluxes,
                   struct t_atmos_for_ocean p_as, struct t_sea_ice p_sea_ice,
                   int edges_block_size, int edges_start_block, int edges_end_block,
                   int edges_start_index, int edges_end_index, int cells_block_size,
                   int cells_start_block, int cells_end_block, int cells_start_index,
                   int cells_end_index);
};

#endif  // SRC_BACKENDS_CPU_TKE_CPU_HPP_
