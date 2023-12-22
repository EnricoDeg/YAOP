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

#ifndef SRC_YAOP_HPP_
#define SRC_YAOP_HPP_

#include <iostream>
#include "src/shared/interface/data_struct.hpp"

/*! \brief YAOP main class, part of the library interface.
 *
 */
class YAOP {
 public:
    /*! \brief YAOP main class constructor called in the model initialization.
     *
     *  Internally the backend is selected based on the library configuration,
     *  some constant parameters needed by the selected scheme are set and all the
     *  internal memory is allocated.
     */
    YAOP(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
        int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
        int l_lc, double clc, double ReferencePressureIndbars, double pi);

    YAOP(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
        int vert_cor_type, float dtime, float OceanReferenceDensity, float grav,
        int l_lc, float clc, float ReferencePressureIndbars, float pi);

    /*! \brief YAOP main class destructor called in the model finalization.
     *
     *  Internally the memory is deallocated.
     */
    ~YAOP();

    /*! \brief YAOP main class time loop calculation of tke scheme.
     *
     *  Internally it is calling the tke scheme implementation of the backend selected
     *  during the configuration.
     *  During the first time step the internal data structures are filled and then the
     *  passed arrays are not used, so the model has to make sure that the memory address
     *  does not change during the time loop.
     */
    void calc_tke(double *depth_CellInterface, double *prism_center_dist_c,
              double *inv_prism_center_dist_c, double *prism_thick_c,
              int *dolic_c, int *dolic_e, double *zlev_i, double *wet_c,
              int *edges_cell_idx, int *edges_cell_blk,
              double *temp, double *salt, double *stretch_c, double *eta_c,
              double *p_vn_x1, double *p_vn_x2, double *p_vn_x3,
              double *tke, double *tke_plc_in, double *hlc_in, double *wlc_in,
              double *u_stokes_in, double *a_veloc_v, double *a_temp_v, double *a_salt_v,
              double *iwe_Tdis, double *cvmix_dummy_1, double *cvmix_dummy_2,
              double *cvmix_dummy_3, double *tke_Tbpr, double *tke_Tspr,
              double *tke_Tdif, double *tke_Tdis, double *tke_Twin,
              double *tke_Tiwf, double *tke_Tbck, double *tke_Ttot,
              double *tke_Lmix, double *tke_Pr, double *stress_xw,
              double *stress_yw, double *fu10, double *concsum,
              int edges_block_size, int edges_start_block, int edges_end_block,
              int edges_start_index, int edges_end_index, int cells_block_size,
              int cells_start_block, int cells_end_block, int cells_start_index,
              int cells_end_index);

    void calc_tke(float *depth_CellInterface, float *prism_center_dist_c,
              float *inv_prism_center_dist_c, float *prism_thick_c,
              int *dolic_c, int *dolic_e, float *zlev_i, float *wet_c,
              int *edges_cell_idx, int *edges_cell_blk,
              float *temp, float *salt, float *stretch_c, float *eta_c,
              float *p_vn_x1, float *p_vn_x2, float *p_vn_x3,
              float *tke, float *tke_plc_in, float *hlc_in, float *wlc_in,
              float *u_stokes_in, float *a_veloc_v, float *a_temp_v, float *a_salt_v,
              float *iwe_Tdis, float *cvmix_dummy_1, float *cvmix_dummy_2,
              float *cvmix_dummy_3, float *tke_Tbpr, float *tke_Tspr,
              float *tke_Tdif, float *tke_Tdis, float *tke_Twin,
              float *tke_Tiwf, float *tke_Tbck, float *tke_Ttot,
              float *tke_Lmix, float *tke_Pr, float *stress_xw,
              float *stress_yw, float *fu10, float *concsum,
              int edges_block_size, int edges_start_block, int edges_end_block,
              int edges_start_index, int edges_end_index, int cells_block_size,
              int cells_start_block, int cells_end_block, int cells_start_index,
              int cells_end_index);

    void calc_vertical_stability();

    void calc_pp();

    void calc_idemix();

 private:
    struct Impl;
    Impl *m_impl;
    bool m_is_struct_init;
};

#endif  // SRC_YAOP_HPP_
