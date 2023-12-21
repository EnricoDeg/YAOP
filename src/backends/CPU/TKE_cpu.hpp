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

#include <iostream>
#include "src/backends/TKE_backend.hpp"
#include "src/backends/CPU/cpu_kernels.hpp"
#include "src/shared/utils.hpp"

/*! \brief TKE CPU backend class, derived from TKE_backend.
 *
 */
template <class T>
class TKE_cpu : public TKE_backend {
 public:
    /*! \brief TKE_cpu class constructor.
    *
    *   It calls the TKE_backend constructor and the internal_fields_malloc method with
    *   specific memory view and policy.
    */
    TKE_cpu(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
            int vert_cor_type, T dtime, T OceanReferenceDensity, T grav,
            int l_lc, T clc, T ReferencePressureIndbars, T pi)
            : TKE_backend(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
                          vert_cor_type, dtime, OceanReferenceDensity, grav,
                          l_lc, clc, ReferencePressureIndbars, pi) {
        // Allocate internal arrays memory and create memory views
        std::cout << "Initializing TKE cpu... " << std::endl;

        this->internal_fields_malloc<T, memview_nms::mdspan, memview_nms::dextents, cpu_mdspan_impl>
                                    (&p_internal_view);
    }

    /*! \brief TKE_cpu class destructor.
    *
    *   It calls internal_fields_free with specific memory view policy.
    */
    ~TKE_cpu() {
        // Free internal arrays memory
        std::cout << "Finalizing TKE cpu... " << std::endl;

        this->internal_fields_free<T, cpu_mdspan_impl>();
    }

 protected:
    /*! \brief CPU implementation of TKE.
    *
    *   It fills the memory view structures during the first call and then compute the
    *   turbulent kinetic energy vertical scheme.
    */
    void calc_impl(t_patch<T> p_patch, t_cvmix<T> p_cvmix,
                   t_ocean_state<T> ocean_state, t_atmo_fluxes atmos_fluxes,
                   t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
                   int edges_block_size, int edges_start_block, int edges_end_block,
                   int edges_start_index, int edges_end_index, int cells_block_size,
                   int cells_start_block, int cells_end_block, int cells_start_index,
                   int cells_end_index) {
        // The pointer to the data should not change inside the time loop
        // structs view are filled only at the first time step
        if (!m_is_view_init) {
            this->fill_struct_memview<T, memview_nms::mdspan, memview_nms::dextents, cpu_memview_policy>
                                     (&p_cvmix_view, &p_cvmix, p_constant.nblocks, p_constant.nlevs, p_constant.nproma);
            this->fill_struct_memview<T, memview_nms::mdspan, memview_nms::dextents, cpu_memview_policy>
                                     (&p_patch_view, &p_patch, p_constant.nblocks, p_constant.nlevs, p_constant.nproma);
            this->fill_struct_memview<T, memview_nms::mdspan, memview_nms::dextents, cpu_memview_policy>
                                     (&ocean_state_view, &ocean_state, p_constant.nblocks, p_constant.nlevs,
                                      p_constant.nproma);
            this->fill_struct_memview<T, memview_nms::mdspan, memview_nms::dextents, cpu_memview_policy>
                                     (&atmos_fluxes_view, &atmos_fluxes, p_constant.nblocks, p_constant.nproma);
            this->fill_struct_memview<T, memview_nms::mdspan, memview_nms::dextents, cpu_memview_policy>
                                     (&p_as_view, &p_as, p_constant.nblocks, p_constant.nproma);
            this->fill_struct_memview<T, memview_nms::mdspan, memview_nms::dextents, cpu_memview_policy>
                                     (&p_sea_ice_view, &p_sea_ice, p_constant.nblocks, p_constant.nproma);
            m_is_view_init = true;
        }

        // over cells
        for (int jb = cells_start_block; jb <= cells_end_block; jb++) {
            int start_index, end_index;
            get_index_range(cells_block_size, cells_start_block, cells_end_block,
                            cells_start_index, cells_end_index, jb, &start_index, &end_index);
            calc_impl_cells<T>(jb, start_index, end_index,
                                    p_patch_view, p_cvmix_view,
                                    ocean_state_view, atmos_fluxes_view,
                                    p_as_view, p_sea_ice_view,
                                    p_internal_view, p_constant,
                                    p_constant_tke);
        }

        // over edges
        for (int jb = edges_start_block; jb <= edges_end_block; jb++) {
            int start_index, end_index;
            get_index_range(edges_block_size, edges_start_block, edges_end_block,
                            edges_start_index, edges_end_index, jb, &start_index, &end_index);
            calc_impl_edges<T>(jb, start_index, end_index,
                            p_patch_view, p_cvmix_view,
                            p_internal_view, p_constant);
        }
    }
};

#endif  // SRC_BACKENDS_CPU_TKE_CPU_HPP_
