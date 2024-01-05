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

#ifndef SRC_BACKENDS_GPU_TKE_GPU_HPP_
#define SRC_BACKENDS_GPU_TKE_GPU_HPP_

#include <iostream>
#include "src/backends/TKE_backend.hpp"
#include "src/backends/GPU/gpu_kernels.hpp"
#include "src/shared/utils.hpp"

/*! \brief TKE GPU backend class, derived from TKE_backend.
 *
 */
template <class T>
class TKE_gpu : public TKE_backend<T> {
 public:
    /*! \brief TKE_gpu class constructor.
    *
    *   It calls the TKE_backend constructor and the internal_fields_malloc method with
    *   specific memory view and policy.
    */
    TKE_gpu(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
            int vert_cor_type, T dtime, T OceanReferenceDensity, T grav,
            int l_lc, T clc, T ReferencePressureIndbars, T pi)
            : TKE_backend<T>(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
                             vert_cor_type, dtime, OceanReferenceDensity, grav,
                             l_lc, clc, ReferencePressureIndbars, pi) {
        // Allocate internal arrays memory and create memory views
        std::cout << "Initializing TKE (GPU)... " << std::endl;

        TKE_backend<T>::template internal_fields_malloc
                                 <memview_nms::mdspan, memview_nms::dextents, gpu_memview_policy>();
    }

    /*! \brief TKE_gpu class destructor.
    *
    *   It calls internal_fields_free with specific memory view policy.
    */
    ~TKE_gpu() {
        // Free internal arrays memory
        std::cout << "Finalizing TKE (GPU)... " << std::endl;

        TKE_backend<T>::template internal_fields_free<gpu_memview_policy>();
    }

 protected:
    /*! \brief Function to launch the kernel.
    *
    *   It is templated with the launch policy and calls the launch function of the launch_policy class.
    */
    template <class launch_policy>
    void launch_kernel(int threadsPerBlock, int blocksPerGrid, void* func, void **args) {
        launch_policy::launch(threadsPerBlock, blocksPerGrid, func, args);
    }

    /*! \brief GPU implementation of TKE.
    *
    *   It fills the memory view structures during the first call and then compute the
    *   turbulent kinetic energy vertical scheme.
    */
    void calc_impl(t_patch_base *p_patch, t_cvmix_base *p_cvmix,
                   t_ocean_state_base *ocean_state, t_atmo_fluxes_base *atmos_fluxes,
                   t_atmos_for_ocean_base *p_as, t_sea_ice_base *p_sea_ice,
                   int edges_block_size, int edges_start_block, int edges_end_block,
                   int edges_start_index, int edges_end_index, int cells_block_size,
                   int cells_start_block, int cells_end_block, int cells_start_index,
                   int cells_end_index) {
        // The pointer to the data should not change inside the time loop
        // structs view are filled only at the first time step
        if (!this->m_is_view_init) {
            int nblocks = this->p_constant.nblocks;
            int nlevs = this->p_constant.nlevs;
            int nproma = this->p_constant.nproma;
            TKE_backend<T>::template fill_struct_memview<gpu_memview_policy>
                                     (p_cvmix, nblocks, nlevs, nproma);
            TKE_backend<T>::template fill_struct_memview<gpu_memview_policy>
                                     (p_patch, nblocks, nlevs, nproma);
            TKE_backend<T>::template fill_struct_memview<gpu_memview_policy>
                                     (ocean_state, nblocks, nlevs, nproma);
            TKE_backend<T>::template fill_struct_memview<gpu_memview_policy>
                                     (atmos_fluxes, nblocks, nproma);
            TKE_backend<T>::template fill_struct_memview<gpu_memview_policy>
                                     (p_as, nblocks, nproma);
            TKE_backend<T>::template fill_struct_memview<gpu_memview_policy>
                                     (p_sea_ice, nblocks, nproma);
            this->m_is_view_init = true;
        }

        // over cells
        for (int jb = cells_start_block; jb <= cells_end_block; jb++) {
            int start_index, end_index;
            get_index_range(cells_block_size, cells_start_block, cells_end_block,
                            cells_start_index, cells_end_index, jb, &start_index, &end_index);
            int threadsPerBlockI = 512;  // too many registers used for 1024
            int blocksPerGridI = (end_index - start_index) / threadsPerBlockI + 1;
            void *args[] = {&jb, &start_index, &end_index,
                            &this->p_patch_view, &this->p_cvmix_view,
                            &this->ocean_state_view, &this->atmos_fluxes_view,
                            &this->p_as_view, &this->p_sea_ice_view,
                            &this->p_internal_view, &this->p_constant,
                            &this->p_constant_tke};
            this->launch_kernel<gpu_launch_policy>(threadsPerBlockI, blocksPerGridI,
                                                  reinterpret_cast<void *>(calc_impl_cells<T>), args);
        }

        // over edges
        for (int jb = edges_start_block; jb <= edges_end_block; jb++) {
            int start_index, end_index;
            get_index_range(edges_block_size, edges_start_block, edges_end_block,
                            edges_start_index, edges_end_index, jb, &start_index, &end_index);
            int threadsPerBlockI = 1024;
            int blocksPerGridI = (end_index - start_index) / threadsPerBlockI + 1;
            void *args[] = {&jb, &start_index, &end_index,
                            &this->p_patch_view, &this->p_cvmix_view,
                            &this->p_internal_view, &this->p_constant};
            this->launch_kernel<gpu_launch_policy>(threadsPerBlockI, blocksPerGridI,
                                                  reinterpret_cast<void *>(calc_impl_edges<T>), args);
        }
    }
};

#endif  // SRC_BACKENDS_GPU_TKE_GPU_HPP_
