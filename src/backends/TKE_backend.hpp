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

#ifndef SRC_BACKENDS_TKE_BACKEND_HPP_
#define SRC_BACKENDS_TKE_BACKEND_HPP_

#include <memory>
#include "src/shared/interface/data_struct.hpp"
#include "src/shared/interface/memview_struct.hpp"

/*! \brief TKE backend class.
 *
 *  It implements the part which is in common for each backend:
 *   - library frontend-backend interface
 *   - library initialization and finalization
 */
class TKE_backend {
 public:
    typedef std::shared_ptr<TKE_backend> Ptr;

    /*! \brief TKE backend class constructor.
    */
    TKE_backend(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
                int l_lc, double clc, double ReferencePressureIndbars, double pi);

    /*! \brief TKE backend class destructor.
    */
    virtual ~TKE_backend() = default;

    /*! \brief TKE backend calculation.
    *
    *  It calls calc_impl which contains the actual implementation.
    */
    void calc(t_patch p_patch, t_cvmix p_cvmix,
              t_ocean_state ocean_state, t_atmo_fluxes atmos_fluxes,
              t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
              int edges_block_size, int edges_start_block, int edges_end_block,
              int edges_start_index, int edges_end_index, int cells_block_size,
              int cells_start_block, int cells_end_block, int cells_start_index,
              int cells_end_index);

 protected:
    /*! \brief Polymorphic function for the actual TKE scheme backend implementation.
    *
    */
    virtual void calc_impl(t_patch p_patch, t_cvmix p_cvmix,
                           t_ocean_state ocean_state, t_atmo_fluxes atmos_fluxes,
                           t_atmos_for_ocean p_as, t_sea_ice p_sea_ice,
                           int edges_block_size, int edges_start_block, int edges_end_block,
                           int edges_start_index, int edges_end_index, int cells_block_size,
                           int cells_start_block, int cells_end_block, int cells_start_index,
                           int cells_end_index) = 0;

    /*! \brief fill a structure of memory views given a structure of pointers about the grid info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    void fill_struct_memview(t_patch_view<memview, dext> *p_patch_view,
                             t_patch *p_patch, int nblocks, int nlevs, int nproma) {
        p_patch_view->depth_CellInterface = memview_policy::memview(p_patch->depth_CellInterface,
                                                                    nblocks, nlevs+1, nproma);
        p_patch_view->prism_center_dist_c = memview_policy::memview(p_patch->prism_center_dist_c,
                                                                    nblocks, nlevs+1, nproma);
        p_patch_view->inv_prism_center_dist_c = memview_policy::memview(p_patch->inv_prism_center_dist_c,
                                                                        nblocks, nlevs+1, nproma);
        p_patch_view->prism_thick_c = memview_policy::memview(p_patch->prism_thick_c, nblocks, nlevs, nproma);
        p_patch_view->dolic_c = memview_policy::memview(p_patch->dolic_c, nblocks, nproma);
        p_patch_view->dolic_e = memview_policy::memview(p_patch->dolic_e, nblocks, nproma);
        p_patch_view->zlev_i = memview_policy::memview(p_patch->zlev_i, nlevs);
        p_patch_view->wet_c = memview_policy::memview(p_patch->wet_c, nblocks, nlevs, nproma);
        p_patch_view->edges_cell_idx = memview_policy::memview(p_patch->edges_cell_idx, 2, nlevs, nproma);
        p_patch_view->edges_cell_blk = memview_policy::memview(p_patch->edges_cell_blk, 2, nlevs, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the sea ice info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    void fill_struct_memview(t_sea_ice_view<memview, dext> *p_sea_ice_view, t_sea_ice *p_sea_ice,
                             int nblocks, int nproma) {
        p_sea_ice_view->concsum = memview_policy::memview(p_sea_ice->concsum, nblocks, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the atmosphere info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    void fill_struct_memview(t_atmos_for_ocean_view<memview, dext> *p_as_view, t_atmos_for_ocean *p_as,
                             int nblocks, int nproma) {
        p_as_view->fu10 = memview_policy::memview(p_as->fu10, nblocks, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the atmosphere fluxes info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    void fill_struct_memview(t_atmo_fluxes_view<memview, dext> *atmos_fluxes_view, t_atmo_fluxes *atmos_fluxes,
                             int nblocks, int nproma) {
        atmos_fluxes_view->stress_xw = memview_policy::memview(atmos_fluxes->stress_xw, nblocks, nproma);
        atmos_fluxes_view->stress_yw = memview_policy::memview(atmos_fluxes->stress_yw, nblocks, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the ocean state info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    void fill_struct_memview(t_ocean_state_view<memview, dext> *ocean_state_view,
                             t_ocean_state *ocean_state, int nblocks, int nlevs, int nproma) {
        ocean_state_view->temp = memview_policy::memview(ocean_state->temp, nblocks, nlevs, nproma);
        ocean_state_view->salt = memview_policy::memview(ocean_state->salt, nblocks, nlevs, nproma);
        ocean_state_view->stretch_c = memview_policy::memview(ocean_state->stretch_c, nblocks, nproma);
        ocean_state_view->eta_c = memview_policy::memview(ocean_state->eta_c, nblocks, nproma);
        ocean_state_view->p_vn_x1 = memview_policy::memview(ocean_state->p_vn_x1, nblocks, nlevs, nproma);
        ocean_state_view->p_vn_x2 = memview_policy::memview(ocean_state->p_vn_x2, nblocks, nlevs, nproma);
        ocean_state_view->p_vn_x3 = memview_policy::memview(ocean_state->p_vn_x3, nblocks, nlevs, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the cvmix info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    void fill_struct_memview(t_cvmix_view<memview, dext> *p_cvmix_view, t_cvmix *p_cvmix,
                             int nblocks, int nlevs, int nproma) {
        p_cvmix_view->tke = memview_policy::memview(p_cvmix->tke, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_plc = memview_policy::memview(p_cvmix->tke_plc, nblocks, nlevs+1, nproma);
        p_cvmix_view->hlc = memview_policy::memview(p_cvmix->hlc, nblocks, nproma);
        p_cvmix_view->wlc = memview_policy::memview(p_cvmix->wlc, nblocks, nlevs+1, nproma);
        p_cvmix_view->u_stokes = memview_policy::memview(p_cvmix->u_stokes, nblocks, nproma);
        p_cvmix_view->a_veloc_v = memview_policy::memview(p_cvmix->a_veloc_v, nblocks, nlevs+1, nproma);
        p_cvmix_view->a_temp_v = memview_policy::memview(p_cvmix->a_temp_v, nblocks, nlevs+1, nproma);
        p_cvmix_view->a_salt_v = memview_policy::memview(p_cvmix->a_salt_v, nblocks, nlevs+1, nproma);
        p_cvmix_view->iwe_Tdis = memview_policy::memview(p_cvmix->iwe_Tdis, nblocks, nlevs+1, nproma);
        p_cvmix_view->cvmix_dummy_1 = memview_policy::memview(p_cvmix->cvmix_dummy_1, nblocks, nlevs+1, nproma);
        p_cvmix_view->cvmix_dummy_2 = memview_policy::memview(p_cvmix->cvmix_dummy_2, nblocks, nlevs+1, nproma);
        p_cvmix_view->cvmix_dummy_3 = memview_policy::memview(p_cvmix->cvmix_dummy_3, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Tbpr = memview_policy::memview(p_cvmix->tke_Tbpr, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Tspr = memview_policy::memview(p_cvmix->tke_Tspr, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Tdif = memview_policy::memview(p_cvmix->tke_Tdif, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Tdis = memview_policy::memview(p_cvmix->tke_Tdis, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Twin = memview_policy::memview(p_cvmix->tke_Twin, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Tiwf = memview_policy::memview(p_cvmix->tke_Tiwf, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Tbck = memview_policy::memview(p_cvmix->tke_Tbck, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Ttot = memview_policy::memview(p_cvmix->tke_Ttot, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Lmix = memview_policy::memview(p_cvmix->tke_Lmix, nblocks, nlevs+1, nproma);
        p_cvmix_view->tke_Pr = memview_policy::memview(p_cvmix->tke_Pr, nblocks, nlevs+1, nproma);
    }

    /*! \brief allocate internal memory and return a 1D memory view object of the allocated memory.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    memview<double, dext<int, 1>> memview_malloc(double *field, int dim1) {
        return memview_policy::memview_malloc(field, dim1);
    }

    /*! \brief allocate internal memory and return a 2D memory view object of the allocated memory.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    memview<double, dext<int, 2>> memview_malloc(double *field, int dim1, int dim2) {
        return memview_policy::memview_malloc(field, dim1, dim2);
    }

    /*! \brief allocate internal memory and return a 3D memory view object of the allocated memory.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    memview<double, dext<int, 3>> memview_malloc(double *field, int dim1, int dim2, int dim3) {
        return memview_policy::memview_malloc(field, dim1, dim2, dim3);
    }

    /*! \brief deallocate internal memory.
    *
    *   It is templated with a memview_policy which defines how to deallocate memory in the actual backend.
    */
    template <typename memview_policy>
    void memview_free(double *field) {
        return memview_policy::memview_free(field);
    }

    /*! \brief fill the internal data structure allocating the arrays and creating memory views.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              class memview_policy>
    void internal_fields_malloc(t_tke_internal_view<memview, dext> *p_internal_view) {
        p_internal_view->tke_old = this->memview_malloc<memview, dext, memview_policy>
                                         (m_tke_old, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->forc_tke_surf_2D = this->memview_malloc<memview, dext, memview_policy>
                                                  (m_forc_tke_surf_2D, p_constant.nproma);
        p_internal_view->dzw_stretched = this->memview_malloc<memview, dext, memview_policy>
                                               (m_dzw_stretched, p_constant.nlevs, p_constant.nproma);
        p_internal_view->dzt_stretched = this->memview_malloc<memview, dext, memview_policy>
                                               (m_dzt_stretched, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->tke_Av = this->memview_malloc<memview, dext, memview_policy>
                                        (m_tke_Av, p_constant.nblocks, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->tke_kv = this->memview_malloc<memview, dext, memview_policy>
                                        (m_tke_kv, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->Nsqr = this->memview_malloc<memview, dext, memview_policy>
                                      (m_Nsqr, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->Ssqr = this->memview_malloc<memview, dext, memview_policy>
                                      (m_Ssqr, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->a_dif = this->memview_malloc<memview, dext, memview_policy>
                                       (m_a_dif, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->b_dif = this->memview_malloc<memview, dext, memview_policy>
                                       (m_b_dif, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->c_dif = this->memview_malloc<memview, dext, memview_policy>
                                       (m_c_dif, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->a_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_a_tri, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->b_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_b_tri, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->c_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_c_tri, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->d_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_d_tri, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->sqrttke = this->memview_malloc<memview, dext, memview_policy>
                                         (m_sqrttke, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->forc = this->memview_malloc<memview, dext, memview_policy>
                                      (m_forc, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->ke = this->memview_malloc<memview, dext, memview_policy>
                                    (m_ke, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->cp = this->memview_malloc<memview, dext, memview_policy>
                                    (m_cp, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->dp = this->memview_malloc<memview, dext, memview_policy>
                                    (m_dp, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->tke_upd = this->memview_malloc<memview, dext, memview_policy>
                                         (m_tke_upd, p_constant.nlevs+1, p_constant.nproma);
        p_internal_view->tke_unrest = this->memview_malloc<memview, dext, memview_policy>
                                            (m_tke_unrest, p_constant.nlevs+1, p_constant.nproma);
    }

    /*! \brief free the internal data structure memory deallocating the arrays.
    *
    *   It is templated with a memview_policy which defines how to deallocate memory in the actual backend.
    */
    template <typename memview_policy>
    void internal_fields_free() {
        this->memview_free<memview_policy>(m_tke_old);
        this->memview_free<memview_policy>(m_dzw_stretched);
        this->memview_free<memview_policy>(m_dzt_stretched);
        this->memview_free<memview_policy>(m_tke_old);
        this->memview_free<memview_policy>(m_tke_Av);
        this->memview_free<memview_policy>(m_tke_kv);
        this->memview_free<memview_policy>(m_Nsqr);
        this->memview_free<memview_policy>(m_Ssqr);
        this->memview_free<memview_policy>(m_a_dif);
        this->memview_free<memview_policy>(m_b_dif);
        this->memview_free<memview_policy>(m_c_dif);
        this->memview_free<memview_policy>(m_a_tri);
        this->memview_free<memview_policy>(m_b_tri);
        this->memview_free<memview_policy>(m_c_tri);
        this->memview_free<memview_policy>(m_d_tri);
        this->memview_free<memview_policy>(m_sqrttke);
        this->memview_free<memview_policy>(m_forc);
        this->memview_free<memview_policy>(m_ke);
        this->memview_free<memview_policy>(m_cp);
        this->memview_free<memview_policy>(m_dp);
        this->memview_free<memview_policy>(m_tke_upd);
        this->memview_free<memview_policy>(m_tke_unrest);
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

#endif  // SRC_BACKENDS_TKE_BACKEND_HPP_
