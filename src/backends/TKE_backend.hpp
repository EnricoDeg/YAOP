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

#if defined CUDA || defined HIP

#else
#include "src/backends/CPU/cpu_memory.hpp"
#endif

/*! \brief TKE backend class.
 *
 *  It implements the part which is in common for each backend:
 *   - library frontend-backend interface
 *   - library initialization and finalization
 */
template <class T>
class TKE_backend {
 public:
    typedef std::shared_ptr<TKE_backend> Ptr;

    /*! \brief TKE backend class constructor.
    */
    TKE_backend(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                             int vert_cor_type, T dtime, T OceanReferenceDensity, T grav,
                             int l_lc, T clc, T ReferencePressureIndbars, T pi) {
        // Fill structures with parameters
        p_constant.nproma = nproma;
        p_constant.nblocks = nblocks;
        p_constant.vert_mix_type = vert_mix_type;
        p_constant.vmix_idemix_tke = vmix_idemix_tke;
        p_constant.vert_cor_type = vert_cor_type;
        p_constant.dtime = dtime;
        p_constant.OceanReferenceDensity = OceanReferenceDensity;
        p_constant.grav = grav;
        p_constant.l_lc = l_lc;
        p_constant.clc = clc;
        p_constant.ReferencePressureIndbars = ReferencePressureIndbars;
        p_constant.pi = pi;
        p_constant.nlevs = nlevs;

        // Internal parameters are set for now to default values
        p_constant_tke.c_k = 0.1;
        p_constant_tke.c_eps = 0.7;
        p_constant_tke.cd = 3.75;
        p_constant_tke.alpha_tke = 30.0;
        p_constant_tke.clc = 0.15;
        p_constant_tke.mxl_min = 1.0e-8;
        p_constant_tke.KappaM_min = 1.0e-4;
        p_constant_tke.KappaH_min = 1.0e-5;
        p_constant_tke.KappaM_max = 100.0;
        p_constant_tke.tke_surf_min = 1.0e-4;
        p_constant_tke.tke_min = 1.0e-6;
        p_constant_tke.tke_mxl_choice = 2;
        p_constant_tke.handle_old_vals = 1;
        p_constant_tke.only_tke = true;
        p_constant_tke.use_Kappa_min = false;
        p_constant_tke.use_ubound_dirichlet = false;
        p_constant_tke.use_lbound_dirichlet = false;

        m_is_view_init = false;
    }

    /*! \brief TKE backend class destructor.
    */
    virtual ~TKE_backend() = default;

    /*! \brief TKE backend calculation.
    *
    *  It calls calc_impl which contains the actual implementation.
    */
    void calc(t_patch_base *p_patch, t_cvmix_base *p_cvmix,
                           t_ocean_state_base *ocean_state, t_atmo_fluxes_base *atmos_fluxes,
                           t_atmos_for_ocean_base *p_as, t_sea_ice_base *p_sea_ice,
                           int edges_block_size, int edges_start_block, int edges_end_block,
                           int edges_start_index, int edges_end_index, int cells_block_size,
                           int cells_start_block, int cells_end_block, int cells_start_index,
                           int cells_end_index) {
        this->calc_impl(p_patch, p_cvmix, ocean_state, atmos_fluxes, p_as, p_sea_ice,
                        edges_block_size, edges_start_block, edges_end_block,
                        edges_start_index, edges_end_index, cells_block_size,
                        cells_start_block, cells_end_block, cells_start_index,
                        cells_end_index);
    }

 protected:
    /*! \brief Polymorphic function for the actual TKE scheme backend implementation.
    *
    */
    virtual void calc_impl(t_patch_base *p_patch, t_cvmix_base *p_cvmix,
                           t_ocean_state_base *ocean_state, t_atmo_fluxes_base *atmos_fluxes,
                           t_atmos_for_ocean_base *p_as, t_sea_ice_base *p_sea_ice,
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
    template <template <class> class memview_policy>
    void fill_struct_memview(t_patch_base *p_patch, int nblocks, int nlevs, int nproma) {
        T *ptr = static_cast<t_patch<T>*>(p_patch)->get_depth_CellInterface();
        this->p_patch_view.depth_CellInterface = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_patch<T>*>(p_patch)->get_prism_center_dist_c();
        this->p_patch_view.prism_center_dist_c = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_patch<T>*>(p_patch)->get_inv_prism_center_dist_c();
        this->p_patch_view.inv_prism_center_dist_c = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_patch<T>*>(p_patch)->get_prism_thick_c();
        this->p_patch_view.prism_thick_c = memview_policy<T>::memview(ptr, nblocks, nlevs, nproma);
        int *ptr_i = static_cast<t_patch<T>*>(p_patch)->get_dolic_c();
        this->p_patch_view.dolic_c = memview_policy<int>::memview(ptr_i, nblocks, nproma);
        ptr_i = static_cast<t_patch<T>*>(p_patch)->get_dolic_e();
        this->p_patch_view.dolic_e = memview_policy<int>::memview(ptr_i, nblocks, nproma);
        ptr = static_cast<t_patch<T>*>(p_patch)->get_zlev_i();
        this->p_patch_view.zlev_i = memview_policy<T>::memview(ptr, nlevs);
        ptr = static_cast<t_patch<T>*>(p_patch)->get_wet_c();
        this->p_patch_view.wet_c = memview_policy<T>::memview(ptr, nblocks, nlevs, nproma);
        ptr_i = static_cast<t_patch<T>*>(p_patch)->get_edges_cell_idx();
        this->p_patch_view.edges_cell_idx = memview_policy<int>::memview(ptr_i, 2, nlevs, nproma);
        ptr_i = static_cast<t_patch<T>*>(p_patch)->get_edges_cell_blk();
        this->p_patch_view.edges_cell_blk = memview_policy<int>::memview(ptr_i, 2, nlevs, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the sea ice info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class> class memview_policy>
    void fill_struct_memview(t_sea_ice_base *p_sea_ice, int nblocks, int nproma) {
        T *ptr = static_cast<t_sea_ice<T>*>(p_sea_ice)->get_concsum();
        this->p_sea_ice_view.concsum = memview_policy<T>::memview(ptr, nblocks, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the atmosphere info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class> class memview_policy>
    void fill_struct_memview(t_atmos_for_ocean_base *p_as, int nblocks, int nproma) {
        T *ptr = static_cast<t_atmos_for_ocean<T>*>(p_as)->get_fu10();
        this->p_as_view.fu10 = memview_policy<T>::memview(ptr, nblocks, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the atmosphere fluxes info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class> class memview_policy>
    void fill_struct_memview(t_atmo_fluxes_base *atmos_fluxes, int nblocks, int nproma) {
        T *ptr = static_cast<t_atmo_fluxes<T>*>(atmos_fluxes)->get_stress_xw();
        this->atmos_fluxes_view.stress_xw = memview_policy<T>::memview(ptr, nblocks, nproma);
        ptr = static_cast<t_atmo_fluxes<T>*>(atmos_fluxes)->get_stress_yw();
        this->atmos_fluxes_view.stress_yw = memview_policy<T>::memview(ptr, nblocks, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the ocean state info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class> class memview_policy>
    void fill_struct_memview(t_ocean_state_base *ocean_state, int nblocks, int nlevs, int nproma) {
        T *ptr = static_cast<t_ocean_state<T>*>(ocean_state)->get_temp();
        this->ocean_state_view.temp = memview_policy<T>::memview(ptr, nblocks, nlevs, nproma);
        ptr = static_cast<t_ocean_state<T>*>(ocean_state)->get_salt();
        this->ocean_state_view.salt = memview_policy<T>::memview(ptr, nblocks, nlevs, nproma);
        ptr = static_cast<t_ocean_state<T>*>(ocean_state)->get_stretch_c();
        this->ocean_state_view.stretch_c = memview_policy<T>::memview(ptr, nblocks, nproma);
        ptr = static_cast<t_ocean_state<T>*>(ocean_state)->get_eta_c();
        this->ocean_state_view.eta_c = memview_policy<T>::memview(ptr, nblocks, nproma);
        ptr = static_cast<t_ocean_state<T>*>(ocean_state)->get_p_vn_x1();
        this->ocean_state_view.p_vn_x1 = memview_policy<T>::memview(ptr, nblocks, nlevs, nproma);
        ptr = static_cast<t_ocean_state<T>*>(ocean_state)->get_p_vn_x2();
        this->ocean_state_view.p_vn_x2 = memview_policy<T>::memview(ptr, nblocks, nlevs, nproma);
        ptr = static_cast<t_ocean_state<T>*>(ocean_state)->get_p_vn_x3();
        this->ocean_state_view.p_vn_x3 = memview_policy<T>::memview(ptr, nblocks, nlevs, nproma);
    }

    /*! \brief fill a structure of memory views given a structure of pointers about the cvmix info.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class> class memview_policy>
    void fill_struct_memview(t_cvmix_base *p_cvmix, int nblocks, int nlevs, int nproma) {
        T *ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke();
        this->p_cvmix_view.tke = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_plc();
        this->p_cvmix_view.tke_plc = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_hlc();
        this->p_cvmix_view.hlc = memview_policy<T>::memview(ptr, nblocks, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_wlc();
        this->p_cvmix_view.wlc = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_u_stokes();
        this->p_cvmix_view.u_stokes = memview_policy<T>::memview(ptr, nblocks, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_a_veloc_v();
        this->p_cvmix_view.a_veloc_v = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_a_temp_v();
        this->p_cvmix_view.a_temp_v = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_a_salt_v();
        this->p_cvmix_view.a_salt_v = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_iwe_Tdis();
        this->p_cvmix_view.iwe_Tdis = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_cvmix_dummy_1();
        this->p_cvmix_view.cvmix_dummy_1 = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_cvmix_dummy_2();
        this->p_cvmix_view.cvmix_dummy_2 = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_cvmix_dummy_3();
        this->p_cvmix_view.cvmix_dummy_3 = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Tbpr();
        this->p_cvmix_view.tke_Tbpr = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Tspr();
        this->p_cvmix_view.tke_Tspr = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Tdif();
        this->p_cvmix_view.tke_Tdif = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Tdis();
        this->p_cvmix_view.tke_Tdis = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Twin();
        this->p_cvmix_view.tke_Twin = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Tiwf();
        this->p_cvmix_view.tke_Tiwf = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Tbck();
        this->p_cvmix_view.tke_Tbck = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Ttot();
        this->p_cvmix_view.tke_Ttot = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Lmix();
        this->p_cvmix_view.tke_Lmix = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
        ptr = static_cast<t_cvmix<T>*>(p_cvmix)->get_tke_Pr();
        this->p_cvmix_view.tke_Pr = memview_policy<T>::memview(ptr, nblocks, nlevs+1, nproma);
    }

    /*! \brief allocate internal memory and return a 1D memory view object of the allocated memory.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              template <class> class memview_policy>
    memview<T, dext<int, 1>> memview_malloc(T *field, int dim1) {
        return memview_policy<T>::memview_malloc(field, dim1);
    }

    /*! \brief allocate internal memory and return a 2D memory view object of the allocated memory.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              template <class> class memview_policy>
    memview<T, dext<int, 2>> memview_malloc(T *field, int dim1, int dim2) {
        return memview_policy<T>::memview_malloc(field, dim1, dim2);
    }

    /*! \brief allocate internal memory and return a 3D memory view object of the allocated memory.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              template <class> class memview_policy>
    memview<T, dext<int, 3>> memview_malloc(T *field, int dim1, int dim2, int dim3) {
        return memview_policy<T>::memview_malloc(field, dim1, dim2, dim3);
    }

    /*! \brief deallocate internal memory.
    *
    *   It is templated with a memview_policy which defines how to deallocate memory in the actual backend.
    */
    template <template <class> class memview_policy>
    void memview_free(T *field) {
        return memview_policy<T>::memview_free(field);
    }

    /*! \brief fill the internal data structure allocating the arrays and creating memory views.
    *
    *   It is templated with a memview class and a dext class which define the memory view implementation
    *   and with a memview_policy which defines how to allocate and deallocate memory in the actual backend
    *   and how to create a memory view object.
    */
    template <template <class ...> class memview,
              template <class, size_t> class dext,
              template <class> class memview_policy>
    void internal_fields_malloc() {
        this->p_internal_view.tke_old = this->memview_malloc<memview, dext, memview_policy>
                                         (m_tke_old, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.forc_tke_surf_2D = this->memview_malloc<memview, dext, memview_policy>
                                                  (m_forc_tke_surf_2D, p_constant.nproma);
        this->p_internal_view.dzw_stretched = this->memview_malloc<memview, dext, memview_policy>
                                               (m_dzw_stretched, p_constant.nlevs, p_constant.nproma);
        this->p_internal_view.dzt_stretched = this->memview_malloc<memview, dext, memview_policy>
                                               (m_dzt_stretched, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.tke_Av = this->memview_malloc<memview, dext, memview_policy>
                                        (m_tke_Av, p_constant.nblocks, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.tke_kv = this->memview_malloc<memview, dext, memview_policy>
                                        (m_tke_kv, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.Nsqr = this->memview_malloc<memview, dext, memview_policy>
                                      (m_Nsqr, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.Ssqr = this->memview_malloc<memview, dext, memview_policy>
                                      (m_Ssqr, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.a_dif = this->memview_malloc<memview, dext, memview_policy>
                                       (m_a_dif, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.b_dif = this->memview_malloc<memview, dext, memview_policy>
                                       (m_b_dif, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.c_dif = this->memview_malloc<memview, dext, memview_policy>
                                       (m_c_dif, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.a_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_a_tri, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.b_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_b_tri, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.c_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_c_tri, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.d_tri = this->memview_malloc<memview, dext, memview_policy>
                                       (m_d_tri, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.sqrttke = this->memview_malloc<memview, dext, memview_policy>
                                         (m_sqrttke, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.forc = this->memview_malloc<memview, dext, memview_policy>
                                      (m_forc, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.ke = this->memview_malloc<memview, dext, memview_policy>
                                    (m_ke, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.cp = this->memview_malloc<memview, dext, memview_policy>
                                    (m_cp, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.dp = this->memview_malloc<memview, dext, memview_policy>
                                    (m_dp, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.tke_upd = this->memview_malloc<memview, dext, memview_policy>
                                         (m_tke_upd, p_constant.nlevs+1, p_constant.nproma);
        this->p_internal_view.tke_unrest = this->memview_malloc<memview, dext, memview_policy>
                                            (m_tke_unrest, p_constant.nlevs+1, p_constant.nproma);
    }

    /*! \brief free the internal data structure memory deallocating the arrays.
    *
    *   It is templated with a memview_policy which defines how to deallocate memory in the actual backend.
    */
    template <template <class> class memview_policy>
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
    struct t_constant<T> p_constant;
    struct t_constant_tke<T> p_constant_tke;

    bool m_is_view_init;

    T *m_tke_old;
    T *m_tke_Av;
    T *m_tke_kv;
    T *m_forc_tke_surf_2D;
    T *m_dzw_stretched;
    T *m_dzt_stretched;
    T *m_Nsqr;
    T *m_Ssqr;
    T *m_a_dif;
    T *m_b_dif;
    T *m_c_dif;
    T *m_a_tri;
    T *m_b_tri;
    T *m_c_tri;
    T *m_d_tri;
    T *m_sqrttke;
    T *m_forc;
    T *m_ke;
    T *m_cp;
    T *m_dp;
    T *m_tke_upd;
    T *m_tke_unrest;

    // Structures with memory views
    struct t_cvmix_view<T, memview_nms::mdspan, memview_nms::dextents> p_cvmix_view;
    struct t_patch_view<T, memview_nms::mdspan, memview_nms::dextents> p_patch_view;
    struct t_ocean_state_view<T, memview_nms::mdspan, memview_nms::dextents> ocean_state_view;
    struct t_atmo_fluxes_view<T, memview_nms::mdspan, memview_nms::dextents> atmos_fluxes_view;
    struct t_atmos_for_ocean_view<T, memview_nms::mdspan, memview_nms::dextents> p_as_view;
    struct t_sea_ice_view<T, memview_nms::mdspan, memview_nms::dextents> p_sea_ice_view;
    struct t_tke_internal_view<T, memview_nms::mdspan, memview_nms::dextents> p_internal_view;
};

#endif  // SRC_BACKENDS_TKE_BACKEND_HPP_
