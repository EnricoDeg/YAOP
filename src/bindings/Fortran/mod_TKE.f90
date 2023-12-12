! Copyright (C) 2023  Enrico Degregori, Wilton Jaciel Loch
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module mod_TKE
    use iso_c_binding
    implicit none
    private

    integer, parameter :: pd = 12
    integer, parameter :: rd = 307
    integer, parameter :: dp = selected_real_kind(pd,rd)

    public :: tke_init_f
    public :: tke_finalize_f
    public :: tke_calc_f

    contains

    subroutine tke_init_f(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke, &
                          vert_cor_type, dtime, OceanReferenceDensity, grav, &
                          l_lc, clc, ReferencePressureIndbars, pi)
        implicit none
        integer, intent(in)  :: nproma
        integer, intent(in)  :: nlevs
        integer, intent(in)  :: nblocks
        integer, intent(in)  :: vert_mix_type
        integer, intent(in)  :: vmix_idemix_tke
        integer, intent(in)  :: vert_cor_type
        real(dp), intent(in) :: dtime
        real(dp), intent(in) :: OceanReferenceDensity
        real(dp), intent(in) :: grav
        integer, intent(in)  :: l_lc
        real(dp), intent(in) :: clc
        real(dp), intent(in) :: ReferencePressureIndbars
        real(dp), intent(in) :: pi

        interface
            subroutine tke_init_c(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke, &
                                  vert_cor_type, dtime, OceanReferenceDensity, grav, &
                                  l_lc, clc, ReferencePressureIndbars, pi) bind(C, name="TKE_Init")
                use iso_c_binding
                implicit none

                integer(c_int), value :: nproma
                integer(c_int), value :: nlevs
                integer(c_int), value :: nblocks
                integer(c_int), value :: vert_mix_type
                integer(c_int), value :: vmix_idemix_tke
                integer(c_int), value :: vert_cor_type
                real(c_double), value :: dtime
                real(c_double), value :: OceanReferenceDensity
                real(c_double), value :: grav
                integer(c_int), value :: l_lc
                real(c_double), value :: clc
                real(c_double), value :: ReferencePressureIndbars
                real(c_double), value :: pi
            end subroutine tke_init_c
        end interface

        CALL tke_init_c(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke, &
                        vert_cor_type, dtime, OceanReferenceDensity, grav, &
                        l_lc, clc, ReferencePressureIndbars, pi)
    end subroutine tke_init_f

    subroutine tke_finalize_f()
        implicit none

        interface
            subroutine tke_finalize_c() bind(C, name="TKE_Finalize")
                use iso_c_binding
                implicit none

            end subroutine tke_finalize_c
        end interface

        CALL tke_finalize_c()
    end subroutine tke_finalize_f

    subroutine tke_calc_f(depth_CellInterface, prism_center_dist_c, &
                          inv_prism_center_dist_c, prism_thick_c, &
                          dolic_c, dolic_e, zlev_i, wet_c, &
                          edges_cell_idx, edges_cell_blk, &
                          temp, salt, stretch_c, eta_c, &
                          p_vn_x1, p_vn_x2, p_vn_x3, &
                          tke, tke_plc_in, hlc_in, wlc_in, &
                          u_stokes_in, a_veloc_v, a_temp_v, a_salt_v, &
                          iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2, &
                          cvmix_dummy_3, tke_Tbpr, tke_Tspr, &
                          tke_Tdif, tke_Tdis, tke_Twin, &
                          tke_Tiwf, tke_Tbck, tke_Ttot, &
                          tke_Lmix, tke_Pr, stress_xw, &
                          stress_yw, fu10, concsum, &
                          edges_block_size, edges_start_block, edges_end_block, &
                          edges_start_index, edges_end_index, cells_block_size, &
                          cells_start_block, cells_end_block, cells_start_index, &
                          cells_end_index)

        implicit none


        real(c_double), intent(in),    TARGET :: depth_CellInterface(:,:,:)
        real(c_double), intent(in),    TARGET :: prism_center_dist_c(:,:,:)
        real(c_double), intent(in),    TARGET :: inv_prism_center_dist_c(:,:,:)
        real(c_double), intent(in),    TARGET :: prism_thick_c( :,:,:)
        integer, intent(in),           TARGET :: dolic_c(:,:)
        integer, intent(in),           TARGET :: dolic_e(:,:)
        real(c_double), intent(in),    TARGET :: zlev_i(:)
        real(c_double), intent(in),    TARGET :: wet_c(:,:,:)
        integer, intent(in),           TARGET :: edges_cell_idx(:,:,:)
        integer, intent(in),           TARGET :: edges_cell_blk(:,:,:)
        real(c_double), intent(in),    TARGET :: temp(:,:,:)
        real(c_double), intent(in),    TARGET :: salt(:,:,:)
        real(c_double), intent(in),    TARGET :: stretch_c(:,:)
        real(c_double), intent(in),    TARGET :: eta_c(:,:)
        real(c_double), intent(in),    TARGET :: p_vn_x1(:,:,:)
        real(c_double), intent(in),    TARGET :: p_vn_x2(:,:,:)
        real(c_double), intent(in),    TARGET :: p_vn_x3(:,:,:)
        real(c_double), intent(inout), TARGET :: tke(:,:,:)
        real(c_double), intent(in),    TARGET :: tke_plc_in(:,:,:)
        real(c_double), intent(in),    TARGET :: hlc_in(:,:)
        real(c_double), intent(in),    TARGET :: wlc_in(:,:,:)
        real(c_double), intent(in),    TARGET :: u_stokes_in(:,:)
        real(c_double), intent(inout), TARGET :: a_veloc_v(:,:,:)
        real(c_double), intent(inout), TARGET :: a_temp_v(:,:,:)
        real(c_double), intent(inout), TARGET :: a_salt_v(:,:,:)
        real(c_double), intent(in),    TARGET :: iwe_Tdis(:,:,:)
        real(c_double), intent(inout), TARGET :: cvmix_dummy_1(:,:,:)
        real(c_double), intent(inout), TARGET :: cvmix_dummy_2(:,:,:)
        real(c_double), intent(inout), TARGET :: cvmix_dummy_3(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Tbpr(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Tspr(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Tdif(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Tdis(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Twin(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Tiwf(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Tbck(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Ttot(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Lmix(:,:,:)
        real(c_double), intent(inout), TARGET :: tke_Pr(:,:,:)
        real(c_double), intent(in),    TARGET :: stress_xw(:,:)
        real(c_double), intent(in),    TARGET :: stress_yw(:,:)
        real(c_double), intent(in),    TARGET :: fu10(:,:)
        real(c_double), intent(in),    TARGET :: concsum(:,:)
        integer, intent(in)                    :: edges_block_size
        integer, intent(in)                    :: edges_start_block
        integer, intent(in)                    :: edges_end_block
        integer, intent(in)                    :: edges_start_index
        integer, intent(in)                    :: edges_end_index
        integer, intent(in)                    :: cells_block_size
        integer, intent(in)                    :: cells_start_block
        integer, intent(in)                    :: cells_end_block
        integer, intent(in)                    :: cells_start_index
        integer, intent(in)                    :: cells_end_index

        type(c_ptr) :: depth_CellInterface_ptr, prism_center_dist_c_ptr, &
                       inv_prism_center_dist_c_ptr, prism_thick_c_ptr, &
                       dolic_c_ptr, dolic_e_ptr, zlev_i_ptr, wet_c_ptr, &
                       edges_cell_idx_ptr, edges_cell_blk_ptr, &
                       temp_ptr, salt_ptr, stretch_c_ptr, eta_c_ptr, &
                       p_vn_x1_ptr, p_vn_x2_ptr, p_vn_x3_ptr, &
                       tke_ptr, tke_plc_in_ptr, hlc_in_ptr, wlc_in_ptr, &
                       u_stokes_in_ptr, a_veloc_v_ptr, a_temp_v_ptr, a_salt_v_ptr, &
                       iwe_Tdis_ptr, cvmix_dummy_1_ptr, cvmix_dummy_2_ptr, &
                       cvmix_dummy_3_ptr, tke_Tbpr_ptr, tke_Tspr_ptr, &
                       tke_Tdif_ptr, tke_Tdis_ptr, tke_Twin_ptr, &
                       tke_Tiwf_ptr, tke_Tbck_ptr, tke_Ttot_ptr, &
                       tke_Lmix_ptr, tke_Pr_ptr, stress_xw_ptr, &
                       stress_yw_ptr, fu10_ptr, concsum_ptr

        interface
            subroutine tke_calc_c(depth_CellInterface_c, prism_center_dist_c_c, &
                                  inv_prism_center_dist_c_c, prism_thick_c_c, &
                                  dolic_c_c, dolic_e_c, zlev_i_c, wet_c_c, &
                                  edges_cell_idx_c, edges_cell_blk_c, &
                                  temp_c, salt_c, stretch_c_c, eta_c_c, &
                                  p_vn_x1_c, p_vn_x2_c, p_vn_x3_c, &
                                  tke_c, tke_plc_in_c, hlc_in_c, wlc_in_c, &
                                  u_stokes_in_c, a_veloc_v_c, a_temp_v_c, a_salt_v_c, &
                                  iwe_Tdis_c, cvmix_dummy_1_c, cvmix_dummy_2_c, &
                                  cvmix_dummy_3_c, tke_Tbpr_c, tke_Tspr_c, &
                                  tke_Tdif_c, tke_Tdis_c, tke_Twin_c, &
                                  tke_Tiwf_c, tke_Tbck_c, tke_Ttot_c, &
                                  tke_Lmix_c, tke_Pr_c, stress_xw_c, &
                                  stress_yw_c, fu10_c, concsum_c, &
                                  edges_block_size_c, edges_start_block_c, edges_end_block_c, &
                                  edges_start_index_c, edges_end_index_c, cells_block_size_c, &
                                  cells_start_block_c, cells_end_block_c, cells_start_index_c, &
                                  cells_end_index_c) bind(C, name="TKE_Calc")
                use iso_c_binding

                implicit none
                type(c_ptr), value    :: depth_CellInterface_c
                type(c_ptr), value    :: prism_center_dist_c_c
                type(c_ptr), value    :: inv_prism_center_dist_c_c
                type(c_ptr), value    :: prism_thick_c_c
                type(c_ptr), value    :: dolic_c_c
                type(c_ptr), value    :: dolic_e_c
                type(c_ptr), value    :: zlev_i_c
                type(c_ptr), value    :: wet_c_c
                type(c_ptr), value    :: edges_cell_idx_c
                type(c_ptr), value    :: edges_cell_blk_c
                type(c_ptr), value    :: temp_c
                type(c_ptr), value    :: salt_c
                type(c_ptr), value    :: stretch_c_c
                type(c_ptr), value    :: eta_c_c
                type(c_ptr), value    :: p_vn_x1_c
                type(c_ptr), value    :: p_vn_x2_c
                type(c_ptr), value    :: p_vn_x3_c
                type(c_ptr), value    :: tke_c
                type(c_ptr), value    :: tke_plc_in_c
                type(c_ptr), value    :: hlc_in_c
                type(c_ptr), value    :: wlc_in_c
                type(c_ptr), value    :: u_stokes_in_c
                type(c_ptr), value    :: a_veloc_v_c
                type(c_ptr), value    :: a_temp_v_c
                type(c_ptr), value    :: a_salt_v_c
                type(c_ptr), value    :: iwe_Tdis_c
                type(c_ptr), value    :: cvmix_dummy_1_c
                type(c_ptr), value    :: cvmix_dummy_2_c
                type(c_ptr), value    :: cvmix_dummy_3_c
                type(c_ptr), value    :: tke_Tbpr_c
                type(c_ptr), value    :: tke_Tspr_c
                type(c_ptr), value    :: tke_Tdif_c
                type(c_ptr), value    :: tke_Tdis_c
                type(c_ptr), value    :: tke_Twin_c
                type(c_ptr), value    :: tke_Tiwf_c
                type(c_ptr), value    :: tke_Tbck_c
                type(c_ptr), value    :: tke_Ttot_c
                type(c_ptr), value    :: tke_Lmix_c
                type(c_ptr), value    :: tke_Pr_c
                type(c_ptr), value    :: stress_xw_c
                type(c_ptr), value    :: stress_yw_c
                type(c_ptr), value    :: fu10_c
                type(c_ptr), value    :: concsum_c
                integer(c_int), value :: edges_block_size_c
                integer(c_int), value :: edges_start_block_c
                integer(c_int), value :: edges_end_block_c
                integer(c_int), value :: edges_start_index_c
                integer(c_int), value :: edges_end_index_c
                integer(c_int), value :: cells_block_size_c
                integer(c_int), value :: cells_start_block_c
                integer(c_int), value :: cells_end_block_c
                integer(c_int), value :: cells_start_index_c
                integer(c_int), value :: cells_end_index_c

            end subroutine tke_calc_c
        end interface

        depth_CellInterface_ptr     = c_loc(depth_CellInterface(1, 1, 1))
        prism_center_dist_c_ptr     = c_loc(prism_center_dist_c(1, 1, 1))
        inv_prism_center_dist_c_ptr = c_loc(inv_prism_center_dist_c(1, 1, 1))
        prism_thick_c_ptr           = c_loc(prism_thick_c(1, 1, 1))
        dolic_c_ptr                 = c_loc(dolic_c(1, 1))
        dolic_e_ptr                 = c_loc(dolic_e(1, 1))
        zlev_i_ptr                  = c_loc(zlev_i(1))
        wet_c_ptr                   = c_loc(wet_c(1, 1, 1))
        edges_cell_idx_ptr          = c_loc(edges_cell_idx(1, 1, 1))
        edges_cell_blk_ptr          = c_loc(edges_cell_blk(1, 1, 1))
        temp_ptr                    = c_loc(temp(1, 1, 1))
        salt_ptr                    = c_loc(salt(1, 1, 1))
        stretch_c_ptr               = c_loc(stretch_c(1, 1))
        eta_c_ptr                   = c_loc(eta_c(1, 1))
        p_vn_x1_ptr                 = c_loc(p_vn_x1(1, 1, 1))
        p_vn_x2_ptr                 = c_loc(p_vn_x2(1, 1, 1))
        p_vn_x3_ptr                 = c_loc(p_vn_x3(1, 1, 1))
        tke_ptr                     = c_loc(tke(1, 1, 1))
        tke_plc_in_ptr              = c_loc(tke_plc_in(1, 1, 1))
        hlc_in_ptr                  = c_loc(hlc_in(1, 1))
        wlc_in_ptr                  = c_loc(wlc_in(1, 1, 1))
        u_stokes_in_ptr             = c_loc(u_stokes_in(1, 1))
        a_veloc_v_ptr               = c_loc(a_veloc_v(1, 1, 1))
        a_temp_v_ptr                = c_loc(a_temp_v(1, 1, 1))
        a_salt_v_ptr                = c_loc(a_salt_v(1, 1, 1))
        iwe_Tdis_ptr                = c_loc(iwe_Tdis(1, 1, 1))
        cvmix_dummy_1_ptr           = c_loc(cvmix_dummy_1(1, 1, 1))
        cvmix_dummy_2_ptr           = c_loc(cvmix_dummy_2(1, 1, 1))
        cvmix_dummy_3_ptr           = c_loc(cvmix_dummy_3(1, 1, 1))
        tke_Tbpr_ptr                = c_loc(tke_Tbpr(1, 1, 1))
        tke_Tspr_ptr                = c_loc(tke_Tspr(1, 1, 1))
        tke_Tdif_ptr                = c_loc(tke_Tdif(1, 1, 1))
        tke_Tdis_ptr                = c_loc(tke_Tdis(1, 1, 1))
        tke_Twin_ptr                = c_loc(tke_Twin(1, 1, 1))
        tke_Tiwf_ptr                = c_loc(tke_Tiwf(1, 1, 1))
        tke_Tbck_ptr                = c_loc(tke_Tbck(1, 1, 1))
        tke_Ttot_ptr                = c_loc(tke_Ttot(1, 1, 1))
        tke_Lmix_ptr                = c_loc(tke_Lmix(1, 1, 1))
        tke_Pr_ptr                  = c_loc(tke_Pr(1, 1, 1))
        stress_xw_ptr               = c_loc(stress_xw(1, 1))
        stress_yw_ptr               = c_loc(stress_yw(1, 1))
        fu10_ptr                    = c_loc(fu10(1, 1))
        concsum_ptr                 = c_loc(concsum(1, 1))

        CALL tke_calc_c(depth_CellInterface_ptr, prism_center_dist_c_ptr, &
                        inv_prism_center_dist_c_ptr, prism_thick_c_ptr, &
                        dolic_c_ptr, dolic_e_ptr, zlev_i_ptr, wet_c_ptr, &
                        edges_cell_idx_ptr, edges_cell_blk_ptr, &
                        temp_ptr, salt_ptr, stretch_c_ptr, eta_c_ptr, &
                        p_vn_x1_ptr, p_vn_x2_ptr, p_vn_x3_ptr, &
                        tke_ptr, tke_plc_in_ptr, hlc_in_ptr, wlc_in_ptr, &
                        u_stokes_in_ptr, a_veloc_v_ptr, a_temp_v_ptr, a_salt_v_ptr, &
                        iwe_Tdis_ptr, cvmix_dummy_1_ptr, cvmix_dummy_2_ptr, &
                        cvmix_dummy_3_ptr, tke_Tbpr_ptr, tke_Tspr_ptr, &
                        tke_Tdif_ptr, tke_Tdis_ptr, tke_Twin_ptr, &
                        tke_Tiwf_ptr, tke_Tbck_ptr, tke_Ttot_ptr, &
                        tke_Lmix_ptr, tke_Pr_ptr, stress_xw_ptr, &
                        stress_yw_ptr, fu10_ptr, concsum_ptr, &
                        edges_block_size, edges_start_block-1, edges_end_block-1, &
                        edges_start_index-1, edges_end_index-1, cells_block_size, &
                        cells_start_block-1, cells_end_block-1, cells_start_index-1, &
                        cells_end_index-1)

    end subroutine tke_calc_f

end module mod_TKE
