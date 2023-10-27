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

program main
    use mod_TKE
    implicit none

    integer, parameter :: pd = 12
    integer, parameter :: rd = 307
    integer, parameter :: dp = selected_real_kind(pd,rd)

    integer :: nproma = 10240
    integer :: nlevs = 56
    integer :: nblocks = 1
    integer :: ntimesteps = 10

    integer :: edges_block_size
    integer :: edges_start_block
    integer :: edges_end_block
    integer :: edges_start_index
    integer :: edges_end_index
    integer :: cells_block_size
    integer :: cells_start_block
    integer :: cells_end_block
    integer :: cells_start_index
    integer :: cells_end_index

    integer  :: vert_mix_type
    integer  :: vmix_idemix_tke
    integer  :: vert_cor_type
    real(dp) :: dtime
    real(dp) :: OceanReferenceDensity
    real(dp) :: grav
    integer  :: l_lc
    real(dp) :: clc
    real(dp) :: ReferencePressureIndbars
    real(dp) :: pi;

    integer :: i, j, k, t

    real(dp), allocatable, dimension(:,:,:) :: tke
    integer, allocatable, dimension(:,:) :: dolic_c

    real(dp), allocatable :: depth_CellInterface(:,:,:)
    real(dp), allocatable :: prism_center_dist_c(:,:,:)
    real(dp), allocatable :: inv_prism_center_dist_c(:,:,:)
    real(dp), allocatable :: prism_thick_c(:,:,:)
    integer, allocatable :: dolic_e(:,:)
    real(dp), allocatable :: zlev_i(:)
    real(dp), allocatable :: wet_c(:,:,:)
    integer, allocatable :: edges_cell_idx(:,:,:)
    integer, allocatable :: edges_cell_blk(:,:,:)
    real(dp), allocatable :: temp(:,:,:)
    real(dp), allocatable :: salt(:,:,:)
    real(dp), allocatable :: stretch_c(:,:)
    real(dp), allocatable :: eta_c(:,:)
    real(dp), allocatable :: tke_plc_in(:,:,:)
    real(dp), allocatable :: hlc_in(:,:)
    real(dp), allocatable :: wlc_in(:,:,:)
    real(dp), allocatable :: u_stokes_in(:,:)
    real(dp), allocatable :: a_veloc_v(:,:,:)
    real(dp), allocatable :: a_temp_v(:,:,:)
    real(dp), allocatable :: a_salt_v(:,:,:)
    real(dp), allocatable :: iwe_Tdis(:,:,:)
    real(dp), allocatable :: cvmix_dummy_1(:,:,:)
    real(dp), allocatable :: cvmix_dummy_2(:,:,:)
    real(dp), allocatable :: cvmix_dummy_3(:,:,:)
    real(dp), allocatable :: tke_Tbpr(:,:,:)
    real(dp), allocatable :: tke_Tspr(:,:,:)
    real(dp), allocatable :: tke_Tdif(:,:,:)
    real(dp), allocatable :: tke_Tdis(:,:,:)
    real(dp), allocatable :: tke_Twin(:,:,:)
    real(dp), allocatable :: tke_Tiwf(:,:,:)
    real(dp), allocatable :: tke_Tbck(:,:,:)
    real(dp), allocatable :: tke_Ttot(:,:,:)
    real(dp), allocatable :: tke_Lmix(:,:,:)
    real(dp), allocatable :: tke_Pr(:,:,:)
    real(dp), allocatable :: stress_xw(:,:)
    real(dp), allocatable :: stress_yw(:,:)
    real(dp), allocatable :: fu10(:,:)
    real(dp), allocatable :: concsum(:,:)

    edges_block_size = nblocks
    edges_start_block = 1
    edges_end_block = nblocks
    edges_start_index = 1
    edges_end_index = nproma
    cells_block_size = nblocks
    cells_start_block = 1
    cells_end_block = nblocks
    cells_start_index = 1
    cells_end_index = nproma

    vert_mix_type = 2
    vmix_idemix_tke = 4
    vert_cor_type = 0
    dtime = 0.0
    OceanReferenceDensity = 1025.022
    grav = 9.80665
    l_lc = 0
    clc = 0.15
    ReferencePressureIndbars = 1035.0*grav*1.0e-4
    pi = 3.14159265358979323846264338327950288

    allocate(depth_CellInterface(nproma,nlevs,nblocks))
    allocate(prism_center_dist_c(nproma,nlevs,nblocks))
    allocate(inv_prism_center_dist_c(nproma,nlevs,nblocks))
    allocate(prism_thick_c(nproma,nlevs,nblocks))
    allocate(dolic_c(nproma,nblocks))
    allocate(dolic_e(nproma,nblocks))
    allocate(zlev_i(nlevs))
    allocate(wet_c(nproma,nlevs,nblocks))
    allocate(edges_cell_idx(nproma,nlevs,nblocks))
    allocate(edges_cell_blk(nproma,nlevs,nblocks))

    allocate(tke(nproma, nlevs, nblocks))
    allocate(tke_plc_in(nproma, nlevs, nblocks))
    allocate(hlc_in(nproma, nblocks))
    allocate(wlc_in(nproma, nlevs, nblocks))
    allocate(u_stokes_in(nproma, nblocks))
    allocate(a_veloc_v(nproma, nlevs, nblocks))
    allocate(a_temp_v(nproma, nlevs, nblocks))
    allocate(a_salt_v(nproma, nlevs, nblocks))
    allocate(iwe_Tdis(nproma, nlevs, nblocks))
    allocate(cvmix_dummy_1(nproma, nlevs, nblocks))
    allocate(cvmix_dummy_2(nproma, nlevs, nblocks))
    allocate(cvmix_dummy_3(nproma, nlevs, nblocks))
    allocate(tke_Tbpr(nproma, nlevs, nblocks))
    allocate(tke_Tspr(nproma, nlevs, nblocks))
    allocate(tke_Tdif(nproma, nlevs, nblocks))
    allocate(tke_Tdis(nproma, nlevs, nblocks))
    allocate(tke_Twin(nproma, nlevs, nblocks))
    allocate(tke_Tiwf(nproma, nlevs, nblocks))
    allocate(tke_Tbck(nproma, nlevs, nblocks))
    allocate(tke_Ttot(nproma, nlevs, nblocks))
    allocate(tke_Lmix(nproma, nlevs, nblocks))
    allocate(tke_Pr(nproma, nlevs, nblocks))

    allocate(temp(nproma, nlevs, nblocks))
    allocate(salt(nproma, nlevs, nblocks))
    allocate(stretch_c(nproma, nblocks))
    allocate(eta_c(nproma, nblocks))

    allocate(stress_xw(nproma, nblocks))
    allocate(stress_yw(nproma, nblocks))

    allocate(fu10(nproma, nblocks))

    allocate(concsum(nproma, nblocks))

    CALL TKE_Init_f(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke, &
                    vert_cor_type, dtime, OceanReferenceDensity, grav, &
                    l_lc, clc, ReferencePressureIndbars, pi)

    ! Fill array
    do k=1,nblocks
        do j=1,nlevs
            do i=1,nproma
                tke(i,j,k) = 1.0_dp * ((i-1) + (j-1) * nproma + (k-1) * nproma * nlevs)
            end do
        end do
    end do

    dolic_c(:,:) = nlevs;
    !$ACC ENTER DATA COPYIN(depth_CellInterface, prism_center_dist_c, inv_prism_center_dist_c, prism_thick_c, &
    !$ACC                   dolic_c, dolic_e, zlev_i, wet_c, edges_cell_idx, edges_cell_blk)
    !$ACC ENTER DATA COPYIN(tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v, a_temp_v, a_salt_v, iwe_Tdis, &
    !$ACC                   cvmix_dummy_1, cvmix_dummy_2, cvmix_dummy_3, tke_Tbpr, tke_Tspr, tke_Tdif, tke_Tdis, &
    !$ACC                   tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr)
    !$ACC ENTER DATA COPYIN(temp, salt, stretch_c, eta_c)
    !$ACC ENTER DATA COPYIN(stress_xw, stress_yw)
    !$ACC ENTER DATA COPYIN(fu10)
    !$ACC ENTER DATA COPYIN(concsum)

    do t=1,ntimesteps
        !$ACC HOST_DATA USE_DEVICE(tke, dolic_c)
        CALL TKE_Calc_f(depth_CellInterface, prism_center_dist_c, &
                        inv_prism_center_dist_c, prism_thick_c, &
                        dolic_c, dolic_e, zlev_i, wet_c, &
                        edges_cell_idx, edges_cell_blk, &
                        temp, salt, stretch_c, eta_c, &
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
        !$ACC END HOST_DATA
        !$ACC WAIT
    end do

    !$ACC UPDATE HOST(tke)

    do k=1,nblocks
        do j=1,nlevs
            do i=1,nproma
                write(*,*) "tke(i,j,k)", tke(i,j,k)
            end do
        end do
    end do

    CALL TKE_Finalize_f()

    !$ACC EXIT DATA DELETE(depth_CellInterface, prism_center_dist_c, inv_prism_center_dist_c, prism_thick_c, &
    !$ACC                  dolic_c, dolic_e, zlev_i, wet_c, edges_cell_idx, edges_cell_blk)
    !$ACC EXIT DATA DELETE(tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v, a_temp_v, a_salt_v, iwe_Tdis, &
    !$ACC                  cvmix_dummy_1, cvmix_dummy_2, cvmix_dummy_3, tke_Tbpr, tke_Tspr, tke_Tdif, tke_Tdis, &
    !$ACC                  tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr)
    !$ACC EXIT DATA DELETE(temp, salt, stretch_c, eta_c)
    !$ACC EXIT DATA DELETE(stress_xw, stress_yw)
    !$ACC EXIT DATA DELETE(fu10)
    !$ACC EXIT DATA DELETE(concsum)

    deallocate(depth_CellInterface, prism_center_dist_c, inv_prism_center_dist_c, prism_thick_c, &
               dolic_c, dolic_e, zlev_i, wet_c, edges_cell_idx, edges_cell_blk)
    deallocate(tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v, a_temp_v, a_salt_v, iwe_Tdis, &
               cvmix_dummy_1, cvmix_dummy_2, cvmix_dummy_3, tke_Tbpr, tke_Tspr, tke_Tdif, tke_Tdis, &
               tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr)
    deallocate(temp, salt, stretch_c, eta_c)
    deallocate(stress_xw, stress_yw)
    deallocate(fu10)
    deallocate(concsum)

end program
