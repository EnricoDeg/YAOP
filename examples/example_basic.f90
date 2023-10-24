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

    integer :: block_size
    integer :: start_index
    integer :: end_index

    integer :: i, j, k, t

    real(dp), allocatable, dimension(:,:,:) :: tke
    integer, allocatable, dimension(:,:) :: dolic_c

    real(dp), allocatable :: depth_CellInterface(:,:,:)
    real(dp), allocatable :: prism_center_dist_c(:,:,:)
    real(dp), allocatable :: inv_prism_center_dist_c(:,:,:)
    real(dp), allocatable :: prism_thick_c(:,:,:)
    integer, allocatable :: dolic_e(:,:)
    real(dp), allocatable :: zlev_i(:)
    real(dp), allocatable :: wet_c(:)
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

    block_size = nproma
    start_index = 1
    end_index = nproma

    allocate(tke(nproma, nlevs, nblocks))
    allocate(dolic_c(nproma,nblocks))

    CALL TKE_Init_f(nproma, nlevs, nblocks, block_size, start_index, end_index)

    ! Fill array
    do k=1,nblocks
        do j=1,nlevs
            do i=1,nproma
                tke(i,j,k) = 1.0_dp * ((i-1) + (j-1) * nproma + (k-1) * nproma * nlevs)
            end do
        end do
    end do

    dolic_c(:,:) = nlevs;
    !$ACC ENTER DATA COPYIN(dolic_c, tke)

    do t=1,ntimesteps
        !$ACC HOST_DATA USE_DEVICE(tke, dolic_c)
        CALL TKE_Calc_f(1, nblocks, depth_CellInterface, prism_center_dist_c, &
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
                        stress_yw, fu10, concsum)
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

    !$ACC EXIT DATA DELETE(dolic_c, tke)

    deallocate(dolic_c)
    deallocate(tke)

end program
