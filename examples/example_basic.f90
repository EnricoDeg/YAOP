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
        CALL TKE_Calc_f(1, nblocks, tke, dolic_c)
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
