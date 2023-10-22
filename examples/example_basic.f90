program main
    use mod_TKE
    implicit none

    integer, parameter :: pd = 12
    integer, parameter :: rd = 307
    integer, parameter :: dp = selected_real_kind(pd,rd)

    integer :: nproma = 6
    integer :: nlevs = 4
    integer :: nblocks = 2
    integer :: block_size
    integer :: start_index
    integer :: end_index

    integer :: i, j, k

    real(dp), allocatable, dimension(:,:,:) :: temperature

    block_size = nproma
    start_index = 1
    end_index = nproma

    allocate(temperature(nproma, nlevs, nblocks))

    CALL TKE_Init_f(nproma, nlevs, nblocks, block_size, start_index, end_index)

    ! Fill array
    do k=1,nblocks
        do j=1,nlevs
            do i=1,nproma
                temperature(i,j,k) = 1.0_dp * ((i-1) + (j-1) * nproma + (k-1) * nproma * nlevs)
            end do
        end do
    end do

    CALL TKE_Calc_f(1,nblocks,temperature)

    do k=1,nblocks
        do j=1,nlevs
            do i=1,nproma
                write(*,*) "temperature(i,j,k)", temperature(i,j,k)
            end do
        end do
    end do

    CALL TKE_Finalize_f()

    deallocate(temperature)

end program
