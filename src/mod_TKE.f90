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

    subroutine tke_init_f(nproma, nlevs, nblocks, block_size, start_index, end_index)
        implicit none
        integer, intent(in) :: nproma
        integer, intent(in) :: nlevs
        integer, intent(in) :: nblocks
        integer, intent(in) :: block_size
        integer, intent(in) :: start_index
        integer, intent(in) :: end_index

        integer :: start_index_c
        integer :: end_index_c

        interface
            subroutine tke_init_c(nproma, nlevs, nblocks, &
                                  block_size, start_index, end_index) bind(C, name="TKE_Init")
                use iso_c_binding
                implicit none

                integer(c_int), value :: nproma
                integer(c_int), value :: nlevs
                integer(c_int), value :: nblocks
                integer(c_int), value :: block_size
                integer(c_int), value :: start_index
                integer(c_int), value :: end_index
            end subroutine tke_init_c
        end interface

        start_index_c = start_index - 1
        end_index_c = end_index - 1
        
        CALL tke_init_c(nproma, nlevs, nblocks, block_size, start_index_c, end_index_c)
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

    subroutine tke_calc_f(start_block, end_block, tke, dolic_c)
        implicit none
        integer, intent(in) :: start_block
        integer, intent(in) :: end_block
        real(c_double), intent(inout) :: tke(:,:,:)
        integer, intent(in) :: dolic_c(:,:)
        
        type(c_ptr) :: tke_ptr
        type(c_ptr) :: dolic_c_ptr
        integer :: start_block_m1, end_block_m1

        interface
            subroutine tke_calc_c(start_block_c, end_block_c, tke_c, dolic_c_c) bind(C, name="TKE_Calc")
                use iso_c_binding
                implicit none
                integer(c_int), value :: start_block_c
                integer(c_int), value :: end_block_c
                type(c_ptr), value :: tke_c
                type(c_ptr), value :: dolic_c_c
            end subroutine tke_calc_c
        end interface

        start_block_m1 = start_block - 1
        end_block_m1 = end_block - 1
        tke_ptr = c_loc(tke(1,1,1))
        dolic_c_ptr = c_loc(dolic_c(1,1))

        CALL tke_calc_c(start_block_m1, end_block_m1, tke_ptr, dolic_c_ptr)
    end subroutine tke_calc_f

end module mod_TKE
