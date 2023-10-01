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

    subroutine tke_init_f(nproma, nlevs, nblocks)
        implicit none
        integer, intent(in) :: nproma
        integer, intent(in) :: nlevs
        integer, intent(in) :: nblocks

        interface
            subroutine tke_init_c(nproma, nlevs, nblocks) bind(C, name="TKE_Init")
                use iso_c_binding
                implicit none

                integer(c_int), value :: nproma
                integer(c_int), value :: nlevs
                integer(c_int), value :: nblocks
            end subroutine tke_init_c
        end interface
        
        CALL tke_init_c(nproma, nlevs, nblocks)
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

    subroutine tke_calc_f(temperature)
        implicit none
        real(c_double), intent(inout) :: temperature(:,:,:)
        
        type(c_ptr) :: temperature_ptr

        interface
            subroutine tke_calc_c(temperature_c) bind(C, name="TKE_Calc")
                use iso_c_binding
                implicit none
                
                type(c_ptr), value :: temperature_c
            end subroutine tke_calc_c
        end interface

        temperature_ptr = c_loc(temperature(1,1,1))

        CALL tke_calc_c(temperature_ptr)
    end subroutine tke_calc_f

end module mod_TKE
