module mod_TKE
    use iso_c_binding
    implicit none
    private
    public :: tke_init_f
    public :: tke_finalize_f

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

end module mod_TKE
