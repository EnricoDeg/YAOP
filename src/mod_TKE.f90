module mod_TKE
    use iso_c_binding
    implicit none
    private
    public :: tke_init_f

    contains

    subroutine tke_init_f()
        use iso_c_binding
        implicit none
        interface
            subroutine tke_init_c() bind(C, name="TKE_Init")
            
            end subroutine tke_init_c
        end interface
        CALL tke_init_c()
    end subroutine tke_init_f

    

end module mod_TKE