PROGRAM main
    USE mod_TKE
    IMPLICIT NONE

    INTEGER :: nproma = 8
    INTEGER :: nlevs = 64
    INTEGER :: nblocks = 300

    CALL TKE_Init_f(nproma, nlevs, nblocks)

END PROGRAM