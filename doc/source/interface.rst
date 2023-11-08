.. _interface:

Interface
=========

This section describes the interface between ICON and the TKE library.

ICON
````
The call to the TKE in ICON is in ``mo_ocean_physics`` module and it is passing several derived data types::

    CALL calc_tke(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum)

In order to have a TKE library external to the ICON code, the derived data types need to be expanded, so the ``calc_tke`` subroutine represents an interface to the TKE library on the ICON side::

   SUBROUTINE calc_tke_vector(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum)
      TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
      TYPE(t_hydro_ocean_state), TARGET, INTENT(IN) :: ocean_state
      TYPE(t_atmos_fluxes), INTENT(IN) :: atmos_fluxes
      TYPE(t_ho_params), INTENT(INOUT) :: params_oce
      REAL(wp), TARGET, INTENT(IN) :: fu10(:,:) ! t_atmos_for_ocean%fu10
      REAL(wp), TARGET, INTENT(IN) :: concsum(:,:) !< sea ice concentration
      LOGICAL, INTENT(IN), OPTIONAL        :: use_acc

      CALL calc_lib_tke_vector(patch_3d%p_patch_1d(1)%depth_CellInterface, patch_3d%p_patch_1d(1)%prism_center_dist_c,  & 
                               patch_3d%p_patch_1d(1)%inv_prism_center_dist_c, patch_3d%p_patch_1d(1)%prism_thick_c, &
                               patch_3d%p_patch_1d(1)%dolic_c, patch_3d%p_patch_1d(1)%dolic_e, &
                               patch_3d%p_patch_1d(1)%zlev_i, patch_3d%wet_c, &
                               ocean_state%p_prog(nold(1))%tracer, ocean_state%p_prog(nold(1))%stretch_c, &
                               ocean_state%p_prog(nold(1))%eta_c, ocean_state%p_diag%p_vn, &
                               params_oce%cvmix_params%tke, params_oce%cvmix_params%tke_plc, &
                               params_oce%cvmix_params%hlc, params_oce%cvmix_params%wlc, params_oce%cvmix_params%u_stokes,  &
                               params_oce%a_veloc_v, params_oce%a_tracer_v, params_oce%cvmix_params%iwe_Tdis, &
                               params_oce%cvmix_params%cvmix_dummy_1, params_oce%cvmix_params%cvmix_dummy_2, &
                               params_oce%cvmix_params%cvmix_dummy_3, params_oce%cvmix_params%tke_Tbpr, &
                               params_oce%cvmix_params%tke_Tspr, params_oce%cvmix_params%tke_Tdif, params_oce%cvmix_params%tke_Tdis, &
                               params_oce%cvmix_params%tke_Twin, params_oce%cvmix_params%tke_Tiwf, params_oce%cvmix_params%tke_Tbck, &
                               params_oce%cvmix_params%tke_Ttot, params_oce%cvmix_params%tke_Lmix, params_oce%cvmix_params%tke_Pr, &
                               atmos_fluxes%stress_xw, atmos_fluxes%stress_yw, fu10, concsum,  & 
                               patch_3d%p_patch_2d(1)%edges%cell_idx, patch_3d%p_patch_2d(1)%edges%cell_blk, &
                               patch_3d%p_patch_2d(1)%edges%in_domain%block_size, &
                               patch_3d%p_patch_2d(1)%edges%in_domain%start_block, &
                               patch_3d%p_patch_2d(1)%edges%in_domain%end_block, &
                               patch_3d%p_patch_2d(1)%edges%in_domain%start_index, &
                               patch_3d%p_patch_2d(1)%edges%in_domain%end_index, &
                               patch_3d%p_patch_2d(1)%cells%All%block_size, & 
                               patch_3d%p_patch_2d(1)%cells%All%start_block, &
                               patch_3d%p_patch_2d(1)%cells%All%end_block, &
                               patch_3d%p_patch_2d(1)%cells%All%start_index, &
                               patch_3d%p_patch_2d(1)%cells%All%end_index, &
                               n_zlev, nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks, &
                               vert_mix_type, vmix_idemix_tke, vert_cor_type, &
                               dtime, OceanReferenceDensity, grav, l_lc, clc, ReferencePressureIndbars, pi)

   END SUBROUTINE calc_tke_vector



.. toctree::
   :maxdepth: 2

   cpp_doxygen_sphinx
