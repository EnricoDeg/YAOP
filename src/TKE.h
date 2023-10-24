#ifdef __cplusplus
extern "C" {
#endif


// Constructor
void TKE_Init(int nproma, int nlevs, int nblocks,
              int block_size, int start_index, int end_index);

// Destructor
void TKE_Finalize();

// Calculation
void TKE_Calc(int start_block, int end_block,
               double *depth_CellInterface, double *prism_center_dist_c,
               double *inv_prism_center_dist_c, double *prism_thick_c,
               int *dolic_c, int *dolic_e, double *zlev_i, double *wet_c,
               int *edges_cell_idx, int *edges_cell_blk,
               double *temp, double *salt, double *stretch_c, double *eta_c,
               double *tke, double *tke_plc_in, double *hlc_in, double *wlc_in,
               double *u_stokes_in, double *a_veloc_v, double *a_temp_v, double *a_salt_v,
               double *iwe_Tdis, double *cvmix_dummy_1, double *cvmix_dummy_2,
               double *cvmix_dummy_3, double *tke_Tbpr, double *tke_Tspr,
               double *tke_Tdif, double *tke_Tdis, double *tke_Twin,
               double *tke_Tiwf, double *tke_Tbck, double *tke_Ttot,
               double *tke_Lmix, double *tke_Pr, double *stress_xw,
               double *stress_yw, double *fu10, double *concsum);


#ifdef __cplusplus
}
#endif
