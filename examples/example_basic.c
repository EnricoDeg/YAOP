/* Copyright (C) 2023  Enrico Degregori, Wilton Jaciel Loch
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include "src/TKE.h"

int main(int argc, char ** argv) {
    int nproma = 10240;
    int nlevs = 56;
    int nblocks = 1;
    int ntimesteps = 10;

    int edges_block_size = nblocks;
    int edges_start_block = 0;
    int edges_end_block = nblocks - 1;
    int edges_start_index = 0;
    int edges_end_index = nproma - 1;
    int cells_block_size = nblocks;
    int cells_start_block = 0;
    int cells_end_block = nblocks - 1;
    int cells_start_index = 0;
    int cells_end_index = nproma - 1;

    int vert_mix_type = 2;
    int vmix_idemix_tke = 4;
    int vert_cor_type = 0;
    double dtime = 0.0;
    double OceanReferenceDensity = 1025.022;
    double grav = 9.80665;
    int l_lc = 0;
    double clc = 0.15;
    double ReferencePressureIndbars = 1035.0*grav*1.0e-4;
    double pi = 3.14159265358979323846264338327950288;

    TKE_Init(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
             vert_cor_type, dtime, OceanReferenceDensity, grav,
             l_lc, clc, ReferencePressureIndbars, pi);

    double *depth_CellInterface = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *prism_center_dist_c = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *inv_prism_center_dist_c = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *prism_thick_c = malloc(nproma * nlevs * nblocks * sizeof(double));
    int *dolic_c = malloc(nproma * nblocks * sizeof(int));
    int *dolic_e = malloc(nproma * nblocks * sizeof(int));
    double *zlev_i = malloc(nlevs * sizeof(double));
    double *wet_c = malloc(nproma * nlevs * nblocks * sizeof(double));
    int *edges_cell_idx = malloc(nproma * nlevs * nblocks * sizeof(int));
    int *edges_cell_blk = malloc(nproma * nlevs * nblocks * sizeof(int));

    double *temp = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *salt = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *stretch_c = malloc(nproma * nblocks * sizeof(double));
    double *eta_c = malloc(nproma * nblocks * sizeof(double));

    double * tke = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_plc_in = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *hlc_in = malloc(nproma * nblocks * sizeof(double));
    double *wlc_in = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *u_stokes_in = malloc(nproma * nblocks * sizeof(double));
    double *a_veloc_v = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *a_temp_v = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *a_salt_v = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *iwe_Tdis = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *cvmix_dummy_1 = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *cvmix_dummy_2 = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *cvmix_dummy_3 = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Tbpr = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Tspr = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Tdif = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Tdis = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Twin = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Tiwf = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Tbck = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Ttot = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Lmix = malloc(nproma * nlevs * nblocks * sizeof(double));
    double *tke_Pr = malloc(nproma * nlevs * nblocks * sizeof(double));

    double *stress_xw = malloc(nproma * nblocks * sizeof(double));
    double *stress_yw = malloc(nproma * nblocks * sizeof(double));

    double *fu10 = malloc(nproma * nblocks * sizeof(double));

    double *concsum = malloc(nproma * nblocks * sizeof(double));

    // fill array
    for (int i = 0; i < nblocks; ++i)
        for (int j = 0; j < nlevs; ++j)
            for (int k = 0; k < nproma; ++k)
                tke[k + j * nproma + i * nproma * nlevs] = 1.0 * (k + j * nproma + i * nproma * nlevs);

    for (int i = 0; i < nblocks; ++i)
        for (int k = 0; k < nproma; ++k)
            dolic_c[k + i * nproma] = nlevs;

    #pragma acc enter data copyin(depth_CellInterface[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(prism_center_dist_c[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(inv_prism_center_dist_c[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(prism_thick_c[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(dolic_c[0:nproma*nblocks-1], dolic_e[0:nproma*nblocks-1])
    #pragma acc enter data copyin(zlev_i[0:nlevs-1], wet_c[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(edges_cell_idx[0:nproma*nlevs*nblocks-1], edges_cell_blk[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(tke[0:nproma*nlevs*nblocks-1], tke_plc_in[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(hlc_in[0:nproma*nblocks-1], wlc_in[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(u_stokes_in[0:nproma*nblocks-1], a_veloc_v[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(a_temp_v[0:nproma*nlevs*nblocks-1], a_salt_v[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(iwe_Tdis[0:nproma*nlevs*nblocks-1], cvmix_dummy_1[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(cvmix_dummy_2[0:nproma*nlevs*nblocks-1], cvmix_dummy_3[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(tke_Tbpr[0:nproma*nlevs*nblocks-1], tke_Tspr[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(tke_Tdif[0:nproma*nlevs*nblocks-1], tke_Tdis[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(tke_Twin[0:nproma*nlevs*nblocks-1], tke_Tiwf[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(tke_Tbck[0:nproma*nlevs*nblocks-1], tke_Ttot[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(tke_Lmix[0:nproma*nlevs*nblocks-1], tke_Pr[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(temp[0:nproma*nlevs*nblocks-1], salt[0:nproma*nlevs*nblocks-1])
    #pragma acc enter data copyin(stretch_c[0:nproma*nblocks-1], eta_c[0:nproma*nblocks-1])
    #pragma acc enter data copyin(stress_xw[0:nproma*nblocks-1], stress_yw[0:nproma*nblocks-1])
    #pragma acc enter data copyin(fu10[0:nproma*nblocks-1])
    #pragma acc enter data copyin(concsum[0:nproma*nblocks-1])

    for (int t = 0; t < ntimesteps; t++) {
      #pragma acc host_data use_device(tke, dolic_c)
      TKE_Calc(depth_CellInterface, prism_center_dist_c,
               inv_prism_center_dist_c, prism_thick_c,
               dolic_c, dolic_e, zlev_i, wet_c,
               edges_cell_idx, edges_cell_blk,
               temp, salt, stretch_c, eta_c,
               tke, tke_plc_in, hlc_in, wlc_in,
               u_stokes_in, a_veloc_v, a_temp_v, a_salt_v,
               iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2,
               cvmix_dummy_3, tke_Tbpr, tke_Tspr,
               tke_Tdif, tke_Tdis, tke_Twin,
               tke_Tiwf, tke_Tbck, tke_Ttot,
               tke_Lmix, tke_Pr, stress_xw,
               stress_yw, fu10, concsum,
               edges_block_size, edges_start_block, edges_end_block,
               edges_start_index, edges_end_index, cells_block_size,
               cells_start_block, cells_end_block, cells_start_index,
               cells_end_index);

      #pragma acc wait
    }

    #pragma acc update host(tke[0:nproma*nlevs*nblocks-1])

    for (int i = 0; i < nblocks; ++i)
        for (int j = 0; j < nlevs; ++j)
            for (int k = 0; k < nproma; ++k)
                printf("tke[%d] = %e\n", k + j * nproma + i * nproma * nlevs, tke[k + j * nproma + i * nproma * nlevs]);

    TKE_Finalize();

    #pragma acc exit data delete(depth_CellInterface, prism_center_dist_c, inv_prism_center_dist_c)
    #pragma acc exit data delete(prism_thick_c, dolic_c, dolic_e, zlev_i, wet_c, edges_cell_idx, edges_cell_blk)
    #pragma acc exit data delete(tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v, a_temp_v, a_salt_v)
    #pragma acc exit data delete(iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2, cvmix_dummy_3, tke_Tbpr, tke_Tspr)
    #pragma acc exit data delete(tke_Tdif, tke_Tdis, tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr)
    #pragma acc exit data delete(temp, salt, stretch_c, eta_c)
    #pragma acc exit data delete(stress_xw, stress_yw)
    #pragma acc exit data delete(fu10)
    #pragma acc exit data delete(concsum)

    free(depth_CellInterface);
    free(prism_center_dist_c);
    free(inv_prism_center_dist_c);
    free(prism_thick_c);
    free(dolic_c);
    free(dolic_e);
    free(zlev_i);
    free(wet_c);
    free(edges_cell_idx);
    free(edges_cell_blk);

    free(tke);
    free(tke_plc_in);
    free(hlc_in);
    free(wlc_in);
    free(u_stokes_in);
    free(a_veloc_v);
    free(a_temp_v);
    free(a_salt_v);
    free(iwe_Tdis);
    free(cvmix_dummy_1);
    free(cvmix_dummy_2);
    free(cvmix_dummy_3);
    free(tke_Tbpr);
    free(tke_Tspr);
    free(tke_Tdif);
    free(tke_Tdis);
    free(tke_Twin);
    free(tke_Tiwf);
    free(tke_Tbck);
    free(tke_Ttot);
    free(tke_Lmix);
    free(tke_Pr);

    free(temp);
    free(salt);
    free(stretch_c);
    free(eta_c);

    free(stress_xw);
    free(stress_yw);

    free(fu10);

    free(concsum);

    return 0;
}
