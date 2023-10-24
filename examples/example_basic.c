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
#include "TKE.h"

int main(int argc, char ** argv) {

    int nproma = 10240;
    int nlevs = 56;
    int nblocks = 1;
    int ntimesteps = 10;

    int block_size = nproma;
    int start_index = 0;
    int end_index = nproma - 1;

    TKE_Init(nproma, nlevs, nblocks, block_size, start_index, end_index);

    double * tke = malloc(nproma * nlevs * nblocks * sizeof(double));
    int *dolic_c = malloc(nproma * nlevs * nblocks * sizeof(int));

    double *depth_CellInterface;
    double *prism_center_dist_c;
    double *inv_prism_center_dist_c;
    double *prism_thick_c;
    int *dolic_e;
    double *zlev_i;
    double *wet_c;
    int *edges_cell_idx;
    int *edges_cell_blk;
    double *temp;
    double *salt;
    double *stretch_c;
    double *eta_c;
    double *tke_plc_in;
    double *hlc_in;
    double *wlc_in;
    double *u_stokes_in;
    double *a_veloc_v;
    double *a_temp_v;
    double *a_salt_v;
    double *iwe_Tdis;
    double *cvmix_dummy_1;
    double *cvmix_dummy_2;
    double *cvmix_dummy_3;
    double *tke_Tbpr;
    double *tke_Tspr;
    double *tke_Tdif;
    double *tke_Tdis;
    double *tke_Twin;
    double *tke_Tiwf;
    double *tke_Tbck;
    double *tke_Ttot;
    double *tke_Lmix;
    double *tke_Pr;
    double *stress_xw;
    double *stress_yw;
    double *fu10;
    double *concsum;

    // fill array
    for (int i=0; i<nblocks; ++i)
        for (int j=0; j<nlevs; ++j)
            for (int k=0; k<nproma; ++k)
                tke[k + j * nproma + i * nproma * nlevs] = 1.0 * (k + j * nproma + i * nproma * nlevs);

    for (int i=0; i<nblocks; ++i)
        for (int k=0; k<nproma; ++k)
            dolic_c[k + i * nproma] = nlevs;

    #pragma acc enter data copyin(tke[0:nproma*nlevs*nblocks-1], dolic_c[0:nproma*nblocks-1])

    for (int t=0; t<ntimesteps; t++) {
      #pragma acc host_data use_device(tke, dolic_c)
      TKE_Calc(0, nblocks-1, depth_CellInterface, prism_center_dist_c,
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
               stress_yw, fu10, concsum);

      #pragma acc wait
    }

    #pragma acc update host(tke[0:nproma*nlevs*nblocks-1])

    for (int i=0; i<nblocks; ++i)
        for (int j=0; j<nlevs; ++j)
            for (int k=0; k<nproma; ++k)
                printf("tke[%d] = %e\n", k + j * nproma + i * nproma * nlevs, tke[k + j * nproma + i * nproma * nlevs]);

    TKE_Finalize();

    #pragma acc exit data delete(tke, dolic_c)
    free(tke);
    free(dolic_c);

    return 0;

}
