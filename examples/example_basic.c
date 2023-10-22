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
      TKE_Calc(0, nblocks-1, tke, dolic_c);
      
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
