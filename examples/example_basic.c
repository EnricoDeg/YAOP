#include <stdio.h>
#include <stdlib.h>
#include "TKE.h"

int main(int argc, char ** argv) {

    int nproma = 6;
    int nlevs = 4;
    int nblocks = 2;
    int block_size = nproma;
    int start_index = 0;
    int end_index = nproma - 1;

    TKE_Init(nproma, nlevs, nblocks, block_size, start_index, end_index);

    double * temperature = malloc(nproma * nlevs * nblocks * sizeof(double));

    // fill array
    for (int i=0; i<nblocks; ++i)
        for (int j=0; j<nlevs; ++j)
            for (int k=0; k<nproma; ++k)
                temperature[k + j * nproma + i * nproma * nlevs] = 1.0 * (k + j * nproma + i * nproma * nlevs);

    TKE_Calc(0,nblocks-1,temperature);

    TKE_Finalize();

    free(temperature);

    return 0;

}
