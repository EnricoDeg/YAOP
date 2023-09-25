#include <stdio.h>
#include "TKE.h"

int main(int argc, char ** argv) {

	int nproma = 600;
	int nlevs = 64;
	int nblocks = 4;

    TKE_Init(nproma, nlevs, nblocks);

    return 0;

}