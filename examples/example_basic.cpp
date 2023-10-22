#include <memory>
#include "TKE.hpp"

int main(int argc, char ** argv) {
    int nproma = 1200;
    int nlevs = 64;
    int nblocks = 2;
    int block_size = nproma;
    int start_index = 0;
    int end_index = nproma - 1;
    std::shared_ptr<TKE> ocean_physics;
    ocean_physics.reset(new TKE(nproma, nlevs, nblocks, 
                                block_size, start_index, end_index));

    ocean_physics.reset();
}
